/*  GATK4 Variant Calling Pipeline 
 *  Usage: nextflow run gencorefacility/variant-calling-pipeline-gatk4 -with-docker gencorefacility/variant-calling-pipeline-gatk4 
 *
 *  Author: Mohammed Khalfan < mkhalfan@nyu.edu >
 *  NYU Center for Genetics and System Biology 2020
 */

// Setting some defaults here,
// can be overridden in config or via command line
params.out = "${params.outdir}/out"
params.tmpdir = "${params.outdir}/gatk_temp"
params.snpeff_data = "${params.outdir}/snpeff_data"

// Print some stuff here
println "reads: $params.reads"
println "ref: $params.ref"
println "output: $params.out"
println "gatk temp dir: $params.tmpdir"
println "snpeff db: $params.snpeff_db"
println "snpeff data: $params.snpeff_data"

// Setup the reference file
ref = file(params.ref)

/* Prepare the fastq read pairs for input.
 * While doing this, count number of input samples
 */
num_samples = 0
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .tap { read_pairs_ch }
    .subscribe({ num_samples += 1 })

process align {
    publishDir "${params.out}/aligned_reads", mode:'copy'
	
    input:
    set pair_id, file(reads) from read_pairs_ch
     
    output:
    set val(pair_id), file("${pair_id}_aligned_reads.sam") \
	into aligned_reads_ch
	
    script:
    readGroup = \
	"@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${pair_id}"
    """
    bwa mem \
	-K 100000000 \
	-v 3 \
	-t ${task.cpus} \
	-Y \
	-R \"${readGroup}\" \
	$ref \
	${reads[0]} \
	${reads[1]} \
	> ${pair_id}_aligned_reads.sam
    """
}

process markDuplicatesSpark {
    publishDir "${params.out}/dedup_sorted", mode:'copy'

    input:
    set val(pair_id), file(aligned_reads) from aligned_reads_ch

    // If we're doing this step, it's the first round
    // so we set val(1) (round = 1).
    output:
    set val(pair_id), \
	val(1), \
	file("${pair_id}_sorted_dedup.bam") \
	into bam_for_variant_calling, \
	sorted_dedup_ch_for_metrics, \
	bam_for_bqsr
    set val(pair_id), \
	file ("${pair_id}_dedup_metrics.txt") \
	into dedup_qc_ch

    script:
    """
    mkdir -p ${params.tmpdir}/${workflow.runName}/${pair_id}
    gatk --java-options "-Djava.io.tmpdir=${params.tmpdir}/${workflow.runName}/${pair_id}" \
	 MarkDuplicatesSpark \
	-I $aligned_reads \
	-M ${pair_id}_dedup_metrics.txt \
	-O ${pair_id}_sorted_dedup.bam 
    rm -r ${params.tmpdir}/${workflow.runName}/${pair_id}
    """ 
}

process getMetrics {
    publishDir "${params.out}/metrics", mode:'copy'

    input:
    set val(pair_id), \
	val(round), \
	file(sorted_dedup_reads) \
	from sorted_dedup_ch_for_metrics

    output:
    set val(pair_id), 
	file("${pair_id}_alignment_metrics.txt"), \
	file("${pair_id}_insert_metrics.txt"), \
	file("${pair_id}_insert_size_histogram.pdf"), \
	file("${pair_id}_depth_out.txt") \
	into metrics_qc_ch

    script:
    """
    java -jar \$PICARD_JAR \
        CollectAlignmentSummaryMetrics \
	R=${params.ref} \
        I=${sorted_dedup_reads} \
	O=${pair_id}_alignment_metrics.txt
    java -jar \$PICARD_JAR \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	OUTPUT=${pair_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${pair_id}_insert_size_histogram.pdf 
    samtools depth -a ${sorted_dedup_reads} > ${pair_id}_depth_out.txt
    """
}

/* Run HaplotypeCaller on the initial clean bam, and the recalibrated bam
 * which we '.mix' in. This channel is automatically closed after
 * we .'take' num_samples number of objects from it. If we don't close
 * it (using .take for example), it will stay open.
 * Have to create the channel before we can mix it in (here), then we 
 * output to this channel in the bqsr process.
 */
recalibrated_bam_ch = Channel.create()
process haplotypeCaller {
    input:
    set val(pair_id), \
	val(round), \
	file(input_bam) \
	from bam_for_variant_calling.mix(recalibrated_bam_ch.take(num_samples))

    output:
    set val(pair_id), val(round), \
	file("${pair_id}_raw_variants_${round}.vcf") \
	into hc_output_ch

    script:
    """
    gatk HaplotypeCaller \
	-R $ref \
	-I $input_bam \
	-O ${pair_id}_raw_variants_${round}.vcf \
    """
}

process selectVariants {
    input:
    set val(pair_id), \
	val(round), \
	file(raw_variants) \
	from hc_output_ch

    output:
    set val(pair_id), \
	val(round), \
	file("${pair_id}_raw_snps_${round}.vcf") \
	into raw_snps_ch
    set val(pair_id), \
	val(round), \
	file("${pair_id}_raw_indels_${round}.vcf") \
	into raw_indels_ch

    script:
    """
    gatk SelectVariants \
	-R $ref \
	-V $raw_variants \
	-select-type SNP \
	-O ${pair_id}_raw_snps_${round}.vcf
    gatk SelectVariants \
        -R $ref \
        -V $raw_variants \
        -select-type INDEL \
        -O ${pair_id}_raw_indels_${round}.vcf
    """
}

process filterSnps {
    publishDir "${params.out}/filtered_snps", mode:'copy'
    
    input:
    set val(pair_id),
	val(round),
	file(raw_snps) from raw_snps_ch

    output:
    set val(pair_id), 
	val(round),
	file("${pair_id}_filtered_snps_${round}.vcf"),
	file("${pair_id}_filtered_snps_${round}.vcf.idx") \
	into filtered_snps_ch_1, filtered_snps_ch_2

    script:
    """
    gatk VariantFiltration \
	-R $ref \
	-V $raw_snps \
	-O ${pair_id}_filtered_snps_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    """
}

process filterIndels {
    publishDir "${params.out}/filtered_indels", mode:'copy'
    
    input:
    set val(pair_id), \
	val(round), \
	file(raw_indels) \
	from raw_indels_ch

    output:
    set val(pair_id), \
	file("${pair_id}_filtered_indels_${round}.vcf"), \
	file("${pair_id}_filtered_indels_${round}.vcf.idx") \
	into filtered_indels_for_recal

    script:
    """
    gatk VariantFiltration \
        -R $ref \
        -V $raw_indels \
        -O ${pair_id}_filtered_indels_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}

/* if round 1 (it[1] == 1), send snps to BQSR for recal, and snps to qc
 * if round 2 (it[1] == 2), send snps to snpeff and qc
 * todo: change this to use the branch operator
 */
filtered_snps_ch_1.filter({it[1] == 1}).tap{filtered_snps_for_recal}.tap{snps_1_qc_ch}
filtered_snps_ch_2.filter({it[1] == 2}).tap{snps_2_qc_ch}.tap{filtered_snps_for_snpeff}

process bqsr{
    publishDir "${params.out}/bqsr", mode:'copy'

    input:
    set val(pair_id), \
	val(round), \
	file(input_bam), \
	val(round), \
	file(filtered_snps), \
	file(filtered_snps_index), \
	file(filtered_indels), \
	file(filtered_indels_index) \
	from bam_for_bqsr
	.join(filtered_snps_for_recal)
	.join(filtered_indels_for_recal)
    
    output:    
    set val(pair_id), \
	file("${pair_id}_recal_data.table"), \
	file("${pair_id}_post_recal_data.table") \
	into analyze_covariates_in_ch
    set val(pair_id), \
	val(new_round), \
	file("${pair_id}_recal.bam") \
	into recalibrated_bam_ch

    // here is where we iterate the round. 
    // keep this dynamic to allow for doing 
    // multiple rounds of bqsr.
    // note: run SelectVariants to exclude 
    // filtered variants from vcf before bqsr!
    // todo: this step can be optimized potentially by 
    // breaking out SelectVariants, BaseRecalibrator,
    // and ApplyBQSR into individual processes 
    script:
    new_round=round + 1
    """
    echo "New Round: " $new_round
    gatk SelectVariants \
	--exclude-filtered \
	-V $filtered_snps \
	-O ${pair_id}_bqsr_snps.vcf
    gatk SelectVariants \
        --exclude-filtered \
        -V $filtered_indels \
        -O ${pair_id}_bqsr_indels.vcf
    gatk BaseRecalibrator \
	-R $ref \
	-I $input_bam \
	--known-sites ${pair_id}_bqsr_snps.vcf \
	--known-sites ${pair_id}_bqsr_indels.vcf \
	-O ${pair_id}_recal_data.table
    gatk ApplyBQSR \
        -R $ref \
        -I $input_bam \
        -bqsr ${pair_id}_recal_data.table \
        -O ${pair_id}_recal.bam
    gatk BaseRecalibrator \
        -R $ref \
	-I ${pair_id}_recal.bam \
        --known-sites ${pair_id}_bqsr_snps.vcf \
	--known-sites ${pair_id}_bqsr_indels.vcf \
	-O ${pair_id}_post_recal_data.table
    """	
}

process analyzeCovariates{
    publishDir "${params.out}/bqsr", mode:'copy'

    input:
    set val(pair_id), file(recal_table), file(post_recal_table) \
	 from analyze_covariates_in_ch

    output:
    set val(pair_id), file("${pair_id}_recalibration_plots.pdf") \
	into analyzed_covariates_ch

    script:
    """
    gatk AnalyzeCovariates \
	-before $recal_table \
	-after $post_recal_table \
	-plots ${pair_id}_recalibration_plots.pdf
    """
}

process snpEff {
    publishDir "${params.out}/snpeff", mode:'copy'    

    input:
    set val(pair_id), \
	val(round), \
	file(filtered_snps), \
	file(filtered_snps_index) \
	from filtered_snps_for_snpeff

    output:
    file '*' into snpeff_out

    script:
    """
    java -jar \$SNPEFF_JAR -v \
	-dataDir $params.snpeff_data \
	$params.snpeff_db \
	$filtered_snps > ${pair_id}_filtered_snps.ann.vcf
    """
}

process qc {
    input:
    set val(pair_id), \
	file("x${pair_id}_dedup_metrics.txt"), \
	file("${pair_id}_alignment_metrics.txt"), \
        file("${pair_id}_insert_metrics.txt"), \
        file("${pair_id}_insert_size_histogram.pdf"), \
        file("${pair_id}_depth_out.txt"), \
	val(round_1), \
        file("${pair_id}_filtered_snps_1.vcf"), \
        file("${pair_id}_filtered_snps_1.vcf.idx"), \
	val(round_2), \
        file("${pair_id}_filtered_snps_2.vcf"), \
        file("${pair_id}_filtered_snps_2.vcf.idx") \
	from dedup_qc_ch
	.join(metrics_qc_ch)
        .join(snps_1_qc_ch)
        .join(snps_2_qc_ch)

    output:
    file ("${pair_id}_report.csv") into qc_output

    script:
    """
    parse_metrics.sh ${pair_id} > ${pair_id}_report.csv
    """
}

/* Process qc above creates a report for each sample.
 * Below we compile these into a single report.
 */ 
qc_output.collectFile(name: "${workflow.runName}_report.csv", keepHeader: true, storeDir: "${params.out}/reports")

