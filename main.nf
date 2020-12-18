/*  GATK4 Variant Calling Pipeline
 *  Usage: nextflow run /path/to/main.nf
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
 * Use the size parameter to not auto-group, and instead
 * use getBaseName() and remove two regexs to get the ID.
 * This is custom for NYU CGSB sequence data file naming format
 * While doing this, count number of input samples
 */
num_samples = 0
Channel
    .fromFilePairs( params.reads, size: -1)
    { file -> file.getBaseName() - ~/${params.fcid}_/ - ~/n0[12]_/ - ~/.fastq/ }
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
    set val(pair_id),
        file("${pair_id}_sorted_dedup.bam"),
        file("${pair_id}_sorted_dedup.bam.bai") \
        into full_bam_bw_ch, downsample_bam_ch, qualimap_ch

    script:
    """
    gatk \
        MarkDuplicatesSpark \
        -I $aligned_reads \
        -M ${pair_id}_dedup_metrics.txt \
        -O ${pair_id}_sorted_dedup.bam
    """
}

process downsample_bam{
    publishDir "${params.out}/downsampled_bam", mode:'copy'

    input:
    set val(pair_id),
        file(bam),
        file(bam_index) \
        from downsample_bam_ch

    output:
    set file("${pair_id}_downsampled.bam"),
        file("${pair_id}_downsampled.bam.csi") into jbrowse_bam_ch

    script:
    """
    java -jar \$SORTSAMREFNAME_JAR \
        --samoutputformat BAM \
        $bam |\
        java -jar \$BIOSTAR_JAR \
        -n 75 \
        --samoutputformat BAM |\
        samtools sort -o ${pair_id}_downsampled.bam
    samtools index -c -m 6 ${pair_id}_downsampled.bam
    """

}

process qualimap{
    input:
    set val(pair_id),
        file(bam),
        file(bam_index) from qualimap_ch

    output:
    file('*') into multiqc_qualimap_ch

    script:
    """
    qualimap BamQC -bam $bam -outdir ${pair_id} -outformat HTML
    """
}

process getMetrics{
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
        into metrics_qc_ch, metrics_multiqc_ch

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
    set val(pair_id),
        val(round),
        file(input_bam) \
        from bam_for_variant_calling.mix(recalibrated_bam_ch.take(num_samples))

    output:
    set val(pair_id), val(round),
        file("${pair_id}_raw_variants_${round}.vcf") \
        into hc_output_ch
    set val(hc_bamout_pair_id),
        file("${pair_id}_hc_bamout_${round}.bam"),
        file("${pair_id}_hc_bamout_${round}.bai") \
        into hc_bam_bw_ch

    script:
    hc_bamout_pair_id = pair_id + "-hc_bamout"
    """
    gatk HaplotypeCaller \
        -R $ref \
        -I $input_bam \
        -bamout ${pair_id}_hc_bamout_${round}.bam \
        -O ${pair_id}_raw_variants_${round}.vcf
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
    set val(pair_id),
        val(round),
        file("${pair_id}_filtered_indels_${round}.vcf"),
        file("${pair_id}_filtered_indels_${round}.vcf.idx") \
        into filtered_indels_ch_1
    file("${pair_id}_filtered_indels_${round}.vcf") \
        into filtered_indels_bzip_tabix_vcf_ch

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

/* if round 1 (it[1] == 1), send snps and indels to BQSR for recal, and snps to qc
 * if round 2 (it[1] == 2), send snps to snpeff and qc
 * todo: change this to use the branch operator
 */
filtered_snps_ch_1.filter({it[1] == 1}).tap{filtered_snps_for_recal}.tap{snps_1_qc_ch}
filtered_snps_ch_2.filter({it[1] == 2}).tap{snps_2_qc_ch}.tap{filtered_snps_for_snpeff}.tap{bcftools_stats_ch}
filtered_indels_ch_1.filter({it[1] == 1}).tap{filtered_indels_for_recal}

process bqsr{
    publishDir "${params.out}/bqsr", mode:'copy'

    input:
    set val(pair_id),
        val(round),
        file(input_bam),
        val(round),
        file(filtered_snps),
        file(filtered_snps_index),
        val(round),
        file(filtered_indels),
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

process make_bw{
    publishDir "${params.out}/bigwig", mode:'copy'

    input:
    set val(id),
        file(bam),
        file(bam_index) \
        from full_bam_bw_ch
        //.mix(hc_bam_bw_ch)

    output:
    file("${id}_coverage.bam.bw") into bw_out_ch

    // Skip haplotypecaller bamout when round 1
    when:
    bam.getName() != "${id}_hc_bamout_1.bam"

    script:
    """
    bamCoverage \
        -p max  \
        --bam $bam \
        -o ${id}_coverage.bam.bw \
        --binSize 1 \
        --ignoreDuplicates \
        --minMappingQuality 20
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
    file ("${pair_id}_filtered_snps.ann.vcf") into snpeff_bzip_tabix_vcf
    file ("${pair_id}.csv") into multiqc_snpeff_csv_ch

    script:
    """
    java -jar \$SNPEFF_JAR -v \
        -dataDir $params.snpeff_data \
        $params.snpeff_db \
        -csvStats ${pair_id}.csv \
        $filtered_snps > ${pair_id}_filtered_snps.ann.vcf

    """
}

process bzip_tabix_vcf{
    input:
    file(vcf) from filtered_indels_bzip_tabix_vcf_ch
        .mix(snpeff_bzip_tabix_vcf)

    output:
    file("*.vcf.gz*") into jbrowse_vcf_ch

    // ignore indels vcf when it's round 1
    when:
    !vcf.getName().endsWith("indels_1.vcf")

    script:
    """
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    """
}

process qc {
    input:
    set val(pair_id), \
        file("${pair_id}_dedup_metrics.txt"), \
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
    val(pair_id) into pair_id_ch

    script:
    """
    parse_metrics.sh ${pair_id} > ${pair_id}_report.csv
    """
}

process bcftools_stats{
    input:
    set val(pair_id),
        val(round),
        file(vcf),
        file(vcf_index) from bcftools_stats_ch

    output:
    file ("${pair_id}_vcf_stats.vchk") into multiqc_bcftools_stats_ch

    script:
    """
    bcftools stats $vcf > ${pair_id}_vcf_stats.vchk
    """
}

process multiqc{
    publishDir "${params.out}/reports", mode:'copy'

    input:
    file(snpeff_csv) from multiqc_snpeff_csv_ch.collect()
    file(bcftools_stats) from multiqc_bcftools_stats_ch.collect()
    file('*') from multiqc_qualimap_ch.collect()
    file('*') from metrics_multiqc_ch.collect()

    output:
    file("*.html") into multiqc_out

    script:
    """
    for f in \$(ls *metrics.txt);do sed -i 's/_sorted_dedup//g' \$f; done
    for f in \$(ls */genome_results.txt);do sed -i 's/_sorted_dedup//g' \$f; done
    for f in \$(ls *.vchk);do sed -i 's/_filtered_snps_2//g' \$f; done
    for f in \$(ls *_snpeff_stats.csv);do sed -i 's/_snpeff_stats//g' \$f; done

    # make the config file for multiqc
    echo "remove_sections:" > multiqc_conf.yml
    echo "    - bcftools-stats_indel_plot" >> multiqc_conf.yml
    echo "    - qualimap-insert-size-histogram" >> multiqc_conf.yml
    echo "    - bcftools-stats_indel_plot" >> multiqc_conf.yml

    multiqc -f -c multiqc_conf.yml --ignore *.run .
    sed -i 's/>Bcftools</>VCF Stats</g' multiqc_report.html
    """
}

/* Process qc above creates a report for each sample.
 * Below we compile these into a single report.
 */
qc_output.collectFile(name: "${workflow.runName}_report.csv", keepHeader: true, storeDir: "${params.out}/reports")

/* Collect all the pair_ids and send them to pair_id_list_ch for use in jbrowse.
 */
pair_id_ch.collectFile(storeDir: "${params.out}/reports") { item ->
       [ "sample_ids.txt", item + '\n' ]
}.tap{pair_id_list_ch}

process jbrowse{

    input:
    file '*' from jbrowse_vcf_ch.collect()
    file '*' from jbrowse_bam_ch.collect()
    file '*' from bw_out_ch.collect()
    file(pair_id_list) from pair_id_list_ch

    when:
    params.do_jbrowse

    script:
    """
    cgsb_upload2jbrowse \
        -p $params.jbrowse_pi \
        -d $params.dataset_name \
        -s $pair_id_list \
        -f . \
        $ref \
        $params.gff
    """
}
