#!/bin/bash
ID=$1

input=${ID}_alignment_metrics.txt
while read line
do
	if [[ $line == PAIR* ]];then
		ALIGNMENT_METRICS=$(echo $line | cut -d' ' -f2,6,7,10,16,18 | tr ' ' ',')
	fi
done < $input

# Match the first line that starts with a number (the stats)
re='^[0-9]+([.][0-9]+)?$'
input=${ID}_insert_metrics.txt
while read line
do
	MEAN_INSERT_SIZE=$(echo $line | cut -d' ' -f6)
	if [[ $line != '#'* && $MEAN_INSERT_SIZE =~ $re ]];then
        	break
	fi
done < $input

input=${ID}_dedup_metrics.txt
while read line
do
        DUP_PCT=$(echo $line | cut -d' ' -f9)
        if [[ $line != '#'* && $DUP_PCT =~ $re ]];then
                break
        fi
done < $input


snps_1=$(grep -v '^#' ${ID}_raw_snps.vcf | wc -l)
snps_2=$(grep 'PASS' ${ID}_filtered_snps.vcf | wc -l)
snps_3=$(grep -v '^#' ${ID}_raw_snps_recal.vcf | wc -l)
snps_4=$(grep 'PASS' ${ID}_filtered_snps_final.vcf | wc -l)
avg_coverage=$(awk '{sum+=$3} END { print sum/NR}' ${ID}_depth_out.txt)

echo "ID,# reads,aligned reads,% aligned,aligned bases,read length,% paired, %dup, mean insert size,# SNPs pre bqsr,# SNPs filtered pre bqsr,# SNPs post bqsr,# SNPs filtered post bqsr,average coverage"
echo "$ID,$ALIGNMENT_METRICS,$DUP_PCT,$MEAN_INSERT_SIZE,$snps_1,$snps_2,$snps_3,$snps_4,$avg_coverage"

