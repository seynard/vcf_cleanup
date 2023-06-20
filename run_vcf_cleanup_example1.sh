#!/bin/bash
### Pipeline to filter vcf file (only done once to prepare list SNP positions) ###
### run command: sbatch -J clean_vcf -o ${dirin}/log/vcf_cleanup.out -e ${dirin}/log/vcf_cleanup.err --wrap="${script}/vcf_cleanup/run_vcf_cleanup.sh" ###

############################################ parameters ################################################## 
username="seynard"
SCRIPTS="bee_small_sample/vcf_cleanup/"
DIRIN="bee_small_sample/"
DIROUT="bee_small_sample/"
VCFIN="corsica_chr16_3_4_Mb.vcf.gz" #vcf avant filtres
limit_allele=3 #accept up to tri-allelic
type="haploid"
#for all filter threshold if x=-999 the value will be assigned eliminate quantile_prob*100% of the data
########################################################################################################## 

mkdir -p ${DIROUT}/log

########################################################################################################## 
# OPTION 1, no idea on how to select thresholds, automatic choice of values 
##########################################################################################################
########################################################################################################## 
# run diagnostic plots
# adjust values 
##########################################################################################################
limit_FS=-999 
limit_SOR=-999 
limit_MQ=-999 
limit_MQRankSum=-999 
limit_ReadPosRankSum=-999 
limit_QUAL=-999 
limit_QD=-999
limit_GQ=-999 
limit_miss=-999 
limit_het=-999 
limit_GQfiltered=-999 
quantile_prob_above_threshold=0.1 
quantile_prob_below_threshold=0.9 
kept_above_threshold="MQ_MQRankSum_ReadPosRankSum_QUAL_QD~GQ~GQ" 
kept_below_threshold="FS_SOR_allele~miss_het~GQfiltered" 
run='diagnostic' 
########################################################################################################## 
sbatch -W -J vcf_diag_${run} -o ${DIROUT}/log/vcf_diag_${run}.o -e ${DIROUT}/log/vcf_diag_${run}.e \
		--wrap="${SCRIPTS}/vcf_cleanup.sh ${username} ${SCRIPTS} ${DIRIN} ${DIROUT} ${VCFIN} \
				${limit_allele} ${limit_FS} ${limit_SOR} ${limit_MQ} ${limit_MQRankSum} ${limit_ReadPosRankSum} \
				${limit_QUAL} ${limit_QD} ${limit_GQ} ${limit_miss} ${limit_het} ${limit_GQfiltered} \
				${quantile_prob_above_threshold} ${quantile_prob_below_threshold} ${kept_above_threshold} \
				${kept_below_threshold} ${run} ${type}" 
########################################################################################################## 
# run filters
# adjust values
########################################################################################################## 
limit_FS=$(cat ${DIROUT}/FS_value.txt)
limit_SOR=$(cat ${DIROUT}/SOR_value.txt)
limit_MQ=$(cat ${DIROUT}/MQ_value.txt) 
limit_MQRankSum=-12.5
limit_ReadPosRankSum=-8
limit_QUAL=$(cat ${DIROUT}/QUAL_value.txt)
limit_QD=$(cat ${DIROUT}/QD_value.txt)
limit_GQ=$(cat ${DIROUT}/GQ_value.txt) 
limit_miss=$(cat ${DIROUT}/miss_value.txt)
limit_het=$(cat ${DIROUT}/het_value.txt)
limit_GQfiltered=$(cat ${DIROUT}/GQfiltered_value.txt)
kept_above_threshold="MQ_QUAL_QD~GQ~GQ" 
kept_below_threshold="FS_SOR_allele~miss_het~GQfiltered" 
quantile_prob_above_threshold=0.1 
quantile_prob_below_threshold=0.9 
run='filter_all' 
########################################################################################################## 
sbatch -W -J vcf_cleanup_${run} -o ${DIROUT}/log/vcf_cleanup_${run}.o -e ${DIROUT}/log/vcf_cleanup_${run}.e \
		--wrap="${SCRIPTS}/vcf_cleanup.sh ${username} ${SCRIPTS} ${DIRIN} ${DIROUT} ${VCFIN} \
				${limit_allele} ${limit_FS} ${limit_SOR} ${limit_MQ} ${limit_MQRankSum} ${limit_ReadPosRankSum} \
				${limit_QUAL} ${limit_QD} ${limit_GQ} ${limit_miss} ${limit_het} ${limit_GQfiltered} \
				${quantile_prob_above_threshold} ${quantile_prob_below_threshold} ${kept_above_threshold} \
				${kept_below_threshold} ${run} ${type}" 

rm ${DIROUT}/*.gt ${DIROUT}/*.gq ${DIROUT}/geno*.txt
##########################################################################################################
