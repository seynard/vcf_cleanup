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
# OPTION 3, we know the values we want to filter for 
##########################################################################################################
########################################################################################################## 
# run diagnostic plots
# fixed values
########################################################################################################## 
limit_FS=61
limit_SOR=4 
limit_MQ=39 
limit_MQRankSum=-12.5
limit_ReadPosRankSum=-8
limit_QUAL=200
limit_QD=20
limit_GQ=10 
limit_miss=0.05
limit_het=0.01
limit_GQfiltered=0.2 
kept_above_threshold="MQ_QUAL_QD~GQ~GQ" 
kept_below_threshold="FS_SOR_allele~miss_het~GQfiltered" 
quantile_prob_above_threshold=0.1 
quantile_prob_below_threshold=0.9 
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
# fixed values 
########################################################################################################## 
limit_FS=61
limit_SOR=4 
limit_MQ=39 
limit_MQRankSum=-12.5
limit_ReadPosRankSum=-8
limit_QUAL=200
limit_QD=20
limit_GQ=10 
limit_miss=0.05
limit_het=0.01
limit_GQfiltered=0.2 
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
