#!/bin/bash
### Pipeline to filter vcf file (only done once to prepare list SNP positions) ###
### run command: sbatch -J clean_vcf -o ${dirin}/log/vcf_cleanup.out -e ${dirin}/log/vcf_cleanup.err --wrap="${script}/vcf_cleanup/run_vcf_cleanup.sh" ###

############################################ modules #####################################################
#module load system/R-3.5.1
#module load bioinfo/bcftools-1.6
#module load bioinfo/tabix-0.2.5
#module load bioinfo/vcftools-0.1.15
#module load bioinfo/samtools-1.8
########################################################################################################## 

############################################ parameters ################################################## 
username="seynard"
SCRIPTS="/work/project/dynagen/seynard/scriptGWAS/vcf_cleanup/"
DIRIN="/genphyse/dynagen/BeeStrong/vcf_SeqApiPop/"
DIROUT="/work/project/dyangen/seynard/GWAS/BeeStrongHav3_1/"
VCFIN="MetaGenotypesCalled870.vcf.gz" #vcf avant filtres
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
rename plot_decision_ plot_decision_diag_${run} ${dir_out}/plot_decision_*.pdf
rename _value _value_diag_${run} ${dir_out}/*_value.txt
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
