#! /bin/bash
### vcf filter (only done once to prepare list SNP positions) ###

############################################ modules #####################################################
module load system/R-3.5.1
module load bioinfo/bcftools-1.6
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/samtools-1.8
########################################################################################################## 

############################################ parameters ################################################## 
username=seynard
SCRIPTS='/work/project/dynagen/seynard/queen_genotype_reconstruction/code/vcfcleanup'
DIRIN='/work/project/dynagen/seynard/queen_genotype_reconstruction/raw_data/The870vcf'
DIROUT='/work/project/dynagen/seynard/queen_genotype_reconstruction/data'
VCFIN='MetaGenotypesCalled870_raw_snps.vcf.gz' #vcf avant filtres
limit_allele=3 #accept up to tri-allelic
#for all filter threshold if x=-999 the value will be assigned eliminate quantile_prob*100% of the data
limit_FS=-999 #Phred-scaled probability that there is strand bias at the site. Often 61 (<=60)
limit_SOR=-999 #Strand Odds Ratio, strand bas. Often 4 (<=3)
limit_MQ=-999 #Root mean square mapping quality over all the reads at the site. Often 39 (=>40)
limit_MQRankSum=-999 #u-based z-approximation from the Rank Sum Test for mapping qualities, compares the mapping qualities of the reads supporting the reference allele and the alternate allele. Often -12.5
limit_ReadPosRankSum=-999 #u-based z-approximation from the Rank Sum Test for site position within reads,compares whether the positions of the reference and alternate alleles are different within the reads. Often -8
limit_QUAL=-999 #Quality of ALT allele. The more samples have the ATL allele, the higher the QUAL score. Often QUAL=200 
limit_QD=-999 #QUAL score normalized by allele depth (AD). Often QD=20
limit_GQ=-999 #Genotype Quality. Often GQ=10
limit_miss=-999 #percent missing data allowed for a SNP. Often Miss=0.05
limit_het=-999 #percent heterozygotes allowed for a SNP. Often Het=0.01
limit_GQfiltered=-999 #limit of missing data due to GQ allowed. Often GQ_p=0.2 (remove markers for which GQ_p% of the values are below limit_GQ)
quantile_prob_above_threshold=0.1 #Quantile probability of draw diagnostic plot
quantile_prob_below_threshold=0.9 #Quantile probability of draw diagnostic plot
kept_above_threshold="MQ_MQRankSum_ReadPosRankSum_QUAL_QD~GQ~GQ" #list of filters for which we keep markers above threshold
kept_below_threshold="FS_SOR_allele~miss_het~GQfiltered" #list of filters for which we keep markers below threshold
run='diagnostic' #type of run we are doing, can be diagnostic, filter_all, filter_sequencial
########################################################################################################## 

mkdir -p ${DIROUT}/log

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
sbatch -W -J vcf_cleanup -o ${DIROUT}/log/vcf_cleanup.o -e ${DIROUT}/log/vcf_cleanup.e \
		--wrap="${SCRIPTS}/vcf_cleanup.sh ${username} ${SCRIPTS} ${DIRIN} ${DIROUT} ${VCFIN} \
		${limit_allele} ${limit_FS} ${limit_SOR} ${limit_MQ} ${limit_MQRankSum} ${limit_ReadPosRankSum} \
		${limit_QUAL} ${limit_QD} ${limit_GQ} ${limit_miss} ${limit_het} ${limit_GQfiltered} \
		${quantile_prob_above_threshold} ${quantile_prob_below_threshold} ${kept_above_threshold} ${kept_below_threshold} ${run}" 


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
limit_GQ=10 
limit_miss=-999 
limit_het=-999 
limit_GQfiltered=-999 
quantile_prob_above_threshold=0.1 
quantile_prob_below_threshold=0.9 
kept_above_threshold="GQ" 
kept_below_threshold="GQfiltered" 
run='diagnostic' 
########################################################################################################## 
sbatch -W -J vcf_cleanup -o ${DIROUT}/log/vcf_cleanup.o -e ${DIROUT}/log/vcf_cleanup.e \
		--wrap="${SCRIPTS}/vcf_cleanup.sh ${username} ${SCRIPTS} ${DIRIN} ${DIROUT} ${VCFIN} \
		${limit_allele} ${limit_FS} ${limit_SOR} ${limit_MQ} ${limit_MQRankSum} ${limit_ReadPosRankSum} \
		${limit_QUAL} ${limit_QD} ${limit_GQ} ${limit_miss} ${limit_het} ${limit_GQfiltered} \
		${quantile_prob_above_threshold} ${quantile_prob_below_threshold} ${kept_above_threshold} ${kept_below_threshold} ${run}" 

########################################################################################################## 
# run filters
# adjust values 
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
sbatch -W -J vcf_cleanup -o ${DIROUT}/log/vcf_cleanup.o -e ${DIROUT}/log/vcf_cleanup.e \
		--wrap="${SCRIPTS}/vcf_cleanup.sh ${username} ${SCRIPTS} ${DIRIN} ${DIROUT} ${VCFIN} \
		${limit_allele} ${limit_FS} ${limit_SOR} ${limit_MQ} ${limit_MQRankSum} ${limit_ReadPosRankSum} \
		${limit_QUAL} ${limit_QD} ${limit_GQ} ${limit_miss} ${limit_het} ${limit_GQfiltered} \
		${quantile_prob_above_threshold} ${quantile_prob_below_threshold} ${kept_above_threshold} ${kept_below_threshold} ${run}" 

