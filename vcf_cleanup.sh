#! /bin/bash
###### Script to clean initial vcf file ######
#Sonia Eynard 2018 - Alain Vignal 2019 - Sonia Eynard 2020

######################
#Variables
######################
username=${1}
SCRIPTS=${2}
DIRIN=${3}
DIROUT=${4}
VCFIN=${5} #vcf before all filters
limit_allele=${6} #accept up to tri-allelic
#for all filter threshold if x=-999 the value will be assigned eliminate quantile_prob*100% of the data 
limit_FS=${7} #Phred-scaled probability that there is strand bias at the site. Often 60
limit_SOR=${8} #Strand Odds Ratio, strand bas. Often 3
limit_MQ=${9} #Root mean square mapping quality over all the reads at the site. Often 40
limit_MQRankSum=${10} #u-based z-approximation from the Rank Sum Test for mapping qualities, compares the mapping qualities of the reads supporting the reference allele and the alternate allele. Often -12.5
limit_ReadPosRankSum=${11} #u-based z-approximation from the Rank Sum Test for site position within reads, compares whether the positions of the reference and alternate alleles are different within the reads. Often -8
limit_QUAL=${12} #Quality of ALT allele. The more samples have the ATL allele, the higher the QUAL score. Often QUAL=200 
limit_QD=${13} #QUAL score normalized by allele depth (AD). Often QD=20
limit_GQ=${14} #Genotype Quality. Often GQ=10
limit_miss=${15} #percent missing data allowed for a SNP. Often Miss=0.05
limit_het=${16} #percent heterozygotes allowed for a SNP. Often Het=0.01
limit_GQfiltered=${17} #limit of missing data due to GQ allowed. Often GQ_p=0.2 (remove markers for which GQ_p% of the values are below limit_GQ)
quantile_prob_above_threshold=${18} #Quantile probability of draw diagnostic plot
quantile_prob_below_threshold=${19} #Quantile probability of draw diagnostic plot
kept_above_threshold="${20}" #list of filters for which we keep markers above threshold
kept_below_threshold="${21}" #list of filters for which we keep markers below threshold
run=${22} #type of run we are doing, can be diagnostic or filter
#define possible filters 'FS','SOR','MQ','MQRankSum','ReadPosRankSum','QUAL','QD','miss','het','GQ','GQfiltered','allele'

vcf_fs_sor_mq=${VCFIN/.vcf.gz/_fs_sor_mq.vcf} #FS_SOR_MQ
vcf_qual=${VCFIN/.vcf.gz/_qual.vcf} #quality_qd
vcf_multi=${VCFIN/.vcf.gz/_multi.vcf} #multiallele
vcf_miss=${VCFIN/.vcf.gz/_miss.vcf} #missing
vcf_het=${VCFIN/.vcf.gz/_het.vcf} #heterozygous
vcf_gq=${VCFIN/.vcf.gz/_gq.vcf} #all filters passed. No limit_GQ set to missing. SNP removed when more than percent_GQfiltered of  GQ < limit_GQ
vcf=${VCFIN/.vcf.gz/_allfilter.vcf} #all filters passed

############################################################################################
#Extract information per markers from vcf to draw diagnostic plots for vcf filtering 
############################################################################################
if [ ${run} == 'diagnostic' ] 
then 
	if [ -f ${DIROUT}/info.txt ]
	then 
		echo 'input files info already exist'
	else 
		bcftools query -f '%CHROM %POS %REF %ALT %FS %SOR %MQ %MQRankSum %ReadPosRankSum %QUAL %QD\n' ${DIRIN}/${VCFIN} > ${DIROUT}/info.txt
		head='CHROM POS REF ALT FS SOR MQ MQRankSum ReadPosRankSum QUAL QD\n'
		sed -i "1s/^/${head}/"  ${DIROUT}/info.txt
	fi
	
	if [ -f ${DIROUT}/geno.txt ]
	then 
		echo 'input files geno already exist'
		if [ -f ${DIROUT}/geno_aa.gt ]
		then 
			echo 'split geno already done'
		else
			split -l 100000 ${DIROUT}/geno.txt ${DIROUT}/geno_
			for file in ${DIROUT}/geno_*; do mv "$file" "${file}.gt"; done
		fi		
	else 
		bcftools query -f '%CHROM %POS [ %GT]\n' ${DIRIN}/${VCFIN} > ${DIROUT}/geno.txt 
		split -l 100000 ${DIROUT}/geno.txt ${DIROUT}/geno_
		for file in ${DIROUT}/geno_*; do mv "$file" "${file}.gt"; done
	fi
	
	if [ -f ${DIROUT}/gq.txt ]
	then 
		echo 'input files GQ already exist'
		if [ -f ${DIROUT}/gq_aa.gq ]
		then 
			echo 'split GQ already done'
		else
			split -l 100000 ${DIROUT}/gq.txt ${DIROUT}/gq_
			for file in ${DIROUT}/gq_*; do mv "$file" "${file}.gq"; done
		fi
		
	else 
		bcftools query -f '%CHROM %POS [ %GQ]\n' ${DIRIN}/${VCFIN}> ${DIROUT}/gq.txt
		split -l 100000 ${DIROUT}/gq.txt ${DIROUT}/gq_
		for file in ${DIROUT}/gq_*; do mv "$file" "${file}.gq"; done
	fi
	wc -l ${DIROUT}/*.txt

	IFS='~'
	read -ra kept_above <<< "${kept_above_threshold}"
	read -ra kept_below <<< "${kept_below_threshold}"
	l_array=$(echo ${#kept_above[@]})
	l_array=$((${l_array}-1))
	for ((i=0; i<=l_array;i++))
	do
		k_above=${kept_above[i]}
		k_below=${kept_below[i]}
		args_above=()
		IFS='_' read -r -a array <<< "${k_above}"
		for j in "${array[@]}"; do v="limit_${j}"; args_above+=${!v}; args_above+='_'; done
		args_below=()
		IFS='_' read -r -a array <<< "${k_below}"
		for j in "${array[@]}"; do v="limit_${j}"; args_below+=${!v}; args_below+='_'; done
		echo ${quantile_prob_above_threshold} 
		echo ${k_above} 
		echo ${args_above}
		echo ${quantile_prob_below_threshold} 
		echo ${k_below} 
		echo ${args_below}
		sbatch --mem=100G -W -J diagnostic_plot_${k_above}_${k_below} -o ${DIROUT}/log/diagnostic_plot_${k_above}_${k_below}.o -e ${DIROUT}/log/diagnostic_plot_${k_above}_${k_below}.e \
				--wrap="Rscript ${SCRIPTS}/combine_filter.r ${DIROUT} ${quantile_prob_above_threshold} ${k_above} ${args_above} ${quantile_prob_below_threshold} ${k_below} ${args_below}"
	done
else
rm ${DIROUT}/*.gt ${DIROUT}/*.gq 
############################################################################################
#Adjust filters 
############################################################################################
	if [ ${limit_FS} = -999 ]; then limit_FS=$(more ${DIROUT}/FS_value.txt); fi
	if [ ${limit_SOR} = -999 ]; then limit_SOR=$(more ${DIROUT}/SOR_value.txt); fi
	if [ ${limit_MQ} = -999 ]; then limit_MQ=$(more ${DIROUT}/MQ_value.txt); fi
	if [ ${limit_MQRankSum} = -999 ]; then limit_MQRankSum=$(more ${DIROUT}/MQRankSum_value.txt); fi
	if [ ${limit_ReadPosRankSum} = -999 ]; then limit_ReadPosRankSum=$(more ${DIROUT}/ReadPosRankSum_value.txt); fi
	if [ ${limit_QUAL} = -999 ]; then limit_QUAL=$(more ${DIROUT}/QUAL_value.txt); fi
	if [ ${limit_QD} = -999 ]; then limit_QD=$(more ${DIROUT}/QD_value.txt); fi
	if [ ${limit_GQ} = -999 ]; then limit_GQ=$(more ${DIROUT}/GQ_value.txt); fi
	if [ ${limit_miss} = -999 ]; then	limit_miss=$(more ${DIROUT}/miss_value.txt); fi
	if [ ${limit_het} = -999 ]; then	limit_het=$(more ${DIROUT}/het_value.txt); fi
	if [ ${limit_GQfiltered} = -999 ]; then limit_GQfiltered=$(more ${DIROUT}/GQfiltered_value.txt); fi
############################################################################################
#Run filters
############################################################################################
	if [ ${run} == 'filter_sequencial' ]
	then
		### option 1 group by group
		#filter on FS, SOR and MQ
		bcftools view -i 'FS<'${limit_FS}' & SOR<'${limit_SOR}' & MQ>'${limit_MQ}'' ${DIRIN}/${VCFIN} > ${DIROUT}/${vcf_fs_sor_mq}
		awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}}{ print $(f["CHROM"]), $(f["POS"]), $(f["FS"]), $(f["SOR"]), $(f["MQ"]) }' ${DIROUT}/info.txt > ${DIROUT}/tmp
		awk -F' ' '$3>='${limit_FS}' || $4>='${limit_SOR}' || $5<='${limit_MQ}'' ${DIROUT}/tmp > ${DIROUT}/fs_sor_mq_removed.txt 
		#filter on QUAL and QD
		bcftools view -i 'QUAL>'${limit_QUAL}' & QD>'${limit_QD}'' ${DIROUT}/${vcf_fs_sor_mq} > ${DIROUT}/${vcf_qual}
		awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}}{ print $(f["CHROM"]), $(f["POS"]), $(f["QUAL"]), $(f["QD"]) }' ${DIROUT}/info.txt > ${DIROUT}/tmp
		awk -F' ' '$3<='${limit_QUAL}' || $4<='${limit_QD}'' ${DIROUT}/tmp > ${DIROUT}/qual_removed.txt
		#filter on number of alleles
		vcftools --vcf ${DIROUT}/${vcf_qual} --positions ${DIROUT}/snp_kept_nb_allele.txt --recode --out ${DIROUT}/${vcf_multi}
		mv ${DIROUT}/${vcf_multi}.recode.vcf ${DIROUT}/${vcf_multi}
		#filter on missing
		bcftools filter -e'F_MISSING < '${limit_miss}'' ${DIROUT}/${vcf_multi} > ${DIROUT}/${vcf_miss}
		awk -F' ' ''${limit_miss}'<=$3' ${DIROUT}/missing.txt > ${DIROUT}/missing_removed.txt 
		#filter on heterozygous
		bcftools +setGT ${DIROUT}/${vcf_miss} -- -t q -i 'GT="het"' -n "./." > ${DIROUT}/tmp.vcf
		bcftools filter -e'F_MISSING < '${limit_het}'' ${DIROUT}/tmp.vcf > ${DIROUT}/${vcf_het}
		awk -F' ' ''${limit_het}'<=$3' ${DIROUT}/heterozygous.txt > ${DIROUT}/heterozygous_removed.txt 
		#filter on GQ 
		vcftools --vcf ${DIROUT}/${vcf_het} --minGQ 10 --recode --out ${DIROUT}/tmp.vcf  
		bcftools filter -e'F_MISSING < '${limit_GQfiltered}'' ${DIROUT}/tmp.vcf.recode.vcf > ${DIROUT}/${vcf_gq}
		awk -F' ' ''${limit_GQfiltered}'<=$3' ${DIROUT}/GQfiltered_below*.txt > ${DIROUT}/gq_removed.txt 
		#final
		cp ${DIROUT}/${vcf_gq} ${DIROUT}/${vcf}
		rm ${DIROUT}/tmp.vcf.recode.vcf ${DIROUT}/tmp.vcf.log ${DIROUT}/tmp ${DIROUT}/tmp.vcf ${DIROUT}/MetaGenotypesCalled870_snps_OK_NoIndel_multi.vcf.log

		#Counting SNPs on the chromosomes in the various vcf files
		zgrep -v '#' ${DIRIN}/${VCFIN} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfInitial
		grep -v '#' ${DIROUT}/${vcf_fs_sor_mq} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfAllFilteredSNPs_fs_sor_mq
		grep -v '#' ${DIROUT}/${vcf_qual} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfAllFilteredSNPs_qual
		grep -v '#' ${DIROUT}/${vcf_multi} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfAllFilteredSNPs_multi
		grep -v '#' ${DIROUT}/${vcf_miss} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfAllFilteredSNPs_miss
		grep -v '#' ${DIROUT}/${vcf_het} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfAllFilteredSNPs_het
		grep -v '#' ${DIROUT}/${vcf_gq} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfAllFilteredSNPsFinal_GqProportion
		grep -v '#' ${DIROUT}/${vcf} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfAllFilteredSNPsFinal
		sbatch -W -J count_phased_geno -o ${DIROUT}/log/count_phased_geno.o -e ${DIROUT}/log/count_phased_geno.e \
			--wrap="python ${SCRIPTS}/count_phased_geno.py ${DIROUT} geno.txt"

		#Compressing vcfs
		bgzip -c ${DIROUT}/${vcf} > ${DIROUT}/${vcf}.gz
		bgzip -c ${DIROUT}/${vcf_fs_sor_mq} > ${DIROUT}/${vcf_fs_sor_mq}.gz
		bgzip -c ${DIROUT}/${vcf_qual} > ${DIROUT}/${vcf_qual}.gz
		bgzip -c ${DIROUT}/${vcf_multi} > ${DIROUT}/${vcf_multi}.gz
		bgzip -c ${DIROUT}/${vcf_miss} > ${DIROUT}/${vcf_miss}.gz
		bgzip -c ${DIROUT}/${vcf_het} > ${DIROUT}/${vcf_het}.gz
		bgzip -c ${DIROUT}/${vcf_gq} > ${DIROUT}/${vcf_gq}.gz

	elif [ ${run} == 'filter_all' ]
	then
		#### option 2 all filter at once
		#list markers pass filters
		kept_above_threshold="$(sed s/~/_/g <<<$kept_above_threshold)"
		kept_below_threshold="$(sed s/~/_/g <<<$kept_below_threshold)"
		IFS='_' read -r -a array <<< "${kept_above_threshold}"
		for j in "${array[@]}"; do v="limit_${j}"; args_above+=${!v}; args_above+='_'; done
		args_below=()
		IFS='_' read -r -a array <<< "${kept_below_threshold}"
		for j in "${array[@]}"; do v="limit_${j}"; args_below+=${!v}; args_below+='_'; done
		echo ${kept_above_threshold} 
		echo ${args_above}
		echo ${kept_below_threshold} 
		echo ${args_below}
		sbatch --mem=20G -W -J filter_${k_above}_${k_below} -o ${DIROUT}/log/filter_${k_above}_${k_below}.o -e ${DIROUT}/log/filter_${k_above}_${k_below}.e \
			--wrap="Rscript ${SCRIPTS}/filter_venn_diag.r ${DIROUT} ${kept_above_threshold} ${args_above} ${kept_below_threshold} ${args_below}"

		#filter
		vcftools --gzvcf ${DIRIN}/${VCFIN} --positions ${DIROUT}/list_kept_test.txt --recode --out ${DIROUT}/${vcf}
		mv ${DIROUT}/${vcf}.recode.vcf ${DIROUT}/${vcf}

		#Counting SNPs on the chromosomes in the various vcf files
		zgrep -v '#' ${DIRIN}/${VCFIN} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfInitial
		grep -v '#' ${DIROUT}/${vcf} | cut -f 1 | sort | uniq -c | awk 'BEGIN {OFS="\t";sum=0}{print $2, $1; sum += $1} END {print "Sum", sum}' > ${DIROUT}/countVcfAllFilteredSNPsFinal
		sbatch -W -J count_phased_geno -o ${DIROUT}/log/count_phased_geno.o -e ${DIROUT}/log/count_phased_geno.e \
			--wrap="python ${SCRIPTS}/count_phased_geno.py ${DIROUT} geno.txt"

		#Compressing vcfs
		bgzip -c ${DIROUT}/${vcf} > ${DIROUT}/${vcf}.gz

	fi
fi




