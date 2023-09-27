# Filter vcf file

<font color='red'>LATEST RELEASED VERSION IS V2. CHECK THE CORRECT BRANCH</font>

Histoiry of the pipeline: written in 2018 by S Eynard, updated by A Vignal in 2019. Version 0 finalised by S Eynard and A Vignal in 2021, and used for SeqApiPop. Version 1, update of version 0, addition 
of an option to take the data type (haploid or diploid) into account when doing filtering as it impacts the filtering done on heterozygotes sites. Also small update on the plotting function.

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [1. Introduction](#1-introduction)
- [2. Filters on annotations in the vcf file](#2-filters-on-annotations-in-the-vcf-file)
	- [2.1. In the INFO field: general SNP quality estimations](#21-in-the-info-field-general-snp-quality-estimations)
	- [2.2. Sample level annotations: genotype quality estimations](#22-sample-level-annotations-genotype-quality-estimations)
- [3. Other filters](#3-other-filters)
- [4. SCRIPTS for filtering](#4-scripts-for-filtering)
	- [4.1. General variables to edit in the calling script run_vcfcleanup.sh](#41-general-variables-to-edit-in-the-calling-script-runvcfcleanupsh)
	- [4.2. Variables to edit for quality filter threshold values :](#42-variables-to-edit-for-quality-filter-threshold-values-)
	- [4.3. Variables to edit for type of run](#43-variables-to-edit-for-type-of-run)
	- [4.4. Parameters used in the study](#44-parameters-used-in-the-study)
		- [4.4.1 Plotting the diagnostic histograms](#441-plotting-the-diagnostic-histograms)
			- [4.4.1.1 Mapping quality metrics: Stand Odds Ratio (SOR)](#4411-mapping-quality-metrics-stand-odds-ratio-sor)
			- [4.4.1.2 Mapping quality metrics: Fisher Strand (FS)](#4412-mapping-quality-metrics-fisher-strand-fs)
			- [4.4.1.3 Mapping quality metrics: Mapping Quality (MQ)](#4413-mapping-quality-metrics-mapping-quality-mq)
			- [4.4.1.4 Genotyping quality metrics: SNP quality (QUAL)](#4414-genotyping-quality-metrics-snp-quality-qual)
			- [4.4.1.5 Genotyping quality metrics: quality by depth (QD)](#4415-genotyping-quality-metrics-quality-by-depth-qd)
			- [4.4.1.6. Individual genotyping metrics: heterozygote calls](#4416-individual-genotyping-metrics-heterozygote-calls)
			- [4.4.1.7. Individual genotyping metrics: missing genotypes](#4417-individual-genotyping-metrics-missing-genotypes)
			- [4.4.1.8. Individual genotyping metrics: genotype quality (QD)](#4418-individual-genotyping-metrics-genotype-quality-qd)
			- [4.4.1.9 Notice on resulting missing data.](#4419-notice-on-resulting-missing-data)
		- [4.4.2 Running the filters: generate Venn diagrams and filtered vcf](#442-running-the-filters-generate-venn-diagrams-and-filtered-vcf)
- [5. results:](#5-results)
<!-- /TOC -->

## 1. Introduction

[On hard filtering variants](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)

Variants were filtered on the INFO field and on samples-level annotations of the vcf. Additionally, we removed the SNP markers having an allele noted *, as observed for indels (InDel), that cannot be managed easily in subsequent analyses.

### R packages required
The scripts will work on a slurm cluster. In theory, they will install the required R packages if missing. If you experience any issue with this step please install the required packages ('data.table','VennDiagram','reshape2','RColorBrewer','grDevices','ggplot2','viridis') before running the script.

### To edit in run_vcfcleanup_3.sh
* run='diagnostic'
Will only make the plots
* run='filter_all'
Will produce a vcf file once all filters are passed
* run='filter_sequential'
Will produce a vcf file at each filtering stage (not usually useful)


## 2. Filters on annotations in the vcf file
The general marker INFO fields (FS, SOR, MQ...) and the sample level annotations that were analysed by plotting their distribution of values in the dataset and/or used for filtering are indicated in bold in the lists 2.1 and 2.2. Other annotations (DP, AC, AF) are indicated for reference.

### 2.1. In the INFO field: general SNP quality estimations
* **FS** = FisherStrand; phred-scaled probability that there is strand bias at the site.
* **SOR** = StrandOddsRatio: another way to estimate strand bias using a test similar to the symmetric odds ratio test.
* **MQ** = RMSMappingQuality: root mean square mapping quality over all the reads at the site.

* ***MQRankSum*** = MappingQualityRankSumTest: compares the mapping qualities of the reads supporting the reference allele and the alternate allele.
* ***ReadPosRankSum*** = ReadPosRankSumTest: compares whether the positions of the reference and alternate alleles are different within the reads.

* **QUAL** = Phred-scaled quality score for the assertion made in ALT. The more samples have the ATL allele, the higher the QUAL score
* DP = In the INFO field: combined depth across samples
* **QD** = QUAL score normalized by allele depth (AD). For a single sample, the HaplotypeCaller calculates the QD by taking QUAL/AD. For multiple samples, HaplotypeCaller and GenotypeGVCFs calculate the QD by taking QUAL/AD of samples with a non hom-ref genotype call. The reason we leave out the samples with a hom-ref call is to not penalize the QUAL for the other samples with the variant call. QD is roughly speaking QUAL / Sum(Sum(AD)). 
* AC = Allele count in genotypes, for each ALT allele, in the same order as listed
* AF = allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes.

### 2.2. Sample level annotations: genotype quality estimations
* GT = Genotype
* AD = Allele Depth. Only informative reads counted. Used for calculating QD (see INFO fields). Sum of AD can be inferior to DP.
* DP = depth for the given sample.
* PGT = Phased Genotype.
* PID = The PID contains the first site in the phased sites. For example, if sites 1,2,3 are phased, the PIDs for all the 3 sites will contain 1.
* **GQ** = Genotyping Quality: difference between the second lowest PL and the lowest PL (which is always 0.
* PL = : normalized Phred-scaled likelihoods of the genotypes considered in the variant record for each sample.
* PS : phase set. A phase set is defined as a set of phased genotypes to which this genotype belongs.

## 3. Other filters
* SNPs with more than 3 alleles were filtered out (**variable limit_allele=3**), can be ajusted to ignore tri-allelic markers directly. 
* If 'type' is set as haploid, for example when honeybee haploid drones are sequenced, SNPs with a high proportion of heterozygote calls (**variable limit_het=0.01**) were filtered out. 
The remaining heterozygote calls (< 1%) are retained and set to missing to missing as they are probably genotypiong errors. Although SNPs with more than 5% of missing data were filtered out, some markers may have more than 5% missing data after complete filtering due to the hetozygote calls set to missing. When diploid individuals are sequenced no filter on heterozygotes calls are applied. 

## 4. SCRIPTS for filtering

* [run_vcfcleanup.sh](Scripts_2_VcfCleanup/run_vcfcleanup.sh), will call the script:
  * [vcf_cleanup.sh](Scripts_2_VcfCleanup/vcf_cleanup.sh), which will call the scripts:
	* For the option 'diagnostic': 
		* [diagnostic.r](Scripts_2_VcfCleanup/diagnostic.r). Will output:
			* histogram and ecdf plots for the distribution of the various quality estimators in the input vcf.
			* Values set for filtering in run_vcfcleanup.sh will be indicated on the plots
	* For the options 'filter_sequential' or 'filter_all': 
		* [filter_list.r](Scripts_2_VcfCleanup/filter_list.r). Will output:
			* Filters will be run, producing lists of SNPs containing values of interest.
		* [filter.r](Scripts_2_VcfCleanup/filter.r). Will output:
			* Venn diagrams
			* the list of SNPs to keep: list_kept.txt.
				* list_kept.txt will be used by vcf_cleanup.sh to produce the filtered vcf.
			* the number of SNPs in the input and the output vcfs will be counted.
		* [count_phased_geno.py](Scripts_2_VcfCleanup/count_phased_geno.py). Will output:
			* count_phased_geno.txt : number of unphased, phased and missing genotype calls for each variant of the input vcf

* Options: 'diagnostic', 'filter_sequential' or 'filter_all'
'diagnostic' will produce diagnostic plots only
'filter_sequential' and 'filter_all' will perform filtering and produce filtered vcfs, 1 for each of our variable of interested one after the other if 'filter_sequential' or 1 fully filtered vcf for 'filter_all'.

* Run: option 1, 2 or 3 (run_vcf_cleanup1.sh, run_vcf_cleanup2.sh, run_vcf_cleanup3.sh)\
Option1: if there is no hypothesis on threshold values for filtering (use of -999, the diagnostic script will provide values for filtering based on quantile distribution estimation).\
Option2: if there is hypothesis on threshold values for filtering but one wants the diagnostic plots to inform on the best proposed value (use of -999 for diagnostic only. The diagnostic script will suggest values for filtering based on quantile distribution estimation).\
Option3: if there is an hypothesis on threshold values that we want to be applied for diagnostic plots and filtering of the vcf file.

### 4.1. General variables to edit in the calling script run_vcf_cleanup.sh
* All editing of paths
  - username=avignal # Deprecated
  - SCRIPTS='~/seqapipopOnHAV3_1/vcf_cleanup_scripts' #path to the other scripts called
  - DIRIN='~/seqapipopOnHAV3_1/combineGVCFs/The870vcf' #path to directory containing the input vcf
  - DIROUT='~/seqapipopOnHAV3_1/vcf_cleanup' #path to output directory
  - VCFIN='MetaGenotypesCalled870_raw_snps.vcf.gz' #name of the vcf file to filter

### 4.2. Variables to edit for quality filter threshold values :
The variables limit_FS, limit_SOR ... limit_het can either be set to a specified value, or set to -999, in which case each filter threshold will be calculated such as a percentage of the data, specified in the variables quantile_prob_above_threshold and quantile_prob_below_threshold, will be kept.

* type=haploid #can be diploid, and therefore will not apply filtering on heterozygotes calls
* limit_allele=3 #accept up to three alleles (edit to 2 or 4)
* limit_FS=61
* limit_SOR=4
* limit_MQ=39
* limit_MQRankSum=-12.5
* limit_ReadPosRankSum=-8
* limit_QUAL=200
* limit_QD=20
* limit_GQ=10
* limit_miss=0.05
* limit_het=0.01
* limit_GQfiltered=0.2
* quantile_prob_above_threshold=0.1 #In this example, any of the variables above for which variants are kept above a threshold is set to -999, will be adjusted to eliminate 10 % of the data.
* quantile_prob_below_threshold=0.9 #In this example, any of the variables above for which variants are kept below a threshold is set to -999, will be adjusted to eliminate 10 % of the data.
* kept_above_threshold="MQ_QUAL_QD\~GQ\~GQ" #variables for which we want to filter (keep SNPs) above a certain threshold
* kept_below_threshold="FS_SOR_allele\~miss_het\~GQfiltered" #variables for which we want to filter keep SNPs) below a certain threshold

For the purpose of the diagnostic, filters need to be applied in a specific order (here defined as i) technical variable such as SOR, FS, number of alleles ... ii) quality variables such as QUAL and QD and iii) sample specific variables such as rate of missing genotypes, heterozygotes calls, and GQ content. 
in:
  - kept_above_threshold="MQ_QUAL_QD ~GQ ~GQ"
  - kept_below_threshold="FS_SOR_allele ~miss_het ~GQfiltered"
Filter joined by a _ are done together. Groups of filters are separated by ~
The '~' separate groups of variables, the '_' separate each variable within a group.

### 4.3. Variables to edit for type of run
* The variable #run can take the values: 'diagnostic', 'filter_all' or 'filter_sequential'
  - diagnostic: will only output the distribution plots for each quality parameter, to help decide on threshold values setting.
  - filter_all: will filter the vcf file using all parameters set by the variables simultaneously, as specified in the kept_above_threshold and kept_below_threshold variables
  - filter_sequential: filter on each parameter set by the variables, but sequencially, and produce a dedicated vcf file for each of the applied filter. 

### 4.4. Parameters used in the study
#### 4.4.1 Plotting the diagnostic histograms
A run with run='diagnostic' will plot histograms and empirical cumulative distribution functions (ECDF), without performing the actual filtering. Once satisfactory filtering values are obtained, a second run='filter_all' will perform the filtering and produce the Venn diagrams.

```bash
#! /bin/bash

#run_vcf_cleanup.sh

#vcf filter (only done once to prepare list SNP positions)

# modules #####################################################
module load system/R-3.5.1
module load bioinfo/bcftools-1.6
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/samtools-1.8
# end modules #################################################

# parameters ##################################################
username=avignal
SCRIPTS='~/seqapipopOnHAV3_1/vcf_cleanup_scripts' #path to scripts
DIRIN='~/combineGVCFs/The870vcf' #path to directory containing the input vcf
DIROUT='~/seqapipopOnHAV3_1/vcf_cleanup' #path to output directory
VCFIN='MetaGenotypesCalled870_raw_snps.vcf.gz' #inputvcf before filters
limit_allele=3
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
quantile_prob_above_threshold=0.1
quantile_prob_below_threshold=0.9
kept_above_threshold="MQ_QUAL_QD~GQ~GQ"
kept_below_threshold="FS_SOR_allele~miss_het~GQfiltered"
type='haploid'
run='diagnostic'
# end parameters ###############################################
sbatch -W -J vcf_diag -o ${DIROUT}/log/vcf_cleanup.o -e ${DIROUT}/log/vcf_cleanup.e \
		--wrap="${SCRIPTS}/vcf_cleanup.sh ${username} ${SCRIPTS} ${DIRIN} ${DIROUT} ${VCFIN} \
		${limit_allele} ${limit_FS} ${limit_SOR} ${limit_MQ} ${limit_MQRankSum} ${limit_ReadPosRankSum} \
		${limit_QUAL} ${limit_QD} ${limit_GQ} ${limit_miss} ${limit_het} ${limit_GQfiltered} \
		${quantile_prob_above_threshold} ${quantile_prob_below_threshold} \
		${kept_above_threshold} ${kept_below_threshold} ${run} ${type}"

# end of file
```
<div style="page-break-after: always"></div>

##### 4.4.1.1 Mapping quality metrics: Stand Odds Ratio (SOR)

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_SOR_Hist_2.pdf.png)

Counts of SNPs according to SOR values. The blue dotted line indicates the threshold retained for filtering the vcf: SOR > 4.

-----------------------

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_SOR_ECDF_2.pdf.png)
Empirical cumulative distribution function of SOR values.  The blue dotted line indicates the threshold used for filtering the vcf: SOR > 4.

-----------------------

<div style="page-break-after: always"></div>

##### 4.4.1.2 Mapping quality metrics: Fisher Strand (FS)

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_FS_Hist_2.png)
Counts of SNPs according to FS values. The blue dotted line indicates the threshold used for filtering the vcf: FS < 61. X axis is on a log scale log(61)=1.785.

-----------------------

<div style="page-break-after: always"></div>

##### 4.4.1.3 Mapping quality metrics: Mapping Quality (MQ)

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_MQ_Hist_2.png)
Counts of SNPs according to MQ values. The blue dotted line indicates the threshold used for filtering the vcf: MQ > 40.

-----------------------

<div style="page-break-after: always"></div>

##### 4.4.1.4 Genotyping quality metrics: SNP quality (QUAL)

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_QUAL_hist_2.png)
Counts of SNPs according to QUAL values. The blue dotted line indicates the threshold used for filtering the vcf: QUAL > 200. X axis is on a log scale log(200)=2.3.

-----------------------

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_QUAL_ECDF_2.png)
Empirical cumulative distribution function of QUAL values. The blue dotted line indicates the threshold used for filtering the vcf: QUAL > 200. X axis is on a log scale log(200)=2.3.

-----------------------

<div style="page-break-after: always"></div>

##### 4.4.1.5 Genotyping quality metrics: quality by depth (QD)

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_QD_hist_2.png)
Counts of SNPs according to QUAL values. The blue dotted line indicates the threshold used for filtering the vcf: QD < 20.

-----------------------

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_QD_ECDF_2.png)
Empirical cumulative distribution function of QUAL. values. The blue dotted line indicates the threshold used for filtering the vcf: QD < 20.

-----------------------

<div style="page-break-after: always"></div>

##### 4.4.1.6. Individual genotyping metrics: heterozygote calls

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_het_hist_2.png)
Counts of SNPs according to proportion of heterozygote genotypes. The blue dotted line indicates the threshold used for filtering the vcf: heterozygote genotypes for a SNP must be < 1%.

-----------------------

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_het_ECDF_2.png)
Empirical cumulative distribution function of heterozygote genotypes. The blue dotted line indicates the threshold used for filtering the vcf: heterozygote genotypes for a SNP must be < 1%.

-----------------------

<div style="page-break-after: always"></div>

##### 4.4.1.7. Individual genotyping metrics: missing genotypes

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_miss_hist_2.png)
Counts of SNPs according to proportion of missing genotypes. The blue dotted line indicates the threshold used for filtering the vcf: missing genotypes for a SNP must be < 5%.

-----------------------

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_miss_ECDF_2.png)
Empirical cumulative distribution function of missing genotypes. The blue dotted line indicates the threshold used for filtering the vcf: missing genotypes for a SNP must be < 5%.

-----------------------

<div style="page-break-after: always"></div>

##### 4.4.1.8. Individual genotyping metrics: genotype quality (QD)
QD is a quality measure of each individual genotype. Our filter removes all SNPs having more than 20 % of genotypes having QD values < 10. In the remaining SNPs, genotypes with QD values < 10 are retained.

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_GQ_GQfiltered1_2.png)
Counts of all individual GQ values for all SNPs. The blue dotted line indicates the threshold used: GQ < 10.

-----------------------

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_GQ_GQfiltered2_2.png)
Empirical cumulative distribution function of GQ values for all SNPs. The blue dotted line indicates the threshold used: GQ < 10.

-----------------------

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_GQ_GQfiltered3_2.png)
SNP counts according to proportion of genotypes with GQ < 10. The blue dotted line indicates the threshold used for filtering the vcf: proportion of genotypes for a SNP with GQ < 10%, must be < 20%.

-----------------------

![](SeqApiPop_2_VcfCleanup.assets/plot_decision_GQ_GQfiltered4_2.png)
Empirical cumulative distribution function of SNP according to their proportion of genotypes with GQ < 10. The blue dotted line indicates the threshold used for filtering the vcf: proportion of genotypes for a SNP with GQ < 10%, must be < 20%.

-----------------------

##### 4.4.1.9 Notice on resulting missing data.
although SNPs with more than 5% of missing data were filtered out, some markers may have more than 5% missing data due to the heterozygote calls still remaining after the heterozygote filter, that were set to missing. These can be removed if needed during further analyses.

#### 4.4.2 Running the filters: generate Venn diagrams and filtered vcf
Running the following script with run='filter_all', will perform the filtering and produce the Venn diagrams and filtered vcf file.

```bash
#! /bin/bash

#run_vcf_cleanup.sh

#vcf filter (only done once to prepare list SNP positions)

# modules #####################################################
module load system/R-3.5.1
module load bioinfo/bcftools-1.6
module load bioinfo/tabix-0.2.5
module load bioinfo/vcftools-0.1.15
module load bioinfo/samtools-1.8
# end modules #################################################

# parameters ##################################################
username=avignal #Deprecated
SCRIPTS='~/seqapipopOnHAV3_1/vcf_cleanup_scripts' #path to scripts
DIRIN='~/combineGVCFs/The870vcf' #path to directory containing the input vcf
DIROUT='~/seqapipopOnHAV3_1/vcf_cleanup' #path to output directory
VCFIN='MetaGenotypesCalled870_raw_snps.vcf.gz' #inputvcf before filters
limit_allele=3
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
quantile_prob_above_threshold=0.1
quantile_prob_below_threshold=0.9
kept_above_threshold="MQ_QUAL_QD~GQ~GQ"
kept_below_threshold="FS_SOR_allele~miss_het~GQfiltered"
type='haploid'
run='filter_all'
# end parameters ##############################################
sbatch -W -J vcf_filter -o ${DIROUT}/log/vcf_cleanup.o -e ${DIROUT}/log/vcf_cleanup.e \
		--wrap="${SCRIPTS}/vcf_cleanup.sh ${username} ${SCRIPTS} ${DIRIN} ${DIROUT} ${VCFIN} \
		${limit_allele} ${limit_FS} ${limit_SOR} ${limit_MQ} ${limit_MQRankSum} ${limit_ReadPosRankSum} \
		${limit_QUAL} ${limit_QD} ${limit_GQ} ${limit_miss} ${limit_het} ${limit_GQfiltered} \
		${quantile_prob_above_threshold} ${quantile_prob_below_threshold} \
		${kept_above_threshold} ${kept_below_threshold} ${run} ${type}"

# end of file
```

## 5. results:

<p align="center">
  <img src="SeqApiPop_2_VcfCleanup.assets/d450393d_2.png" />
</p>

**Filters on mapping quality (MQ) and strand bias (FS and SOR) metrics.**

FS (FisherStrand): phred-scaled probability that there is strand mapping bias at the site; SOR (StrandOddsRatio): strand bias mapping estimate; MQ (RMSMappingQuality): root mean square mapping quality over all the reads at the site. The intercept of the 3 filters gives 10,057,214 SNPs, used for further filtering

-----------------------

<div style="page-break-after: always"></div>

<p align="center">
  <img src="SeqApiPop_2_VcfCleanup.assets/b5ea5cde_2.png" />
</p>

**Filters on genotyping quality.**

Int1 is the intersect of the mapping quality filters. QUAL: Phred-scaled quality score for the assertion made in ALT: the more samples have the ATL allele, the higher the QUAL score. QD: quality score normalized by allele depth in which only informative reads are counted. The intercept of the 2 filters with the previous mapping quality filters gives 8,175,852 SNPs, used for further filtering (see figure 3).

-----------------------

<div style="page-break-after: always"></div>

<p align="center">
  <img src="SeqApiPop_2_VcfCleanup.assets/89e031ea_2.png" />
</p>

**Filters individual genotyping quality.**

Int2 is the intersect of the previous filters. Filters are (i) het: proportion of heterozygote calls less than 1% for a SNP, as haploid drones were sequenced, the remaining heterozygote calls were set to missing; (ii) allele: less than 4 alleles for a SNP; (iii) miss: less than 5% missing data; (iv) GCfiltered: SNPs are removed if more than 20% samples have a genotyping quality (GQ) under 10. Note: although SNPs with more than 5% of missing data were filtered out, some markers may have more than 5% missing data due to the heterozygote calls that were set to missing.

