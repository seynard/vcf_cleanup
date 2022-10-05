#module load system/R-3.5.1
#R
######## functions ######## 
installpackages<-function(package_name){if(package_name %in% rownames(installed.packages()) == FALSE) {install.packages(package_name, repos = "http://cran.us.r-project.org")}}
run_filter<-function(filter_list,filter_name,filter_value,threshold_direction,infile,previous_markers_name,previous_markers){
	f_list<-list()
	name_filter<-vector()
	for(i in 1:length(filter_list)){
		if(filter_list[i]=='allele'){
			f_list[[i]]<-paste0(infile$CHROM[nchar(infile$ALT)<=filter_value[filter_name==filter_list[i]] & grepl('\\*',infile$ALT)==F],'_',infile$POS[nchar(infile$ALT)<=filter_value[filter_name==filter_list[i]] & grepl('\\*',infile$ALT)==F])
		}else if(filter_list[i]%in%filter_name & filter_list[i]!='allele'){
			if(threshold_direction[filter_name==filter_list[i]]=='below'){
				name_filter[i]<-paste0(filter_name[filter_name==filter_list[i]],'<',filter_value[filter_name==filter_list[i]])
				f_list[[i]]<-paste0(infile$CHROM[infile[,filter_list[i]]<filter_value[filter_name==filter_list[i]]],'_',infile$POS[infile[,filter_list[i]]<filter_value[filter_name==filter_list[i]]])
			}else if (threshold_direction[filter_name==filter_list[i]]=='above'){
				f_list[[i]]<-paste0(infile$CHROM[infile[,filter_list[i]]>filter_value[filter_name==filter_list[i]]],'_',infile$POS[infile[,filter_list[i]]>filter_value[filter_name==filter_list[i]]])
				name_filter[i]<-paste0(filter_name[filter_name==filter_list[i]],'>',filter_value[filter_name==filter_list[i]])
			}
		}
	}
	f_list<-f_list[lapply(f_list,length)>0] 
	return(f_list)
}
plot_venn<-function(filter_list,infile,previous_markers_name,previous_markers){
	f_list<-infile
	name_filter<-name_filter[!is.na(name_filter)]
	if(previous_markers_name!='none'){f_list[[length(f_list)+1]]<-(previous_markers)
		name_filter<-c(name_filter,previous_markers_name)}
	if('allele'%in%filter_list){name_filter<-c('allele',name_filter)}
	myCol<-brewer.pal(length(f_list),"Set2")
	plot<-venn.diagram(
		x=f_list,
		category.names=name_filter,
		filename=NULL,
		output=TRUE,
		cex=0.5,
		lwd=2,
		col=myCol,
		fill=myCol,
		alpha=0.3,
		cat.cex=1,
		margin=c(0.3,0.3,0.3,0.3))
	pdf(paste0(dir,'/venn',paste0(name_filter,collapse='_'),'.pdf'),width=10,height=10)
	grid.draw(plot)
	dev.off()
}
make_list<-function(infile){
	f_list<-infile
	int<-Reduce(intersect,f_list)
	return(int)
}
###########################  
######## parameters ######## 
args<-commandArgs(TRUE)
dir<-args[1]
kept_above_threshold<-args[2]
args_above=args[3]
kept_below_threshold<-args[4]
args_below=args[5]
###########################  
package_list<-c('data.table','VennDiagram','reshape2','RColorBrewer','grDevices')
for(i in 1:length(package_list)){
	installpackages(package_list[i])
	library(package_list[i],character.only=T)}

ABOVE<-data.frame(filt=unlist(strsplit(kept_above_threshold,'_')),Filt=unlist(strsplit(args_above,'_')),info='above')
BELOW<-data.frame(filt=unlist(strsplit(kept_below_threshold,'_')),Filt=unlist(strsplit(args_below,'_')),info='below')
filter_info<-rbind(ABOVE,BELOW)
filter_info<-filter_info[!duplicated(filter_info),]
print(filter_info)
filter_info$Filt<-as.numeric(as.character(filter_info$Filt))

info<-fread(paste0(dir,'/info.txt'),header=T,data.table=F)
if('miss'%in%filter_info$filt){miss<-fread(paste0(dir,'/missing.txt'),header=T,data.table=F)
	info$miss<-miss$Miss}
if('het'%in%filter_info$filt){het<-fread(paste0(dir,'/heterozygous.txt'),header=T,data.table=F)
	info$het<-het$Het}
if('GQ'%in%filter_info$filt & 'GQfiltered'%in%filter_info$filt){gq<-fread(paste0(dir,'/GQfiltered.txt'),header=T,data.table=F)
	info$GQfiltered<-gq[,3]}
for(i in 5:ncol(info)){info[,i][info[,i]=='.']<-0; info[,i]<-as.numeric(info[,i])}	
#first venn (mapping quality filter)
filter_list<-c('FS','SOR','MQ','MQRankSum','ReadPosRankSum')
f_list1<-run_filter(filter_list,filter_info$filt,filter_info$Filt,filter_info$info,info,'none',NA)
plot_venn1<-plot_venn(f_list1,info,'none',NA)
int1<-make_list(f_list1)
#second venn quality
filter_list<-c('QUAL','QD')
f_list2<-run_filter(filter_list,filter_info$filt,filter_info$Filt,filter_info$info,info,'int1',int1)
plot_venn2<-plot_venn(f_list2,info,'int1',int1)
int2<-make_list(f_list2)
#third venn allele
filter_list<-c('allele','miss','het','GQfiltered')
f_list3<-run_filter(filter_list,filter_info$filt,filter_info$Filt,filter_info$info,info,'int2',int2)
plot_venn3<-plot_venn(f_list3,info,'int2',int2)
int3<-make_list(f_list3)
list_final<-colsplit(int3,'.1_',c('CHROM','POS'))
list_final$CHROM<-paste0(list_final$CHROM,'.1')
write.table(list_final,paste0(dir,'/list_kept.txt'),col.names=F,row.names=F,quote=F,sep='\t')

print('count initial VCF')
print(data.table(table(info$CHROM)))
print(nrow(info))

list<-colsplit(int1,'_',c('NC','CHROM','POS'))
list$CHROM<-paste0(list$NC,'_',list$CHROM)
list$NC<-NULL
print('count VCF GATK filters')
print(data.table(table(list$CHROM)))
print(length(int1))

list<-colsplit(int2,'_',c('NC','CHROM','POS'))
list$CHROM<-paste0(list$NC,'_',list$CHROM)
list$NC<-NULL
print('count VCF quality filters')
print(data.table(table(list$CHROM)))
print(length(int2))

list<-colsplit(int3,'_',c('NC','CHROM','POS'))
list$CHROM<-paste0(list$NC,'_',list$CHROM)
list$NC<-NULL
print('count VCF allele and biological filters')
print(data.table(table(list$CHROM)))
print(length(int3))
