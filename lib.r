#R

installpackages<-function(package_name){
	if(package_name %in% rownames(installed.packages()) == FALSE) {
		install.packages(package_name,repos="http://cran.us.r-project.org")
		}
	}

package_list<-c('data.table','ggplot2','reshape2','viridis','venn','RColorBrewer','grDevices','grid','ggpolypath')
for(i in 1:length(package_list)){
	installpackages(package_list[i])
	library(package_list[i],character.only=T)}
