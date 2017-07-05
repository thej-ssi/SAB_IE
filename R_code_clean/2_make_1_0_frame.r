########################
## Make 1 and 0 frame ##
########################

#### make SNP matrix in correct format, with only 1 and 0 #####
# Make frame that contains 0 if equal to reference, and 1 otherwise
# remove sites that does not contain real SNPs (all samples have the same base which is different from reference)

#Import data file format from 1_vs5_clean_up_data.r
gwas_raw_data <- read.table("Outputs/bestsnp_ALL_newRownames_newColnames.tsv", head=T, stringsAsFactors=F)
gwas_raw_data$Position <- sub("^", "SNP", gwas_raw_data$Position)
set.seed(1)


test<-gwas_raw_data
ref<-test[,2]
t<-data.frame(ref=ref)
for(i in 3:ncol(test)){
	v<-test[,i]
	v[v==ref]<-0
	v[!v==0]<-1
	t<-cbind(t, as.numeric(v))
}

row.names(t)<-test$Position
names(t)<-names(test)[-1]



### Remove rows that does not contain SNPs anymore, are 116 SNP, were all samples have same base but they all differ from reference

gwas.data<-t[,-1]
dim(gwas.data)
#[1] 120985    241
gwas.data$snp_counts<-apply(gwas.data, 1, sum)
max(gwas.data$snp_counts)
#[1] 241
min(gwas.data$snp_counts)
#[1] 1 #no with 0 (if no SNPs, and equal to ref, which should not exists) 
gwas.data.sorted<-gwas.data[!gwas.data$snp_counts==241,]  # all same as ref is not a SNP, so remove cases were all are diff (1) OBS dependent on set!!!
gwas.data.sorted<-gwas.data.sorted[,-ncol(gwas.data.sorted)]
nrow(gwas.data.sorted)
# [1] 120850

write.table(gwas.data.sorted, "newoutput/bestsnp_ALL_newRownames_newColnames_0and1.tsv", quote=F, row.names=T, col.names=T, sep="\t")




