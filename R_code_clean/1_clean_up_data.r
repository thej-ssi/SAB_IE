################################################
## change row and column names ++ of SNP data ##
################################################

# Load bestsnp.tsv file from nasp run with standard settings
setwd("S:/MPV/Projekter/10 Staphylococcus projects/10_5 S. aureus virulence/10_5_3_SAB vs SA-IE/Data_for_github_upload")
d<-read.table("Outputs/bestsnp_all.tsv", h=T, comment.char="")


#Test number of different bases at each snp-position
bases<-d[,250:253]
zeros<-c()
for(i in 1:nrow(bases)){
	temp<-bases[i,]
	t<-sum(temp ==0)
	zeros<-c(zeros, t)
}

#How many different bases at each SNP position?
sum(zeros==0)
#[1] 65 SNP positions where all four bases are found in at least one sample
sum(zeros==1)
#[1] 3805 SNP positions where 3 different bases are found
sum(zeros==2)
#[1] 116985 SNP positions where 2 different bases are found
sum(zeros==3)
#[1] 130 SNP positions where only one type of base is found. Irrelevant SNPs where all samples have the same base, but they all differ from the reference sequence
sum(zeros==4)
#[1] 0

# percent of SNPs that have more than 2 bases:
(65+3805)/(65+3805+116985)*100
#3.2%

# rename rows to the SNPs position in reference sequence
library(reshape)
position<-colsplit(d$LocusID, split="::", names=c("gene", "position"))
dd<-d[,2:243] # remove position and non-important last columns
new<-cbind(position$position, dd) # add position

# rename columns to contain only sample name
id<-names(dd)[2:ncol(dd)] # first is position and reference
new_id<-colsplit(id, split="_S", names=c("ID", "others"))
colnames_vec<-c("Position","Reference", as.character(new_id$ID))
colnames(new)<-colnames_vec


#print reformatted data table
write.table(new, "newoutput/bestsnp_ALL_newRownames_newColnames.tsv", quote=F, row.names=F, col.names=T, sep="\t")



