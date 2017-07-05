#############################
## normal fisher analysis  ##
#############################


gwas.data.sorted <- read.table("/Users/beritlilje/Documents/Projects/SAB_vs_IE/SAB_IE_vs4_onlyLeftSide/SNP_analysis/NASP/NASP_output/bestsnp_ALL_newRownames_newColnames_0and1.tsv", head=T, stringsAsFactors=F)

# set number of IE and SAB samples, OBS! assumes that data is sorted with all IE samples first, followed by all SAB-only samples
ie_counts<-120
sab_counts<-121

m<-as.matrix(gwas.data.sorted)
mt<-t(m)

groups<-c(rep("ie", ie_counts), rep("sab", sab_counts))  
data.frame(rownames(mt), groups) # to see if groups are correct


snps <- mt
phen <- factor(groups)


pval <- apply(snps, 2, function(e) fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)
min(pval)
#[1] 0.001818291
length(which(pval < 0.05))
# [1] 609
pval.corrected <- p.adjust(pval, method="fdr")
min(pval.corrected)
#[1] 1


gwas.data.sorted$pval<-pval
gwas.data.sorted$fdr<-pval.corrected
f<-gwas.data.sorted[sort.list(gwas.data.sorted$pval),]
f[1:10, (ncol(f)-1):ncol(f)]
#                  pval fdr
#SNP833306  0.001818291   1
#SNP833307  0.001818291   1
#SNP2747511 0.002315315   1
#SNP2577546 0.002475605   1
#SNP2827154 0.002838438   1
#SNP2827157 0.002838438   1
#SNP2827181 0.002838438   1
#SNP141816  0.002901996   1
#SNP645920  0.003570054   1
#SNP986298  0.003570054   1

# investigate SB IE proportions for top hit:
table(factor(f[1,1:(ncol(f)-2)],levels=c(0,1)), phen)
#     ie sab
#  0  94 112
#  1  26   9




##########################################
## Plot SNPs combinations on power plot ##
##########################################



gwas.data.sorted <- read.table("/Users/beritlilje/Documents/Projects/SAB_vs_IE/SAB_IE_vs4_onlyLeftSide/SNP_analysis/NASP/NASP_output/bestsnp_ALL_newRownames_newColnames_0and1.tsv", head=T, stringsAsFactors=F)
m<-as.matrix(gwas.data.sorted)
ie_counts<-120
sab_counts<-121


# Find the SAB vs IE counts that are present in this dataset
ie<-m[,1:ie_counts]
sab<-m[,(ie_counts+1):(ie_counts+sab_counts)]

ie_sum<-apply(ie, 1, sum)
sab_sum<-apply(sab, 1, sum)



########## The power plot #########

number_of_SNPs<-nrow(gwas.data.sorted)
# 120850 SNPs


p_value<-data.frame(t=(1:(sab_counts+1)))

# calculate p-values for all combinations af SAB and IE samples:
for(i in 0:ie_counts){
  p_string=c()
  for(ii in 0:sab_counts){
    p<-fisher.test(matrix(c(i,ii, ie_counts-i, sab_counts-ii), nrow=2))$p.value
    p_string<-c(p_string, p)
  }
  p_value<-cbind(p_value, p_string)
}


p_value<-p_value[,-1]
dim(p_value)

row.names(p_value)<-(paste("SAB_SNP", 0:sab_counts, sep="_"))
names(p_value)<-(paste("IE_SNP", 0:ie_counts, sep="_"))


##### Bin power plot to, larger than one, bwtween 0.05 and 1, and smaller than 0.05
fdr_frame<-p_value*number_of_SNPs # correct for mulitple testing
fdr_frame[fdr_frame>=1]<-2
fdr_frame[fdr_frame<0.05]<-0
fdr_frame[fdr_frame>0 & fdr_frame<2]<-1


fdr_frame_new<-fdr_frame

# add all our combinations (does not have to do with p-values)
for(i in 1:length(ie_sum)){
  ie_temp<-as.numeric(ie_sum[i])+1 # since the first column is the 0 column, 0 hits should get index 1 etc. 
  sab_temp<-as.numeric(sab_sum[i])+1
  fdr_frame_new[sab_temp,ie_temp]<-3
}




library(gplots)
my_palette <- colorRampPalette(c("red", "orange", "light yellow", "green"))(n = 4) # Frame consist of: 0=sign, 1=0.05 to 1, 2=not sign and 3 =our combinataions
heatmap.2(as.matrix(fdr_frame_new), Rowv=FALSE, Colv=FALSE, col=my_palette,dendrogram="none", scale="none", trace="none", keysize = 0.1, key=F, cexRow=0.5, cexCol=0.5, margins=c(5,10), colsep=1:120, rowsep=1:121)




