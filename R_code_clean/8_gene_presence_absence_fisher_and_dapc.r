############################
## Gene presence abscence ##
############################


# Assmeble data with spades
# Run Prokka with default settings
# Run Roary with options: -e, -n, -s, -z turned on 



# get a list of accessory genes, in terminal:
grep "gene" accessory.tab > accessory_genes.txt


# in R
d<-read.table("Outputs/Roary_output/gene_presence_absence.Rtab", h=T)
accessory<-read.table("Outputs/Roary_output/accessory_genes.txt")

# set number of IE and SAB samples, OBS! assumes that data is sorted with all IE samples first, followed by all SAB-only samples
ie_counts<-120
sab_counts<-121
row.names(d)<-d[,1]
d<-d[,-1]

#### only keep accessory ####
library(reshape)
a_genes<-colsplit(accessory$V2, split="=", names=c("not_relevant", "gene_name"))
d_accessory<-d[row.names(d) %in% a_genes$gene_name,]


nrow(d)
#[1] 5508
nrow(d_accessory)
#[1] 3326
nrow(accessory)
#[1] 3326

# calcualte 
m<-as.matrix(d_accessory)
mt<-t(m)

groups<-c(rep("ie", ie_counts), rep("sab", sab_counts))  
data.frame(rownames(mt), groups) # to see if groups are correct

snps <- mt
phen <- factor(groups)


##### fisher test #####
pval <- apply(snps, 2, function(e) fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)
min(pval)
#[1]  0.004938278
length(which(pval < 0.05))
#[1] 40
pval.corrected <- p.adjust(pval, method="fdr")
min(pval.corrected)
#[1] 1
f<-data.frame(gene=d_accessory, pval=pval, fdr=pval.corrected)


f_sorted<-f[sort.list(f$pva,),]
head(f_sorted) 



####### DAPC #######

library(adegenet)

v<-c(1,2,3,4,5,10,15,20,25,30,40,50,60,70, 80)
cross<-xvalDapc(snps, phen, n.pca=v, n.rep=100) # no PC producing better hits than random, so can stop after this step.

#$`Mean Successful Assignment by Number of PCs of PCA`
#        1         2         3         4         5        10        15        20 
#0.5283333 0.4875000 0.4900000 0.4775000 0.4812500 0.4975000 0.4662500 0.4591667 
#       25        30        40        50        60        70        80 
#0.4616667 0.4933333 0.4662500 0.4483333 0.4458333 0.5070833 0.5070833 
#
#$`Number of PCs Achieving Highest Mean Success`
#[1] "1"
#
#$`Root Mean Squared Error by Number of PCs of PCA`
#        1         2         3         4         5        10        15        20 
#0.4831178 0.5198157 0.5167003 0.5302973 0.5256775 0.5100041 0.5416186 0.5471200 
#       25        30        40        50        60        70        80 
#0.5459446 0.5157249 0.5427393 0.5590481 0.5613723 0.5020963 0.5015081 
#
#$`Number of PCs Achieving Lowest MSE`
#[1] "1"

# all PCs are at around 0.5, which means no signal 



