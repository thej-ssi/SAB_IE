#####################################################
## PCA, DAPC and population stratification on SNPs ##
#####################################################


### DAPC results ###
library(adegenet)

gwas.data.sorted <- read.table("/Users/beritlilje/Documents/Projects/SAB_vs_IE/SAB_IE_vs4_onlyLeftSide/SNP_analysis/NASP/NASP_output/bestsnp_ALL_newRownames_newColnames_0and1.tsv", head=T, stringsAsFactors=F)

# set number of IE and SAB samples, OBS! assumes that data is sorted with all IE samples first, followed by all SAB-only samples
ie_counts<-120
sab_counts<-121


m<-as.matrix(gwas.data.sorted)
mt<-t(m)

groups<-c(rep("ie", ie_counts), rep("sab", sab_counts))  
data.frame(rownames(mt), groups) # to see if groups are correct


#### DAPC #####
snps <- mt
phen <- factor(groups)

241/3 # to find maximal value for dapc (number of PCs should not be higher than 1/3 of number of samples)
#[1] 80.33333 # number of maximum PC, sohuld be n/3
v<-c(1,2,3,4,5,10,15,20,25,30,40,50,60,70, 80)
set.seed(1)
cross<-xvalDapc(snps, phen, n.pca=v, n.rep=100) 

cross
#$`Median and Confidence Interval for Random Chance`
#     2.5%       50%     97.5% 
#0.4273760 0.4937672 0.5518595 
#
#$`Mean Successful Assignment by Number of PCs of PCA`
#        1         2         3         4         5        10        15        20 
#0.5275000 0.5008333 0.5070833 0.4925000 0.4900000 0.4725000 0.4287500 0.4729167 
#       25        30        40        50        60        70        80 
#0.4729167 0.4612500 0.5029167 0.4695833 0.5170833 0.5179167 0.5145833 
#
#$`Number of PCs Achieving Highest Mean Success`
#[1] "1"
#
#$`Root Mean Squared Error by Number of PCs of PCA`
#        1         2         3         4         5        10        15        20 
#0.4815340 0.5053876 0.5031671 0.5145116 0.5170697 0.5357368 0.5783567 0.5345851 
#       25        30        40        50        60        70        80 
#0.5334798 0.5478335 0.5060914 0.5389515 0.4912251 0.4915077 0.4966380 
#
#$`Number of PCs Achieving Lowest MSE`
#[1] "1"

# Don't have to run this step, as no optimal number of PCs were better than random.
dapc1<-dapc(snps, phen, n.pca=1, n.da=1) # PC=1 from crossvalidation, only 1 group, so only one axes to retain in discriminatory analysis
scatter(dapc1)  



####### PCA #######
cc_type<-read.table("/Users/beritlilje/Documents/Projects/SAB_vs_IE/SAB_IE_vs5_focus_on_whole_set/Final_figures/meta_data_sheet.txt", h=T)
row.names(cc_type)<-cc_type$Project_ID
cc_type_ordered<-cc_type[row.names(snps),] # get data in same order 
cc_type_ordered$CC_type<-factor(cc_type_ordered$CC_type)

pca1 <- dudi.pca(snps, scale=FALSE) # chose 2 axes, only choose scale=TRUE if you are measuring 2 different things, e.g. hight and weight

s.class(pca1$li, fac=phen, col=funky(2), xlim=c(-150,150), ylim=c(-100,100))
add.scatter.eig(pca1$eig[1:50],2,1,2, ratio=.3)

s.class(pca1$li, fac=cc_type_ordered$CC_type, col=funky(23), xlim=c(-150,150), ylim=c(-100,100)) #,xlim=c(-40,60), ylim=c(-40,60))
add.scatter.eig(pca1$eig[1:50],2,1,2, ratio=.3)


# find contributions (for axes in PCA plot):
round(pca1$eig[1]/sum(pca1$eig)*100, 1)
#[1] 40.6
round(pca1$eig[2]/sum(pca1$eig)*100, 1)
#[1] 24.3



###########################################
##### With population stratification ######
###########################################


# correct for CC type
set.seed(1)
snps.corrected <- apply(snps, 2, function(e) residuals(lm(e~cc_type_ordered$CC_type))) # choose number of groups that was seen in PCA+hclust

#### Make PCA plots ####
pca_corrected <- dudi.pca(snps.corrected, scale=FALSE) # chose 5 axes

s.class(pca_corrected$li, fac=phen, col=funky(2), xlim=c(-40,40), ylim=c(-40,40))
add.scatter.eig(pca_corrected$eig[1:50],3,1,2, ratio=.3)

s.class(pca_corrected$li, fac=cc_type_ordered$CC_type, col=funky(23),xlim=c(-40,40), ylim=c(-40,40))
add.scatter.eig(pca_corrected$eig[1:50],3,1,2, ratio=.3)



#### DAPC ####
set.seed(1)
cross<-xvalDapc(snps.corrected, phen, n.pca=v, n.rep=100)
#$`Median and Confidence Interval for Random Chance`
#     2.5%       50%     97.5% 
#0.4356749 0.5020661 0.5601584 
#
#$`Mean Successful Assignment by Number of PCs of PCA`
#        1         2         3         4         5        10        15        20 
#0.5095833 0.4970833 0.5150000 0.5391667 0.5362500 0.5154167 0.5166667 0.5058333 
#       25        30        40        50        60        70        80 
#0.5141667 0.4916667 0.5262500 0.5558333 0.5295833 0.5262500 0.4825000 
#
#$`Number of PCs Achieving Highest Mean Success`
#[1] "50"
#
#$`Root Mean Squared Error by Number of PCs of PCA`
#        1         2         3         4         5        10        15        20 
#0.4917196 0.5054392 0.4892951 0.4685972 0.4710914 0.4911190 0.4914194 0.5026665 
#       25        30        40        50        60        70        80 
#0.4942025 0.5137350 0.4825964 0.4540329 0.4786592 0.4832795 0.5241726 
#
#$`Number of PCs Achieving Lowest MSE`
#[1] "50"
dapc2<-dapc(snps.corrected, phen, n.pca=50,n.da=1) # not necessary, as one can see from PCA plot that no PCs will show values above the CI
scatter(dapc2)  



#### T-test on stratefied data (can't use Fisher test anymore) ######
pval <- apply(snps.corrected, 2, function(e) t.test(e[1:ie_counts], e[sab_counts:(ie_counts+sab_counts)])$p.value)
min(pval)
#[1] 0.01686422
pval.corrected <- p.adjust(pval, method="fdr")
min(pval.corrected)
#[1] 0.8191223

# no reason to make plot, nothing significant, but can list the top hits anyway
out<-data.frame(colnames(snps.corrected), pval, pval.corrected) # obs, here corrected is mulitple testing correction, all are corrected for pop stratification
f<-out[sort.list(out$pval),]
f[1:10,]

#colnames.snps.corrected.       pval pval.corrected
#              SNP2285256 0.01686422      0.8191223
#              SNP1629428 0.02106782      0.8191223
#               SNP974358 0.02150012      0.8191223
#              SNP1687148 0.02150012      0.8191223
#              SNP1823878 0.02150012      0.8191223
#              SNP2052218 0.02150012      0.8191223
#              SNP2365453 0.02150012      0.8191223
#              SNP1940061 0.02150012      0.8191223
#               SNP954134 0.02303450      0.8191223
#              SNP2210911 0.02633183      0.8191223




#####  mann-Whitney-U test (non-parametric t-test) ####
pval <- apply(snps.corrected, 2, function(e) wilcox.test(e[1:ie_counts], e[sab_counts:(ie_counts+sab_counts)])$p.value)
min(pval)
#[1] 0.001409708

pval.corrected <- p.adjust(pval, method="fdr")
min(pval.corrected) 
#[1] 0.9114139


##### Yet another method (something from a dapc tutorial) #####
set.seed(1)
pval2 <- numeric(0)
for(i in 1:ncol(snps.corrected)){
	foo <- suppressWarnings(glm(phen ~ snps.corrected[,i], family="binomial")) 
	ANOVA <- anova(foo, test="Chisq")
	pval2[i] <- ANOVA$"Pr(>Chi)"[2]
}

min(pval2)
#[1] 0.003002037

pval2.corrected <- p.adjust(pval2, method="fdr")
min(pval2.corrected)
#[1] 0.7446625



