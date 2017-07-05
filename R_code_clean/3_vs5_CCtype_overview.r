 ###########################
 ## sab ie CC types stats ##
 ###########################

library(plyr)

#d<-read.table("/Users/beritlilje/Documents/Projects/SAB_vs_IE/SAB_IE_vs5_focus_on_whole_set/Final_figures/meta_data_sheet.txt", sep="\t", h=T)
d<-read.table("meta_data_sheet.txt", sep="\t", h=T)

### Make dataframe with CC type for SAB and IE samples
sab<-d[d$IE==0,]
ie<-d[d$IE==1,]

sab_cc<-ddply(sab, .(final_CC_type), nrow)
names(sab_cc)<-c("CC", "SAB")
ie_cc<-ddply(ie, .(final_CC_type), nrow)
names(ie_cc)<-c("CC", "IE")

m<-merge(sab_cc, ie_cc, by="CC", all=T)
m[is.na(m)]<-0
m$total<-m$SAB+m$IE


### calculate percentage for each CC type ###
sab_total<-sum(m$SAB)
ie_total<-sum(m$IE)
pvals<-c()
for(i in 1:nrow(m)){
	sab<-m[i,2]
	ie<-m[i,3]
	mat<-matrix(c(sab, ie, sab_total-sab, ie_total-ie), nrow=2, byrow=T)
	pvals<-c(pvals, fisher.test(mat)$p.value)
}



### add p-values to dataframe ###
m$pvals<-pvals
m$pvals_round<-round(m$pvals, 2)


m$sab_percent<-round(m$SAB/sum(m$SAB)*100,1)
m$ie_percent<-round(m$IE/sum(m$IE)*100,1)
ms<-m[sort.list(m$total, decreasing=T),]
ms[,c(1,2,7,3,8,6)]

#       CC SAB sab_percent IE ie_percent pvals_round
#       45  28        23.1 24       20.0        0.64
#       30  22        18.2 23       19.2        0.87
#       15  16        13.2 15       12.5        1.00
#        1  10         8.3  9        7.5        1.00
#        5  10         8.3  7        5.8        0.62
#        8   4         3.3 10        8.3        0.11
#       12   2         1.7  6        5.0        0.17
#       22   5         4.1  3        2.5        0.72
#       25   3         2.5  4        3.3        0.72
#       59   5         4.1  2        1.7        0.45
#        7   3         2.5  2        1.7        1.00
#      188   0         0.0  4        3.3        0.06
#        9   2         1.7  2        1.7        1.00
#       20   2         1.7  1        0.8        1.00
#      509   2         1.7  1        0.8        1.00
#       97   1         0.8  2        1.7        0.62
#      121   2         1.7  0        0.0        0.50
#      182   1         0.8  1        0.8        1.00
#      573   1         0.8  1        0.8        1.00
#        6   0         0.0  2        1.7        0.25
#       50   0         0.0  1        0.8        0.50
#       72   1         0.8  0        0.0        1.00
#singleton   1         0.8  0        0.0        1.00

# copy to word document, add paranthesis etc. repeat for CA set.




