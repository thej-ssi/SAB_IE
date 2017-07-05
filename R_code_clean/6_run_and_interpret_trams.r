#############################
## Run and interpret TRAMS ##
#############################

# download trams, e.g. https://sourceforge.net/projects/strathtrams/
python TRAMS.py -s Outputs/bestsnp_All_newRownames_newColnames.tsv -r Outputs/NC_021554.gbk

# choose which columns to save for further analysis:
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$6,$7,$13,$18,$23,$28,$33,$38,$43,$48,$53,$58,$63,$68,$73,$78,$83,$88,$93,$98,$103,$108,$113,$118,$123,$128,$133,$138,$143,$148,$153,$158,$163,$168,$173,$178,$183,$188,$193,$198,$203,$208,$213,$218,$223,$228,$233,$238,$243,$248,$253,$258,$263,$268,$273,$278,$283,$288,$293,$298,$303,$308,$313,$318,$323,$328,$333,$338,$343,$348,$353,$358,$363,$368,$373,$378,$383,$388,$393,$398,$403,$408,$413,$418,$423,$428,$433,$438,$443,$448,$453,$458,$463,$468,$473,$478,$483,$488,$493,$498,$503,$508,$513,$518,$523,$528,$533,$538,$543,$548,$553,$558,$563,$568,$573,$578,$583,$588,$593,$598,$603,$608,$613,$618,$623,$628,$633,$638,$643,$648,$653,$658,$663,$668,$673,$678,$683,$688,$693,$698,$703,$708,$713,$718,$723,$728,$733,$738,$743,$748,$753,$758,$763,$768,$773,$778,$783,$788,$793,$798,$803,$808,$813,$818,$823,$828,$833,$838,$843,$848,$853,$858,$863,$868,$873,$878,$883,$888,$893,$898,$903,$908,$913,$918,$923,$928,$933,$938,$943,$948,$953,$958,$963,$968,$973,$978,$983,$988,$993,$998,$1003,$1008,$1013,$1018,$1023,$1028,$1033,$1038,$1043,$1048,$1053,$1058,$1063,$1068,$1073,$1078,$1083,$1088,$1093,$1098,$1103,$1108,$1113,$1118,$1123,$1128,$1133,$1138,$1143,$1148,$1153,$1158,$1163,$1168,$1173,$1178,$1183,$1188,$1193,$1198,$1203,$1208,$1213}' Outputs/Trams_output/annotation_bestsnp_All_newRownames_newColnames.tsv > Outputs/Trams_output/annotation_bestsnp_All_selectedColumns.txt

# open output file in excel, replace "SNP type" with isolate name, remove "reference" and the first row, rename file to annotation_bestsnp_All_selectedColumns_cleaned.txt




R

####################
### Prepare data ###
####################

#d<-read.table("annotation_bestsnp_All_selectedColumns.txt", sep="\t", h=T)
d<-read.table("Outputs/Trams_output/annotation_bestsnp_All_selectedColumns_cleaned.txt", sep="\t", h=T)
rem<-c("intergenic", "misc_binding") # remove intergenic regions:
dd<-d[!d$Feature %in% rem,]

ie_counts<-120
sab_counts<-121





################################
## All samples, all mutations ##
################################

#counts number of mutations within each gene per sample
all_mut_frame<-matrix(rep(0,(ie_counts+sab_counts)), nrow=1) # initialize matrix
tags<-unique(dd$Locus.tag)

for(id in tags){
 temp<-dd[dd$Locus.tag==id,]
 temp_counts<-apply(temp[,6:ncol(temp)], 2, function(e) sum(e != ""))
 all_mut_frame<-rbind(all_mut_frame, temp_counts)
}

all_mut_frame<-all_mut_frame[-1,]
row.names(all_mut_frame)<-tags



#### remove genes (rows) withour any mutations ####
all_mut_data_frame<-as.data.frame(all_mut_frame)
all_mut_data_frame$sum<-apply(all_mut_data_frame, 1, sum)
all_mut_sorted<-all_mut_data_frame[all_mut_data_frame$sum>0,] 
all_mut_sorted<-all_mut_sorted[,-ncol(all_mut_sorted)]

nrow(all_mut_data_frame)
#[1] 2238
nrow(all_mut_sorted)
#[1] 2233




######## Calculate P values ########


##### P-value for distributions (Mann-Whitney)
groups<-c(rep("ie", ie_counts), rep("sab", sab_counts))  
data.frame(colnames(all_mut_sorted), groups) # to see if groups are correct

snps <- all_mut_sorted
phen <- factor(groups)

p_list <- apply(snps, 1, function(e) wilcox.test(e~phen)$p.value)

pval.corrected <- p.adjust(p_list, method="fdr")
min(pval.corrected)
# [1] 1

# to see top10 hits:
allMut_p<-data.frame(row.names(all_mut_sorted), p_list)
allMut_p$fdr<-p.adjust(allMut_p$p_list, method="fdr")
fs<-allMut_p[sort.list(allMut_p$p_list),]
head(fs)
#              row.names.all_mut_sorted.     p_list fdr
#CA347_RS13060             CA347_RS13060 0.01514491   1
#CA347_RS06035             CA347_RS06035 0.03887100   1
#CA347_RS12985             CA347_RS12985 0.04134306   1
#CA347_RS10405             CA347_RS10405 0.04783807   1
#CA347_RS10410             CA347_RS10410 0.05239829   1
#CA347_RS01195             CA347_RS01195 0.05387098   1


### P-value for proportions (Fisher test) ###
s<-all_mut_sorted
s[s>0]<-1
data.frame(colnames(s), phen)

pval <- apply(s, 1, function(e) fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)

pval.corrected <- p.adjust(pval, method="fdr")
min(pval.corrected)
# [1] 0.9677821







###################################
## All samples, nonsyn mutations ##
###################################

nonsyn_frame<-matrix(rep(0,(ie_counts+sab_counts)), nrow=1)
tags<-unique(dd$Locus.tag)
for(id in tags){
 temp<-dd[dd$Locus.tag==id,]
 temp_counts<-apply(temp[,6:ncol(temp)], 2, function(e) sum(e=="nonsynonymous"))
 nonsyn_frame<-rbind(nonsyn_frame, temp_counts)
}

nonsyn_frame<-nonsyn_frame[-1,]
row.names(nonsyn_frame)<-tags



#### remove genes (rows) withour any mutations ####
nonsyn_data_frame<-as.data.frame(nonsyn_frame)
nonsyn_data_frame$sum<-apply(nonsyn_data_frame, 1, sum)
nonsyn_sorted<-nonsyn_data_frame[nonsyn_data_frame$sum>0,] 
nonsyn_sorted<-nonsyn_sorted[,-ncol(nonsyn_sorted)]

nrow(nonsyn_data_frame) # one per gene. so should give the same numbers for all, nonsende, and non-syn
#[1] 2238
nrow(nonsyn_sorted)
#[1] 2202


######## Calculate P values ########


##### P-value for distributions (Mann-Whitney)

snps <- nonsyn_sorted

p_list <- apply(snps, 1, function(e) wilcox.test(e~phen)$p.value)

pval.corrected <- p.adjust(p_list, method="fdr")
min(pval.corrected)
# [1] 1




### Fisher test ####
s<-nonsyn_sorted
s[s>0]<-1

pval <- apply(s, 1, function(e) fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)

pval.corrected <- p.adjust(pval, method="fdr")
min(pval.corrected)
# [1] 1




########################################
#### All samples, nonsense mutations ###
########################################


nonsense_frame<-matrix(rep(0,(ie_counts+sab_counts)), nrow=1)
tags<-unique(dd$Locus.tag)
for(id in tags){
 temp<-dd[dd$Locus.tag==id,]
 temp_counts<-apply(temp[,6:ncol(temp)], 2, function(e) sum(e=="nonsense"))
 nonsense_frame<-rbind(nonsense_frame, temp_counts)
}

nonsense_frame<-nonsense_frame[-1,]
row.names(nonsense_frame)<-tags


#### remove genes (rows) withour any mutations ####
nonsense_data_frame<-as.data.frame(nonsense_frame)
nonsense_data_frame$sum<-apply(nonsense_data_frame, 1, sum)
nonsense_sorted<-nonsense_data_frame[nonsense_data_frame$sum>0,] #try 0 and 3
nonsense_sorted<-nonsense_sorted[,-ncol(nonsense_sorted)]

nrow(nonsense_data_frame)
#[1] 2238
nrow(nonsense_sorted)
#[1] 409



####### P value ######

### Fisher test ####

s<-nonsense_sorted
s[s>0]<-1

pval <- apply(s, 1, function(e) fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)

pval.corrected <- p.adjust(pval, method="fdr")
min(pval.corrected)
# [1] 1



