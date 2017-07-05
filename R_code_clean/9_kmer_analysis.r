####################
## K-mer analysis ##
####################

# samples need to have identifier that distinghuish the two groups in the beginning of the sample name.
python kmer_state_reverseComplement.py -i Outputs/Assembled_files/ -o . -p SAB,IE -k 30
# the script outputs a "kmer_state_samples_list.txt" and a "kmer_state_dict.txt", the latter us used for further analysis.
# this script counts k-mers once per sample, and looks trough both forward and reverse compliment of reads, k-mer does not span contig-junctions.

d<-read.table("kmer_state_dict.txt")
names(d)<-c("kmer", "SAB","IE")
nrow(d)
# [1] 11184722

d$pattern<-paste(d$SAB, d$IE, sep="_")

nrow(d)
#[[1] 11184722
length(unique(d$pattern))
# [1] 3553

d_prop<-d[,2:4]
d_prop_unique<-d_prop[!duplicated(d_prop),]
nrow(d_prop_unique)
[1] 3553




pvals<-c()
for(i in 1:nrow(d_prop_unique)){
	sab_wKmer<-d_prop_unique[i,1]
	ie_wKmer<-d_prop_unique[i,2]
	sab_woKmer<-121-sab_wKmer
	ie_woKmer<-120-ie_wKmer
	temp<-matrix(c(sab_wKmer,ie_wKmer,sab_woKmer,ie_woKmer),byrow=T, nrow=2)
	p<-fisher.test(temp)$p.value
	pvals<-c(pvals, p)
}

d_prop_unique$pval<-pvals

### merge with original set, corect for mulitple testing
d_prop_unique_test<-d_prop_unique[,c(3,4)]
m<-merge(d, d_prop_unique_test)
m$fdr<-p.adjust(m$pval, method="fdr")
min(m$fdr)
#[1] 1
min(m$pval)
#[1] 0.0004913731

ms<-m[sort.list(m$pval),]

        pattern                           kmer SAB  IE         pval fdr
1446311 102_117 ATAGTGAGAATCATTATCAATTAGGTAACA 102 117 0.0004913731   1
1762970  113_94 TCCATCTTCAAATTGATAAGGAACATTTTT 113  94 0.0008011367   1
1458403 103_117 AGTGAGAATCATTATCAATTAGGTAACACA 103 117 0.0009181122   1
1458404 103_117 TGATAGTGAGAATCATTATCAATTAGGTAA 103 117 0.0009181122   1
1458405 103_117 GATAGTGAGAATCATTATCAATTAGGTAAC 103 117 0.0009181122   1
1458406 103_117 TAGTGAGAATCATTATCAATTAGGTAACAC 103 117 0.0009181122   1
1458407 103_117 GTGAGAATCATTATCAATTAGGTAACACAC 103 117 0.0009181122   1
5163772    1_12 CCTTTTCGAAATTAATTTTGAAAACTCGTC   1  12 0.0013673628   1

>fasta1
ATAGTGAGAATCATTATCAATTAGGTAACA
>fasta2
AGTGAGAATCATTATCAATTAGGTAACACA
>fasta3
TGATAGTGAGAATCATTATCAATTAGGTAA
>fasta4
GATAGTGAGAATCATTATCAATTAGGTAAC
>fasta5
TAGTGAGAATCATTATCAATTAGGTAACAC
>fasta6
GTGAGAATCATTATCAATTAGGTAACACAC

>All_sign1
TGATAGTGAGAATCATTATCAATTAGGTAACACAC

fasta3          TGATAGTGAGAATCATTATCAATTAGGTAA-----
fasta1          --ATAGTGAGAATCATTATCAATTAGGTAACA---
fasta4          -GATAGTGAGAATCATTATCAATTAGGTAAC----
fasta6          -----GTGAGAATCATTATCAATTAGGTAACACAC
fasta2          ----AGTGAGAATCATTATCAATTAGGTAACACA-
fasta5          ---TAGTGAGAATCATTATCAATTAGGTAACAC--
                     *************************     



one of the samples,SAB_nonCA_CC30_6_NODE_13 matches fasta 6, fasta 2 and fasta 5
another sample, SAB_nonCA_CC45_NODE_17 matches fasta 3 and 4
