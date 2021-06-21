#Statistical tests for comparison of DEGs across different life stage transitions and comparisons among different gene categories (chemosensory and RPL)

#Analysis of the proportion and magnitude of stage bias across the four developmental/sex comparisons (Figure 2)
#DEGs are padj<0.05 and FC > |1| 

#Descriptive statistics for the proportion of DEGs in different comparisons
#Clopper-Pearson 95% confidence intervals for proportion of genes that were differentially expressed for each comparison (3 developmental, 1 sex)
#install.packages('DescTools')
library(DescTools) # version 0.99.40; needed for binomial confidence intervals

#all tissues
minorPropAll_FC<- BinomCI(114,114+6905, conf.level = 0.95, method="clopper-pearson")
majorPropAll_FC <- BinomCI(1455,1455+5807, conf.level = 0.95, method="clopper-pearson")
completePropAll_FC <- BinomCI(1553,1553+6385,conf.level=0.95, method="clopper-pearson")
sexPropAll_FC <- BinomCI(369,369+7774,conf.level=0.95,method="clopper-pearson")

PropAll_FC <- as.data.frame(rbind(minorPropAll_FC,majorPropAll_FC,completePropAll_FC,sexPropAll_FC))
PropAll_FC$Transition <- c('minor', 'major', 'complete','sex')
write.csv(PropAll_FC, quote = FALSE, file="propDEG_withCI_FC_Alltissues.csv") 

#heads only
minorPropHead_FC<- BinomCI(142,142+6266, conf.level = 0.95, method="clopper-pearson")
majorPropHead_FC <- BinomCI(850,850+5413, conf.level = 0.95, method="clopper-pearson")
completePropHead_FC <- BinomCI(1710,1710+4847,conf.level=0.95, method="clopper-pearson")
sexPropHead_FC <- BinomCI(264,264+6362,conf.level=0.95,method="clopper-pearson")

PropHead_FC <- as.data.frame(rbind(minorPropHead_FC,majorPropHead_FC,completePropHead_FC,sexPropHead_FC))
PropHead_FC$Transition <- c('minor', 'major', 'complete','sex')
write.csv(PropHead_FC, quote = FALSE, file="propDEG_withCI_FC_HeadsOnly.csv") 

#Fisher's exact tests to determine if proportions differ significantly among comparisons
#install.packages('RVAideMemoire')
library("RVAideMemoire") # version 0.9-78; needed for post-hoc pairwise tests

#all tissues
DEGallFC<-as.table(matrix(c(114,1455,1553,369,6905,5807,6385,7774),ncol=2,dimnames=list(c("minor","major", "complete","sex"),c("DEG","notDEG"))))
DEGallFC #print table to check numbers

#global fisher's exact test
fisher.test(DEGallFC,simulate.p.value=TRUE) #pvalue must be simulated

#posthoc pairwise comparisons
fisher.multcomp(DEGallFC) # p adjusted with 'fdr' method, which is the same as the Benjamini & Hochberg method

#heads only
DEGheadFC<-as.table(matrix(c(142,850,1710,264,6266,5413,4847,6362),ncol=2,dimnames=list(c("minor","major", "complete","sex"),c("DEG","notDEG"))))
DEGheadFC #print table to check values

#global fisher's exact test
fisher.test(DEGheadFC,simulate.p.value=TRUE) 

#pairwise comparisons
fisher.multcomp(DEGheadFC) 

#Descriptive statistics for the magnitude of DEGs in different comparisons

#all tissues
AllFC<-read.delim(file="DEGs_AllTissues_absFC_minFC1.txt",header=TRUE)

#median and interquartile range for each comparison
median(AllFC$AbsFC[AllFC$comp=="aMinor"])
quantile(AllFC$AbsFC[AllFC$comp=="aMinor"])
median(AllFC$AbsFC[AllFC$comp=="bMajor"])
quantile(AllFC$AbsFC[AllFC$comp=="bMajor"])
median(AllFC$AbsFC[AllFC$comp=="Complete"])
quantile(AllFC$AbsFC[AllFC$comp=="Complete"])
median(AllFC$AbsFC[AllFC$comp=="Sex"])
quantile(AllFC$AbsFC[AllFC$comp=="Sex"])

#Non-parametric rank-sum tests to determine whether magnitude of stage-bias differs among comparisons
#install.packages("FSA")
library(FSA) #v0.8.32- needed for post-hoc Dunn tests

#Kruskal-Wallis rank sum tests to determine whether there are differences among groups
kruskal.test(AllFC$AbsFC~AllFC$comp) 

#posthoc pairwise tests
dunnTest(AllFC$AbsFC~AllFC$comp, method="bh")

#heads only
HeadFC<-read.delim(file="DEGs_HeadOnly_absFC_minFC1.txt", header=TRUE)

#median and interquartile range for each comparison
median(HeadFC$AbsFC[HeadFC$comp=="aMinor"])
quantile(HeadFC$AbsFC[HeadFC$comp=="aMinor"])
median(HeadFC$AbsFC[HeadFC$comp=="bMajor"])
quantile(HeadFC$AbsFC[HeadFC$comp=="bMajor"])
median(HeadFC$AbsFC[HeadFC$comp=="Complete"])
quantile(HeadFC$AbsFC[HeadFC$comp=="Complete"])
median(HeadFC$AbsFC[HeadFC$comp=="Sex"])
quantile(HeadFC$AbsFC[HeadFC$comp=="Sex"])


#Kruskal-Wallis rank sum tests to determine whether there are differences among groups
kruskal.test(HeadFC$AbsFC~HeadFC$comp)

#posthoc pairwise tests
dunnTest(HeadFC$AbsFC~HeadFC$comp, method="bh")

#Comparing absFC between chemosensory genes and RPLs across all transitions and within each type (Figure 4)
#use nonparametric wilcoxon rank stu
gene<-read.csv(file="ChemoVsRPL.csv", header=TRUE)

#global test 
wilcox.test(gene$absFC[gene$GeneType=="chemo"],gene$absFC[gene$GeneType=="RPL"])

#subset data to look at each transition independently
minor<-subset(gene,gene$Transition=="aMinor")
major<-subset(gene,gene$Transition=="bMajor")
complete<-subset(gene,gene$Transition=="Complete")
sex<-subset(gene,gene$Transition=="Sex")

#minor
wilcox.test(minor$absFC[minor$GeneType=="chemo"],minor$absFC[minor$GeneType=="RPL"])

#major
wilcox.test(major$absFC[major$GeneType=="chemo"],major$absFC[major$GeneType=="RPL"])

#complete
wilcox.test(complete$absFC[complete$GeneType=="chemo"],complete$absFC[complete$GeneType=="RPL"])

#sex
wilcox.test(sex$absFC[sex$GeneType=="chemo"],sex$absFC[sex$GeneType=="RPL"])

