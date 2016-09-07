#PS5 RNA-seq Differential Expression Analysis 
#September 3, 2016 
#Using the raw read counts tsv file from an RNA-seq experiment



#Question 1
#set working directory 
getwd()
setwd("/Users/jsmith/Documents/bi623/160902_class15_PS5/")

#read in the tsv file. 
Counts <- read.table("Gacu_gut_counts.tsv", header = TRUE, sep = "\t",
                     row.names = 1) 

#Characterize the dataset
head(counts)
dim(counts)

#Use apply function to find sums of each library.
sum_counts <- apply(Counts, MARGIN = 2, FUN = sum )s

#Read in the file with the raw fastq read counts. 
countsFastq <- read.table("fastqReads.txt", header = FALSE, sep = ":")

#make a vector to hold only the read counts
fastqReads <- countsFastq[,2]

#Divide the aligned read counts by the fastq read counts
readPorportion <- sum_counts/fastqReads  

AveReadPorportion <- mean(readPorportion) * 100
#61.33% of raw reads aligned to the reference genome. 


#Question 2
#Evaluating the effect of normalization on the count data

#download the binary files for the packages 
source("https://bioconductor.org/biocLite.R")

#install the EdgeR package
biocLite("edgeR")

#load the package into R
library(edgeR)

#create a PDF file. 
pdf(file = "Normalized_Unnormalized_ReadCounts.pdf")
#Parameters for plots 
op <- par(mfrow = c(1,2), pty = "s")
#plot the unnormalized read counts for library 1A_02 values 
#vs. library 1A_03 values.
plot(Counts$X1A_02, Counts$X1A_03, col = "darkblue", pch = 19, cex=0.35, main = "Unnormalized Read Counts",
     xlim = c(0,50000), ylim = c(0,50000), ylab = "Read Counts Sample 1A_03",
     xlab = "Read Counts Sample 1A_02", cex.axis = 0.75)
#Draw a line through the origin (intercept=1) with slope=1. 
abline(a = 1, b = 1, col = "red", cex=0.8)

#store raw counts are a "DGEList" object. 
dge <-DGEList(counts = Counts)

#Calcultate TMM Normalization factors to dge.
#cpm = counts per million 
TMMCPM <- as.data.frame(cpm(dge, normalized.lib.sizes = TRUE))

#Plot the normalized read counts. 
plot(TMMCPM$X1A_02, TMMCPM$X1A_03, ylab = "Read Counts Sample 1A_03",
     xlab = "Read Counts Sample 1A_02", main = "Normalized Read Counts",
     col = "cornflowerblue", pch = 19, cex=0.35,  cex.axis = 0.75, 
     xlim = c(0,50000), ylim = c(0,50000))
abline(a=1, b=1, col = "red", cex=0.8)
#clsoe the pdf file.
dev.off()

#Question 3
#Testing and visualizing effects of Population and Treatment on 
#differential expression using the negative binomial generalized
#linear model approach in edgeR


#define a new dge list
dge.er <- DGEList(counts=Counts)

#Keep genes with at least two reads per million (cpm) in at least 8 samples. 
#rowSums is to calculate sums and means from rows and columns.
keep.dge.er <-rowSums(cpm(dge.er) >= 1) >8 
 
#index notation to subset the dge.er gene list for only those
#genes that have more than 2 counts per million
dge.er <- dge.er[keep.dge.er,]

#calculate the TMM normalization factors 
dge.er <- calcNormFactors(dge.er)

#read in the metadata file.  
Targets.er <- read.delim(file = "Gacu_gut_metadata.tsv", header =TRUE, sep = "\t")

#Set up the model by combining our population and treatment variablesinto one.
PopTrt <- paste(Targets.er$Population, Targets.er$Treatment, sep = ".")
#paste function concatentates vectors after converting to character. 
PopTrt

#now use the model.matrix() function to specify the experimental design. 
design.er <- model.matrix(~0 + PopTrt)
design.er

#This created four treatment groups.  
colnames(design.er)

#Estimate the dispersion parameter of each gene. 
#Estimate common, trended and tagwise dispersions for desing.er. 
#this is necessary for modeling the negative binomial distribution. 
dge.er <- estimateGLMCommonDisp(dge.er, design.er)
dge.er <- estimateGLMTrendedDisp(dge.er, design.er)
dge.er <- estimateGLMTagwiseDisp(dge.er, design.er)

#Fit the full generalized linear model, incorporating dispersion estimates. 
fit.er <- glmFit(dge.er, design.er)
fit.er

#Population Contrast 
#To test for the Population effect, contrast the average 
#of Bt.CV and Bt.GF vs. the average of RS.CV and RS.GF. 
lrtPop <- glmLRT(fit.er, contrast = c(-0.5, -0.5, 0.5, 0.5))

#topTags function extracts top DE genes in a data frame. 
topTags(lrtPop, n=16877)

#create a variable to hold the topTags objects.
topLrtPop<- topTags(lrtPop, n=16877)

#Create a subset of the topTags for those genes with FDR <= 0.05
subsetLRTpop <- subset(topLrtPop$table, topLrtPop$table$FDR <= 0.05) 

#check that the subset worked 
head(subsetLRTpop)
max(subsetLRTpop$FDR)

#write all of the test results to a file. 
write.table(topTags(lrtPop, n=16877), 'edgeR_PopLRT.txt', sep = '\t')

#Summarize how many genes are differentially expressed by population
summary(de.er <- decideTestsDGE(lrtPop, p=0.05))

#diffeneretiall expressed tags(reads) have the row names
#added from the original fit model (fit.er)  
detags.er <- rownames(fit.er)[as.logical(de.er)]

#create a pdf file.  
pdf(file = "Pop_GermTreatment_Contrast.pdf")
#set graphing parameters
op <- par(mfrow = c(1,2), pty ="s")
#visualize the differtial expression patterns with a "smear plot"
plotSmear(lrtPop, de.tags = detags.er,  xlim = c(0,14), ylim = c(-5,5), 
          main = "Differentially Expressed Genes \n Between Boot Lake and Rabbit Slough \n Stickleback Populations", 
          cex=0.3, pch = 19, cex.main = 1.0, cex.axis = 0.8)
#abline h option is the Y-values for horizontal lines
abline(h = c(-2,2), col = "blue")

#Treatment Contrast 
#to test for the treatment (microbiota) effects we contrast the
#average of Bt.CV and RS.CV the average of Bt.GF and RS.GF

lrtTrt <- glmLRT(fit.er, contrast = c(-.5, .5, -.5, .5))

#look at the genes with the lowest FDR corrected p-values. 
topTags(lrtTrt)

#create a variable to hold the topTag object
topLrtTrt <- topTags(lrtTrt, n=16877)
dim(topLrtTrt)

#subset the topTag object to those genes with FDR <= 0.05
subsetLrtTrt <- subset(topLrtTrt$table, subset = topLrtTrt$table$FDR <= 0.05)

#write the results to a file.
write.table(topTags(lrtTrt, n=16877), 'edgeR_TrtLRT.txt', 
            sep = '\t')

#Create a smear plot for the Treatment comparison
summary(deTrt.er <- decideTestsDGE(lrtTrt, p=0.05))

#pull out the row names from the original model fit
detagsTrt.er <- (rownames(fit.er) [as.logical(deTrt.er)]) 

#Create a Plot Smear graph of the differentially expressed genes between
#conventional and germ-free treatments
plotSmear(lrtTrt, de.tags = detagsTrt.er, xlim = c(0,14), ylim = c(-5,5), 
          main = "Differentially Expressed Genes \n between Convential and Germ-Free \n Stickleback Treatments", 
          cex=0.3, pch = 19, cex.main = 1.0, cex.axis = 0.8, mar = c(5,5,4,2) + 0.1)
abline(h = c(-2,2), col = "blue")
dev.off()


#Question 4 
#Testing effects of Population and Treatment on differential expression 
#using the general linear model approach in limma (with voom())

#create a DGE list
dge.lv <- DGEList(counts = Counts)
#keep only genes with more than 2 counts per million in 8 or more samples.
keep.dge.lv <- rowSums(cpm(dge.lv) >= 1) >= 8 
#use index notation to subset the original DGE list with genes to keep.
dge.lv <- dge.lv[keep.dge.lv,]

#calcualte TMM Normalization 
dge.lv <- calcNormFactors(dge.lv)

#read in the metadata 
Targets.lv <- read.delim('Gacu_gut_metadata.tsv', header = TRUE, sep = '\t')
#define contrast groups 
LVPopTrt <- paste(Targets.lv$Population, Targets.lv$Treatment, sep = ".") 
#use model matrix to specify experimental design
desing.lv <- model.matrix(~0 + LVPopTrt)

#print the column names of desing.lv 
colnames(desing.lv)

#use Voom to generate "precision weights" from the mean-variance relationships
v.lv <- voom(dge.lv, desing.lv)

#fit the full general linear model, which includes the voom() precision weights
fit.lv <- lmFit(v.lv, desing.lv)

#compute Bayes statistics for Differential Expression 
fit.lv <- eBayes(fit.lv)


#Population Contrast
#compute contrasts from linear model fit
fit.lv.pop <- contrasts.fit(fit.lv, 
                            contrast = c(-.5, -.5, .5, .5))
#this function contrasts.fit re-orients the fitted model object from the 
#coefficients of the original design.

#Ebayes statistics for Differential Expression 
fit.lv.pop <- eBayes(fit.lv.pop)
fit.lv.pop

#look at all the genes with the lowest FDR-corrected p-values
topTable(fit.lv.pop)

#Write results table to a file.
write.table(topTable(fit.lv.pop, n=16877), 'limma_popLRT.txt', sep = '\t')

#summarize how many genes are differentially expressed by Populaton 
#using an FDR = 0.05. 
subsetPop <- subset(topTable(fit.lv.pop, n=16877), subset = adj.P.Val <= 0.05)

#check that the subset was correct 
head(subsetFDR)
max(subsetFDR$adj.P.Val)

#Treatment Contrast 
#To test for the Treatment effect we contrast the average of Bt.CV and RS.CV vs.
#the average of Bt.GF and RS.GF
fit.lv.Trt <- contrasts.fit(fit.lv, contrast = c(-.5,.5,-.5,.5))
fit.lv.Trt<- eBayes(fit.lv.Trt)

#look at genes with the lowest FDR-corrected p-values
topTable(fit.lv.Trt)

#write resutls to a file.
write.table(topTable(fit.lv.Trt, n=16877), 'limma_TrtLRT.txt', sep ='\t')

#Summarize how many genes are differentiall expressed by treatment,
#using an FDR = 0.05 and subset and parse topTable. 
subsetTrt <- subset(topTable(fit.lv.Trt, n=16877), subset = adj.P.Val <= 0.05)

#check that subset worked
head(subsetTrt)
#no genes fall below 0.05 FDR 

#Question 5
#Compare the overlap in differentially expressed genes between 
#edgeR and limma-voom

#install package for Venn Diagram
install.packages('VennDiagram')

#Load the library 
library(VennDiagram)

#Make a Venn Diagram to show overlap between EdgeR and Limma approaches.
#Show overlap of significant genes (FDR = 0.05)

#Variables to use Limma:
#subsetPop, subsetTrt 
namesLVpop <- rownames(subsetPop)
namesLVtrt <- rownames(subsetTrt) 

#Variable to use EdgeR:
#subsetLRTpop, subsetLrtTrt
namesLrtPop <- rownames(subsetLRTpop)
namesLrtTrt <- rownames(subsetLrtTrt)


#Characterize the population test.
#find the set of significant DE genes that overlaps between edgeR and limmaVoom approaches
vennPop <- venn.diagram(
  list(Limma = namesLVpop, EdgeR = namesLrtPop), 
  "VennDiagramPopulation.tiff", fill = c("dodgerblue3", "magenta4"), 
  main = "Comparison of Differentially Expressed Genes \n Identified Between Populations" )


#Characterize the treatment test. 
#find the set of significant DE genes between EdgeR and Limma

vennTrt <- venn.diagram(list(Limma = namesLVtrt, EdgeR = namesLrtTrt),
                        "VennDiagramTreatment.tiff", fill = c("dodgerblue3", "magenta4"),
                        main = "Comparison of Differentially Expressed Genes \n Identified Between Treatments")






