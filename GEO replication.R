####################################################################################
# Replication analysis for EWAS of infant sex in cord blood; GEO dataset GSE129841 #
####################################################################################

# This R script accompanies the manuscript "Associations between infant sex and DNA methylation across umbilical cord blood, artery, and placenta samples" and provides code for replication analysis for EWAS of infant sex in cord blood using the GEO dataset GSE129841, including downloading GEO data, performing analysis of differentially methylated positinos, and performing analysis of differentially methylated regions.

# Author: Anne Bozack
# Date: 9/16/2020

# load required libraries

library(GEOquery)
library(Harman)
library(DMRcate)
library(limma)
library(ENmix)

# set working directory

path <- '/Users/annebozack/Documents/GEO_analysis'
setwd(path)


################
# loading data #
################

# download GEO object

geoMat <- getGEO("GSE129841")


# extract pheno data and betas values

GRset <- getGenomicRatioSetFromGEO("GSE129841")

pheno <- pData(GRset)
beta <- getBeta(GRset)

dim(pheno)
# 114  41

dim(beta)
# 434506  114


##########################################
# removing cross-reactive and SNP probes #
##########################################

# removing additional non-specific probes (Chen et al. Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray. Epigenetics 2013, 8:203ñ9.)

chen <- read.csv('48639-non-specific-probes-Illumina450k.csv') 

beta <- beta[!(rownames(beta) %in% chen$TargetID),]

dim(beta)
# 409741    114


# removing additional probes <= 10 base pairs from SNPs with minor allele frequencies >= 0.05 (Okamura et al. Lists of HumanMethylation450 BeadChip probes with nucleotide-variant information obtained from the Phase 3 data of the 1000 Genomes Project. Genom Data 2016, 7:67–69.)

probesDrop <- read.csv('SNPprobes_le10bp_ge5MAF_Okamura.csv')[,2]

beta <- beta[!(rownames(beta) %in% probesDrop),]

dim(beta)
# 390729    114


############
# ReFactor #
############

# performing ReFactor to estimate components representing cell type proportions (https://github.com/cozygene/refactor; Rahmani et al. Sparse PCA corrects for cell type heterogeneity in epigenome-wide association studies. Nature Methods 2016, 13:443–445)

# saving beta values in tsp format

betaRefact <- cbind(rownames(beta), data.frame(beta))

colnames(betaRefact)[1] <- 'ID'

write.table(betaRefact, file <- 'GEO_beta_refactor.txt', sep = "\t", col.names = TRUE, row.names = FALSE)


# run refactor in python using 6 comonents

# cd /Users/annebozack/Documents/GEO_analysis

# python refactor.py --datafile GEO_beta_refactor.txt --k 6 --out refactor_GEO


# adding refactor components to pheno dataframe

refactComp <- as.matrix(read.table('refactor_GEO.out.components.txt', sep="\t", header=FALSE))

refactComp <- data.frame(refactComp)
colnames(refactComp) <- c('refact1', 'refact2', 'refact3', 'refact4', 'refact5', 'refact6')

pheno <- cbind(pheno, refactComp)

table(colnames(beta) %in% rownames(pheno))
# TRUE 
# 114 


######################################
# converting Beta-values to M-values #
######################################

table(beta < 0)
# FALSE     TRUE 
# 44542912   194

table(beta < 1)
# FALSE     TRUE 
# 6 44543100 

beta[beta < 0] <- 0
beta[beta > 1] <- 1

beta <- shiftBetas(beta, shiftBy=1e-4)

GEO_M <- B2M(beta)

# Saving M values and pheno dataframe
save(GEO_M, pheno, file  = 'GEO_GSE129841_filteredM_pheno.Rdata')


#################
# DMRs: DMRcate #
#################

# performing DMRcate adjusted for refactor components

colnames(pheno)[39] <- 'sex'
pheno$sex_dummy[pheno$sex == 'male'] <- 0
pheno$sex_dummy[pheno$sex == 'female'] <- 1
pheno$sex_dummy <- factor(pheno$sex_dummy, levels = c(0,1))

design <- model.matrix(~ sex_dummy + refact1 + refact2 + refact3 + refact4 + refact5 + refact6, data = pheno)

myannotation <- cpg.annotate("array", GEO_M, what="M", arraytype = "450K", analysis.type="differential", design=design, coef=2, fdr = 0.01)

dmroutput <- dmrcate(myannotation, lambda = 1000, C = 2, min.cpgs = 10)

GEO.results.ranges <- extractRanges(dmroutput, genome = 'hg19')


# regions with mean methylation higher among females

GEO.results.ranges <- data.frame(GEO.results.ranges)
table(GEO.results.ranges$meandiff > 0)
# FALSE  TRUE 
#   3    44 


# regions annotated to SNOR genes

GEO.results.ranges$SNOR <- NA
for (i in 1:nrow(GEO.results.ranges)){
	if (grepl('SNOR', GEO.results.ranges$overlapping.genes[i])) {
		GEO.results.ranges$SNOR[i] <- 1
	} else {
		GEO.results.ranges$SNOR[i] <- 0
	}
}

table(GEO.results.ranges$SNOR)
#   0   1 
#  16  31 


###############
# DMPs: limma #
###############

# performing limma adjusted for refactor components

var <- model.matrix(~ sex_dummy + refact1 + refact2 + refact3 + refact4 + refact5 + refact6, data = pheno)

fit <- lmFit(GEO_M, var)
fit <- eBayes(fit)

geoProbeSexAadj <- topTable(fit, adjust = 'BH', coef = 2, number = Inf, confint = T)
geoProbeSexAadj$adj.P.Val.bonf <- topTable(fit,adjust="bonferroni",coef=2, number = Inf)$adj.P.Val

table(geoProbeSexAadj$adj.P.Val < 0.05)
# FALSE   TRUE 
# 386555   4174 

table(geoProbeSexAadj$adj.P.Val.bonf < 0.05)
#  FALSE   TRUE 
# 390215    514   


# probes with methylation higher among females

table(geoProbeSexAadj$logFC[geoProbeSexAadj$adj.P.Val.bonf < 0.05] > 0)

