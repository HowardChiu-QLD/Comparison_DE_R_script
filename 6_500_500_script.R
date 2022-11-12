#### IFN646 Biomedical Data Science - Group Project (Group 8)
### Gene-level differential expression analysis using DESeq2, edgeR and limma-voom
### This script is for 6_500_500

## Load libraries
library(DESeq2)
library(limma)
library(edgeR)
## gplots is used to plot a Venn Diagram of overlapping genes between the three different tools
# install.packages('gplots')
library(gplots)

## Load data and the metadata file
data <- read.table("data/6_500_500.tsv", header=T, row.names=1)
meta <- read.table("meta/6_500_500_meta.txt", header=T, row.names=1)

## Check that sample names match in both files
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

## Check datatype
class(data)
class(meta)
str(meta)
## Change "sampletype" to factor and check it again
meta$sampletype <- as.factor(meta$sampletype) 
str(meta)

##### DESeq2 #####
## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
## Run analysis
dds <- DESeq(dds)

## Set threshold
p.threshold <- 0.05

## Identify a list of genes showing differential expression by DESeq2
resultsNames(dds)
contrast.deseq2 <- list("sampletype_ex_vs_control")
deseq2_results <- results(dds, contrast=contrast.deseq2)
deseq2_results$threshold <- as.logical(deseq2_results$padj < p.threshold)
genes.deseq2 <- row.names(deseq2_results)[which(deseq2_results$threshold)]
## Check the results (genes showing DE by DESeq2)
length(genes.deseq2)
genes.deseq2

## Create g1 to g1000
genes1000 <- row.names(data)[which((nchar(row.names(data)) < 5)|(row.names(data) == 'g1000'))]
genes1000

## Find those between g1 and g1000 with DE
genes1000.deseq2 <- genes.deseq2[which((nchar(genes.deseq2) < 5)|(genes.deseq2 == 'g1000'))]
genes1000.deseq2

## Find those between g1 and g1000 without DE
genes1000.deseq2_noDE <- genes1000[!genes1000 %in% genes1000.deseq2]
genes1000.deseq2_noDE

## Create g1001 to g9987
genes_outside <- row.names(data)[!row.names(data) %in% genes1000]

## Find those between g1001 and g9987 with DE
genes_outside.deseq2 <- genes.deseq2[!genes.deseq2 %in% genes1000.deseq2]
genes_outside.deseq2

## Find those between g1001 and g9987 without DE
genes_outside.deseq2_noDE <- genes_outside[!genes_outside %in% genes_outside.deseq2]

## Calculate TP (genes between g1 and g1000 showing differential expression by DESeq2)
TP_deseq2 <- length(genes1000.deseq2)
TP_deseq2
## Calculate FP (genes between g1001 and g9987 showing differential expression by DESeq2)
FP_deseq2 <- length(genes_outside.deseq2)
FP_deseq2 
## Calculate TN (genes between g1001 and g9987 without showing differential expression by DESeq2)
TN_deseq2 <- length(genes_outside.deseq2_noDE)
TN_deseq2
## Calculate FN (genes between g1 and g1000 without showing differential expression by DESeq2)
FN_deseq2 <- length(genes1000.deseq2_noDE)
FN_deseq2

## Calculate True Positive Rate (TPR = TP/(TP+FN)) (= Recall)
TPR_deseq2 <- TP_deseq2/(TP_deseq2+FN_deseq2)
TPR_deseq2  ### TPR= recall
## Calculate False Positive Rate (FPR = FP/(TN+FP))
FPR_deseq2 <- FP_deseq2/(TN_deseq2+FP_deseq2)
FPR_deseq2
## Calculate Accuracy (= (TP+TN)/(TP+FP+FN+TN))
accuracy_deseq2 <- (TP_deseq2+TN_deseq2)/(TP_deseq2+FP_deseq2+FN_deseq2+TN_deseq2)
accuracy_deseq2
## Calculate Precision (= TP/(TP+FP))
precision_deseq2 <- TP_deseq2/(TP_deseq2+FP_deseq2)
precision_deseq2
## Calculate AUC score (= 1/2-(FPR/2)+(TPR/2))
auc_deseq2 <- 1/2-(FPR_deseq2/2)+(TPR_deseq2/2)
auc_deseq2


##### EdgeR #####
dge <- DGEList(counts = data, group = meta$sampletype)
## Normalize by total count
dge <- calcNormFactors(dge)

## Design matrix
design.mat <- model.matrix(~ 0 + dge$samples$group)
colnames(design.mat) <- levels(dge$samples$group)
design.mat

## Estimate dispersion parameter for GLM
dge <- estimateGLMCommonDisp(dge, design.mat)
dge <- estimateGLMTrendedDisp(dge, design.mat, method="power")
dge <- estimateGLMTagwiseDisp(dge,design.mat)

## Model fitting
fit.edgeR <- glmFit(dge, design.mat)
## Identify a list of genes showing differential expression by EdgeR
# Differential expression #### ex - control
contrasts.edgeR <- makeContrasts(ex - control, levels=design.mat)
lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)
# Access results tables
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]
## Check the results (genes showing DE by DESeq2)
length(genes.edgeR)
genes.edgeR

## Find those between g1 and g1000 with DE
genes1000.edgeR <- genes.edgeR[which((nchar(genes.edgeR) < 5)|(genes.edgeR == 'g1000'))]
genes1000.edgeR

## Find those between g1 and g1000 without DE
genes1000.edgeR_noDE <- genes1000[!genes1000 %in% genes1000.edgeR]
genes1000.edgeR_noDE

## Find those between g1001 and g9987 with DE
genes_outside.edgeR <- genes.edgeR[!genes.edgeR %in% genes1000.edgeR]
genes_outside.edgeR

## Find those between g1001 and g9987 without DE
genes_outside.edgeR_noDE <- genes_outside[!genes_outside %in% genes_outside.edgeR]

## Calculate TP (genes between g1 and g1000 showing differential expression by EdgeR)
TP_edgeR <- length(genes1000.edgeR)
TP_edgeR
## Calculate FP (genes between g1001 and g9987 showing differential expression by EdgeR)
FP_edgeR <- length(genes_outside.edgeR)
FP_edgeR 
## Calculate TN (genes between g1001 and g9987 without showing differential expression by EdgeR)
TN_edgeR <- length(genes_outside.edgeR_noDE)
TN_edgeR
## Calculate FN (genes between g1 and g1000 without showing differential expression by EdgeR)
FN_edgeR <- length(genes1000.edgeR_noDE)
FN_edgeR

## Calculate True Positive Rate (TPR = TP/(TP+FN)) (= Recall)
TPR_edgeR <- TP_edgeR/(TP_edgeR+FN_edgeR)
TPR_edgeR  ### TPR= recall
## Calculate False Positive Rate (FPR = FP/(TN+FP))
FPR_edgeR <- FP_edgeR/(TN_edgeR+FP_edgeR)
FPR_edgeR
## Calculate Accuracy (= (TP+TN)/(TP+FP+FN+TN))
accuracy_edgeR <- (TP_edgeR+TN_edgeR)/(TP_edgeR+FP_edgeR+FN_edgeR+TN_edgeR)
accuracy_edgeR
## Calculate Precision (= TP/(TP+FP))
precision_edgeR <- TP_edgeR/(TP_edgeR+FP_edgeR)
precision_edgeR
## Calculate AUC score (= 1/2-(FPR/2)+(TPR/2))
auc_edgeR <- 1/2-(FPR_edgeR/2)+(TPR_edgeR/2)
auc_edgeR


##### limma-voom #####
# Create design matrix
design <- model.matrix(~ meta$sampletype)

# Apply voom transformation
nf <- calcNormFactors(data)
v <- voom(data, design, lib.size=colSums(data)*nf, 
          normalize.method="quantile", plot=FALSE)

# Usual limma pipeline
fit.voom <- lmFit(v, design)
fit.voom <- eBayes(fit.voom)
voom_results <- topTable(fit.voom, coef=2, adjust="BH", number = nrow(data))
voom_results$threshold <- as.logical(voom_results$adj.P.Val < p.threshold)
genes.limmavoom <- row.names(voom_results)[which(voom_results$threshold)]
length(genes.limmavoom)

## Find those between g1 and g1000 with DE
genes1000.limmavoom <- genes.limmavoom[which((nchar(genes.limmavoom) < 5)|(genes.limmavoom == 'g1000'))]
genes1000.limmavoom

## Find those between g1 and g1000 without DE
genes1000.limmavoom_noDE <- genes1000[!genes1000 %in% genes1000.limmavoom]
genes1000.limmavoom_noDE

## Find those between g1001 and g9987 with DE
genes_outside.limmavoom <- genes.limmavoom[!genes.limmavoom %in% genes1000.limmavoom]
genes_outside.limmavoom

## Find those between g1001 and g9987 without DE
genes_outside.limmavoom_noDE <- genes_outside[!genes_outside %in% genes_outside.limmavoom]

## Calculate TP (genes between g1 and g1000 showing differential expression by limma-voom)
TP_limmavoom <- length(genes1000.limmavoom)
TP_limmavoom
## Calculate FP (genes between g1001 and g9987 showing differential expression by limma-voom)
FP_limmavoom <- length(genes_outside.limmavoom)
FP_limmavoom 
## Calculate TN (genes between g1001 and g9987 without showing differential expression by limma-voom)
TN_limmavoom <- length(genes_outside.limmavoom_noDE)
TN_limmavoom
## Calculate FN (genes between g1 and g1000 without showing differential expression by limma-voom)
FN_limmavoom <- length(genes1000.limmavoom_noDE)
FN_limmavoom

## Calculate True Positive Rate (TPR = TP/(TP+FN)) (= Recall)
TPR_limmavoom <- TP_limmavoom/(TP_limmavoom+FN_limmavoom)
TPR_limmavoom  ### TPR= recall
## Calculate False Positive Rate (FPR = FP/(TN+FP))
FPR_limmavoom <- FP_limmavoom/(TN_limmavoom+FP_limmavoom)
FPR_limmavoom
## Calculate Accuracy (= (TP+TN)/(TP+FP+FN+TN))
accuracy_limmavoom <- (TP_limmavoom+TN_limmavoom)/(TP_limmavoom+FP_limmavoom+FN_limmavoom+TN_limmavoom)
accuracy_limmavoom
## Calculate Precision (= TP/(TP+FP))
precision_limmavoom <- TP_limmavoom/(TP_limmavoom+FP_limmavoom)
precision_limmavoom
## Calculate AUC score (= 1/2-(FPR/2)+(TPR/2))
auc_limmavoom <- 1/2-(FPR_limmavoom/2)+(TPR_limmavoom/2)
auc_limmavoom


##### Comparison #####
## Create tables to show the results of the three different tools
# Set the option to avoid representing as scientific notation
options(scipen=999)
# Compare TP, FP, TN, FN
de <- matrix(c(TP_deseq2,FP_deseq2,TN_deseq2,FN_deseq2,TP_edgeR,FP_edgeR,TN_edgeR,FN_edgeR,TP_limmavoom,FP_limmavoom,TN_limmavoom,FN_limmavoom), ncol=4, byrow=TRUE)
rownames(de) <- c("DESeq2","EdgeR","limma-voom")
colnames(de) <- c("TP","FP","TN", "FN")
de <- as.table(de)
de
# Compare TPR(Recall), FPR, Accuracy, Precision, AUC
de_eval <- matrix(c(TPR_deseq2,FPR_deseq2,accuracy_deseq2,precision_deseq2,auc_deseq2,TPR_edgeR,FPR_edgeR,accuracy_edgeR,precision_edgeR,auc_edgeR,TPR_limmavoom,FPR_limmavoom,accuracy_limmavoom,precision_limmavoom,auc_limmavoom), ncol=5, byrow=TRUE)
rownames(de_eval) <- c("DESeq2","EdgeR","limma-voom")
colnames(de_eval) <- c("TPR(Recall)", "FPR", "Accuracy", "Precision", "AUC score")
de_eval <- as.table(de_eval)
de_eval

## Plot a Venn Diagram of overlapping genes between the three different tools
venn(list(edgeR = genes.edgeR, DESeq2 = genes.deseq2, limmavoom = genes.limmavoom))
bythree <- intersect(intersect(genes.edgeR,genes.deseq2),genes.limmavoom)
bythree
length(bythree)


## Plot a Venn Diagram of overlapping genes in the golden set)
venn(list(edgeR = genes1000.edgeR, DESeq2 = genes1000.deseq2, limmavoom = genes1000.limmavoom))
## Plot a Venn Diagram of overlapping genes NOT in the golden set)
venn(list(edgeR = genes_outside.edgeR, DESeq2 = genes_outside.deseq2, limmavoom = genes_outside.limmavoom))
# identified by all of these three tools and is in the golden set)
bythree1000 <- intersect(intersect(genes1000.edgeR,genes1000.deseq2),genes1000.limmavoom)
bythree1000
# is NOT in the golden set and is identified by all of these three tools)
bythree_outside <- intersect(intersect(genes_outside.edgeR,genes_outside.deseq2),genes_outside.limmavoom)
bythree_outside
# is in the golden set and is NOT identified by all of these three tools
notByAll<- genes1000[!genes1000 %in% intersect(intersect(genes1000.edgeR,genes1000.deseq2),genes1000.limmavoom)]
notByAll
# TP by three tools
TP_bythree <- length(bythree1000)
TP_bythree
# FP by three tools
FP_bythree <- length(bythree_outside)
FP_bythree
# FN by three tools
FN_bythree <- length(notByAll)
FN_bythree
# TPR by three tools
TPR_bythree <- TP_bythree/(TP_bythree+FN_bythree)
TPR_bythree
# Precision by three tools
precision_bythree <- TP_bythree/(TP_bythree+FP_bythree)
precision_bythree

###### Try two tools: D&E (DnE) (DESeq2 + EdgeR)
# intersection and is in the golden set
EandD1000 <- intersect(genes1000.edgeR,genes1000.deseq2)
EandD1000
length(EandD1000)
# intersection and is NOT in the golden set
EandD_outside <- intersect(genes_outside.edgeR,genes_outside.deseq2)
EandD_outside
length(EandD_outside)
# union and is in the golden set
EorD <- union(genes1000.edgeR,genes1000.deseq2)
EorD
length(EorD)
# E-D (found by EdgeR and not by DESeq2) and is in the golden set
EnotD <- setdiff(genes1000.edgeR,genes1000.deseq2)
EnotD
length(EnotD)
# D-E (found by DESeq2 and not by EdgeR) and is in the golden set
DnotE <- setdiff(genes1000.deseq2,genes1000.edgeR)
DnotE
length(DnotE)
# not found by DESeq2 or not by EdgeR and is in the golden set
notE_or_notD <- genes1000[!genes1000 %in% intersect(genes1000.edgeR,genes1000.deseq2)]
notE_or_notD
length(notE_or_notD)
# TP by D&E
TP_byDnE <- length(EandD1000)
TP_byDnE
# FP by D&E
FP_byDnE <- length(EandD_outside)
FP_byDnE
# FN by D&E
FN_byDnE <- length(notE_or_notD)
FN_byDnE
# TPR by D&E
TPR_byDnE <- TP_byDnE/(TP_byDnE+FN_byDnE)
TPR_byDnE
# Precision by D&E
precision_byDnE <- TP_byDnE/(TP_byDnE+FP_byDnE)
precision_byDnE

###### Try two tools: D&L (DnL) (DESeq2 + Limmavoom)
# intersection and is in the golden set
LandD1000 <- intersect(genes1000.limmavoom,genes1000.deseq2)
LandD1000
length(LandD1000)
# intersection and is NOT in the golden set
LandD_outside <- intersect(genes_outside.limmavoom,genes_outside.deseq2)
LandD_outside
length(LandD_outside)
# union and is in the golden set
LorD <- union(genes1000.limmavoom,genes1000.deseq2)
LorD
length(LorD)
# L-D (found by Limmavoom and not by DESeq2) and is in the golden set
LnotD <- setdiff(genes1000.limmavoom,genes1000.deseq2)
LnotD
length(LnotD)
# D-L (found by DESeq2 and not by Limmavoom) and is in the golden set
DnotL <- setdiff(genes1000.deseq2,genes1000.limmavoom)
DnotL
length(DnotL)
# not found by DESeq2 or not by Limmavoom and is in the golden set
notL_or_notD <- genes1000[!genes1000 %in% intersect(genes1000.limmavoom,genes1000.deseq2)]
notL_or_notD
length(notL_or_notD)
# TP by D&L
TP_byDnL <- length(LandD1000)
TP_byDnL
# FP by D&L
FP_byDnL <- length(LandD_outside)
FP_byDnL
# FN by D&L
FN_byDnL <- length(notL_or_notD)
FN_byDnL
# TPR by D&L
TPR_byDnL <- TP_byDnL/(TP_byDnL+FN_byDnL)
TPR_byDnL
# Precision by D&L
precision_byDnL <- TP_byDnL/(TP_byDnL+FP_byDnL)
precision_byDnL


###### Try two tools: L&E (LnE) (Limmavoom + EdgeR)
# intersection and is in the golden set
EandL1000 <- intersect(genes1000.edgeR,genes1000.limmavoom)
EandL1000
length(EandL1000)
# intersection and is NOT in the golden set
EandL_outside <- intersect(genes_outside.edgeR,genes1000.limmavoom)
EandL_outside
length(EandL_outside)
# union and is in the golden set
EorL <- union(genes1000.edgeR,genes1000.limmavoom)
EorL
length(EorL)
# E-L (found by EdgeR and not by Limmavoom) and is in the golden set
EnotL <- setdiff(genes1000.edgeR,genes1000.limmavoom)
EnotL
length(EnotL)
# D-E (found by Limmavoom and not by EdgeR) and is in the golden set
LnotE <- setdiff(genes1000.limmavoom,genes1000.edgeR)
LnotE
length(LnotE)
# not found by Limmavoom or not by EdgeR and is in the golden set
notE_or_notL <- genes1000[!genes1000 %in% intersect(genes1000.edgeR,genes1000.limmavoom)]
notE_or_notL
length(notE_or_notL)
# TP by L&E
TP_byLnE <- length(EandL1000)
TP_byLnE
# FP by L&E
FP_byLnE <- length(EandL_outside)
FP_byLnE
# FN by L&E
FN_byLnE <- length(notE_or_notL)
FN_byLnE
# TPR by L&E
TPR_byLnE <- TP_byLnE/(TP_byLnE+FN_byLnE)
TPR_byLnE
# Precision by L&E
precision_byLnE <- TP_byLnE/(TP_byLnE+FP_byLnE)
precision_byLnE

# Compare TPR and Precision between every two-tool pair and the three-tool result
tprp <- matrix(c(TPR_byDnE,precision_byDnE,TPR_byDnL,precision_byDnL,TPR_byLnE,precision_byLnE,TPR_bythree,precision_bythree), ncol=2, byrow=TRUE)
rownames(tprp) <- c("by DESeq2 & EdgeR","by DESeq2 & Limmavoom","by Limmavoom & EdgeR", "by all of these three tools")
colnames(tprp) <- c("TPR(Recall)", "Precision")
tprp <- as.table(tprp)
tprp
