library("DESeq2")
library("stringr")
library("ggplot2")
library("ggpubr")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("tidyverse")
library("tibble")
library("AnnotationDbi")
library("tximeta")
library("org.Hs.eg.db")
library("ggbeeswarm")
library("biomaRt")
library("genefilter")
library("IHW")
library("PoiClaClu")
library("tidyr")


################
# Read data in #
################

#locate salmon folder Note: tximeta requires that the entire output directory of Salmon / alevin is present and unmodified in order to identify the provenance of the reference transcripts.
dir <- "../RNAseq.Batch1"
list.files(dir)
list.files(file.path(dir, "quants"))

coldata <- read.csv("../sample_info_w_covariates.csv", row.names=1, stringsAsFactors=FALSE)
coldata

coldata$names <- coldata$run
file.exists(coldata$files)

se <- tximeta(coldata, type = "salmon", skipSeqinfo = TRUE)
dim(se)
head(rownames(se))

#for gene level analysis we summarise the transcripts
gse <- summarizeToGene(se)
dim(gse)
head(rownames(gse))

#removing version names from ensembl IDs
rownames(gse) <- gsub("\\.[0-9]*$", "", rownames(gse))

#tximeta has created an object gse with three matrices: “counts” - the estimated fragment counts for each gene and sample, “abundance” - the estimated transcript abundances in TPM, and “length” - the effective gene lengths which include changes in length due to biases as well as due to transcript usage. The names of the assays can be examined with assayNames, and the assays themselves are stored as assays (a list of matrices). The first matrix in the list can be pulled out via assay.
gse
assayNames(gse)

#look at length assay
head(assay(gse), 3)

#total counts assay
#The counts are the first matrix, so we can examine them with just assay:
colSums(assay(gse))

#GRanges of the genes (from the left-most position of all the transcripts to the right-most position of all the transcripts).
rowRanges(gse)

#colData for the SummarizedExperiment reflects the data.frame that was provided to the tximeta function for importing the quantification data
colData(gse)

#define levels
gse$condition <- factor(gse$condition, levels = c("Control","FECD"))
gse$batch <- factor(gse$batch, levels = c(1))
gse$diagnosis <- factor(gse$diagnosis, levels = c("Control","FECD_E", "FECD_NE")) # FECD_E = expanded, FECD_NE = non-expanded
gse$sex <- factor(gse$sex, levels = c ('M', 'F'))
gse$ethnicity <- factor(gse$ethnicity, levels = c('C', 'BA')) # C = Caucasian, BA = Black African
levels(gse$condition)
levels(gse$batch)
levels(gse$diagnosis)

#check the millions of fragments that could be mapped by Salmon to the genes
round(colSums(assay(gse)) / 1e6, 1)

#pick design for later contrast (cluster or group contrasts)
dds <- DESeqDataSet(gse, design = ~ diagnosis + sex + ethnicity)

#normalise counts
dds <- estimateSizeFactors(dds)

#Which transformation to choose? The rlog tends to work well on small datasets (n < 30)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$diagnosis, rld$name, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#Poisson Distance (Witten 2011), implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples. 
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$diagnosis, dds$name, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

#PCA
pcaData <- plotPCA(rld, intgroup = c("diagnosis", "name", 'ethnicity', 'sex'), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggplot2)
svg('../outputs/pca_cov.svg')
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = ethnicity, shape = diagnosis)) +
  geom_point(size =3) +
  scale_shape_manual(values = c(4,8,15,16,17,18,21,22,3,42))  +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("rlog PCA") +
  labs(col="Sample", shape="Diagnosis") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_line(colour = "grey", linetype = "dashed"),
        axis.ticks.y = element_blank()) + 
  theme(legend.position="bottom", legend.box = "vertical", legend.box.just = "left")
p
dev.off()
        
#get loadings
ntop <- 500
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
rld_matrix <- t(assay(rld)[select, ])

rld_pca <- prcomp(rld_matrix) 
load <- as.data.frame(rld_pca$rotation)

intgroup = c("diagnosis", "name", 'ethnicity', 'sex')

intgroup.df <- as.data.frame(colData(rld)[, intgroup,
        drop = FALSE])
group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
        } else {
        colData(rld)[[intgroup]]
        }

d <- data.frame(rld_pca$x, group = group, intgroup.df, name = colnames(rld))
write.csv(d, '../outputs/PCs_cov.csv')

pc1 <- load[order(load$PC1, decreasing = TRUE),]
pc1 <- subset(pc1, select = PC1)
write.csv(pc1, '../outputs/pc1_loading_genes_cov.csv', row.names = TRUE)

pc2 <- load[order(load$PC2, decreasing = TRUE),]
pc2 <- subset(pc2, select = PC2)
write.csv(pc2, '../outputs/pc2_loading_genes_cov.csv', row.names = TRUE)


#########################
# Control vs Expansion+ #
#########################

dds$diagnosis <- relevel(dds$diagnosis, ref = "Control")
dds <- DESeq(dds)
resultsNames(dds)

res.IHW.cfe <- results(dds, name="diagnosis_FECD_E_vs_Control", filterFun=ihw, alpha = 0.05)
summary(res.IHW.cfe)
res.LFC.cfe <- lfcShrink(dds, coef="diagnosis_FECD_E_vs_Control", type="apeglm")
res.IHW.cfe$shrunkLFC <- res.LFC.cfe$log2FoldChange

#convert IHW results to dataframe
res.IHW.cfedf <- as.data.frame(res.IHW.cfe)
columns(org.Hs.eg.db)
ens.str <- substr(rownames(res.IHW.cfedf), 1, 15)
res.IHW.cfedf$symbol <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
res.IHW.cfedf$genename <- mapIds(org.Hs.eg.db,
                               keys=ens.str,
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first")
res.IHW.cfedf$entrez <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
res.IHW.cfedf <- merge(as.data.frame(res.IHW.cfedf), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res.IHW.cfedf)[1] <- "Gene"
head(res.IHW.cfedf)
write.csv(res.IHW.cfedf, file="../outputs/res.IHW.cfedf_cov.csv")

#remove NA gene names
res.IHW.cfedf.rmNA <- replace_na(res.IHW.cfedf, list(padj=1))
res.IHW.cfedf.rmNA <- res.IHW.cfedf.rmNA %>% drop_na(symbol)
## Write results
write.csv(res.IHW.cfedf.rmNA, file="../outputs/res.IHW.cfedf_cov.rmNA.csv")


#########################
# Control vs Expansion- #
#########################

res.IHW.cfne <- results(dds, name="diagnosis_FECD_NE_vs_Control", filterFun=ihw, alpha = 0.05) # ref level is still control
summary(res.IHW.cfne)
res.LFC.cfne <- lfcShrink(dds, coef="diagnosis_FECD_NE_vs_Control", type="apeglm")
res.IHW.cfne$shrunkLFC <- res.LFC.cfne$log2FoldChange

#convert IHW results to dataframe
res.IHW.cfnedf <- as.data.frame(res.IHW.cfne)
columns(org.Hs.eg.db)
ens.str <- substr(rownames(res.IHW.cfnedf), 1, 15)
res.IHW.cfnedf$symbol <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
res.IHW.cfnedf$genename <- mapIds(org.Hs.eg.db,
                               keys=ens.str,
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first")
res.IHW.cfnedf$entrez <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
res.IHW.cfnedf <- merge(as.data.frame(res.IHW.cfnedf), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res.IHW.cfnedf)[1] <- "Gene"
head(res.IHW.cfnedf)
write.csv(res.IHW.cfnedf, file="../outputs/res.IHW.cfnedf_cov.csv")

res.IHW.cfnedf.rmNA <- replace_na(res.IHW.cfnedf, list(padj=1))
res.IHW.cfnedf.rmNA <- res.IHW.cfnedf.rmNA %>% drop_na(symbol)
## Write results
write.csv(res.IHW.cfnedf.rmNA, file="../outputs/res.IHW.cfnedf_cov.rmNA.csv")


############################
# Expansion+ vs Expansion- #
############################

dds$diagnosis <- relevel(dds$diagnosis, ref = "FECD_E") # now Exp+ is ref level
dds <- DESeq(dds)
resultsNames(dds)

res.IHW.efne <- results(dds, name="diagnosis_FECD_NE_vs_FECD_E", filterFun=ihw, alpha = 0.05)
summary(res.IHW.efne)
res.LFC.efne <- lfcShrink(dds, coef="diagnosis_FECD_NE_vs_FECD_E", type="apeglm")
res.IHW.efne$shrunkLFC <- res.LFC.efne$log2FoldChange

#convert IHW results to dataframe
res.IHW.efnedf <- as.data.frame(res.IHW.efne)
columns(org.Hs.eg.db)
ens.str <- substr(rownames(res.IHW.efnedf), 1, 15)
res.IHW.efnedf$symbol <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
res.IHW.efnedf$genename <- mapIds(org.Hs.eg.db,
                               keys=ens.str,
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first")
res.IHW.efnedf$entrez <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
res.IHW.efnedf <- merge(as.data.frame(res.IHW.efnedf), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(res.IHW.efnedf)[1] <- "Gene"
head(res.IHW.efnedf)
write.csv(res.IHW.efnedf, file="../outputs/res.IHW.efnedf_cov.csv")

res.IHW.efnedf.rmNA <- replace_na(res.IHW.efnedf, list(padj=1))
res.IHW.efnedf.rmNA <- res.IHW.efnedf.rmNA %>% drop_na(symbol)
## Write results
write.csv(res.IHW.efnedf.rmNA, file="../outputs/res.IHW.efnedf_cov.rmNA.csv")


capture.output(sessionInfo(), file = '../outputs/deseq2_session_info.txt')