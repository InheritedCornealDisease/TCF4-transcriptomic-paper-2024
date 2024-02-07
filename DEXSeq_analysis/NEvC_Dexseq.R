library(data.table) 
library(DEXSeq)
library(BiocParallel)
library(org.Hs.eg.db)
library(tidyr)
library(dplyr)


#################
# read files in #
#################

inDir <- "../DEXseq_fromNAS/counts/NE vs C"
outdir <- '../outputs/dexseq_covariates/NEvC'
countFiles <- list.files(inDir, pattern = "clean.txt$", full.names=TRUE)
basename(countFiles)

gfDir <- "../ReferenceGenomes"
flattenedFile <- list.files(gfDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)

allSamples <- read.csv("../sample_info_w_cov.csv")
sampleTable <- subset(allSamples, allSamples$condition %in% c('NE_FECD', 'control'))

NEvCdxd=DEXSeqDataSetFromHTSeq(
  countFiles, sampleData = sampleTable, design= ~ sample + exon + sex:exon + ethnicity:exon + condition:exon, flattenedfile = flattenedFile)

NEvCdxd = estimateSizeFactors(NEvCdxd)


###################################
# specify model, run DEU analysis #
###################################

formulaFullModel    =  ~ sample + exon + sex:exon + ethnicity:exon + condition:exon
formulaReducedModel =  ~ sample + exon + sex:exon + ethnicity:exon 

BPPARAM = MulticoreParam(4)

NEvCdxd = estimateDispersions(NEvCdxd, formula = formulaFullModel, BPPARAM=BPPARAM)

NEvCdxd = testForDEU(NEvCdxd, 
    reducedModel = formulaReducedModel, 
        fullModel = formulaFullModel, BPPARAM=BPPARAM)

NEvCdxd = estimateExonFoldChanges(NEvCdxd, fitExpToVar = "condition", BPPARAM=BPPARAM)

NEvCdxr1 = DEXSeqResults(NEvCdxd)
NEvCdxr1
mcols(NEvCdxr1)$description

table(NEvCdxr1$padj<0.05)
table(tapply(NEvCdxr1$padj<0.05,NEvCdxr1$groupID,any)) 


#########################
# plots, results tables #
#########################

# plot TCF4
pdf(file.path(outdir, 'TCF4_NEvC_cov_splicing.pdf', fsep = '/'))
plotDEXSeq(NEvCdxr1, "ENSG00000196628.20", expression=FALSE, norCounts=TRUE,
   legend=TRUE, splicing=TRUE)
dev.off() 

#Visualising full results in HTML
DEXSeqHTML(NEvCdxr1, FDR=0.05, color=c("#FF000080", "#0000FF80"), path = file.path(outdir, 'DEXSeqReport_NEvC_cov', fsep = '/'))

#making result tables
#adding columns with gene names
NEvCdxr1.df <- as.data.frame(NEvCdxr1)
columns(org.Hs.eg.db)
ens.str <- substr(rownames(NEvCdxr1.df), 1, 15)
NEvCdxr1.df$symbol <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
NEvCdxr1.df$genename <- mapIds(org.Hs.eg.db,
                               keys=ens.str,
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first")
NEvCdxr1.df$entrez <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
names(NEvCdxr1.df)[1] <- "Gene"                   
head(NEvCdxr1.df)
setDT(NEvCdxr1.df)

#exporting all results
NEvCdxr1.df.rmNA <- replace_na(NEvCdxr1.df, list(padj=1))
head(NEvCdxr1.df)
NEvCdxr1.df.rmNA.pord <- NEvCdxr1.df.rmNA[order(NEvCdxr1.df.rmNA$padj), ]
head(NEvCdxr1.df.rmNA.pord)
fwrite(NEvCdxr1.df.rmNA.pord, file.path(outdir, "DEXSeq_NEvC_all.csv", fsep = '/'))

#exporting just sig events
NEvCrdxr1.df.rmNA.pord <- NEvCdxr1.df.rmNA.pord %>% filter(padj<=0.05)
head(NEvCrdxr1.df.rmNA.pord)
fwrite(NEvCrdxr1.df.rmNA.pord, file.path(outdir, "DEXSeq_NEvC_sig.csv", fsep = '/'))