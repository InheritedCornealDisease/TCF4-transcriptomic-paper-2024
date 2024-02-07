library(data.table) 
library(DEXSeq)
library(BiocParallel)
library(org.Hs.eg.db)
library(tidyr)
library(dplyr)


#################
# read files in #
#################

inDir <- "../DEXseq_fromNAS/counts/NE vs E"
outdir <- '../outputs/dexseq_covariates/EvNE'
countFiles = list.files(inDir, pattern= "clean.txt$", full.names=TRUE)
basename(countFiles)

gfDir = "../ReferenceGenomes"
flattenedFile = list.files(gfDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)

allSamples <- read.csv("../sample_info_w_cov.csv")

sampleTable <- subset(allSamples, allSamples$condition %in% c('E_FECD', 'NE_FECD'))

EvNEdxd=DEXSeqDataSetFromHTSeq(
  countFiles, sampleData = sampleTable, design= ~ sample + exon + sex:exon + ethnicity:exon + condition:exon, flattenedfile = flattenedFile)

EvNEdxd = estimateSizeFactors(EvNEdxd)


###################################
# specify model, run DEU analysis #
###################################

formulaFullModel    =  ~ sample + exon + sex:exon + ethnicity:exon + condition:exon
formulaReducedModel =  ~ sample + exon + sex:exon + ethnicity:exon 

EvNEdxd = estimateDispersions(EvNEdxd, formula = formulaFullModel, BPPARAM=BPPARAM)

EvNEdxd = testForDEU( EvNEdxd, 
    reducedModel = formulaReducedModel, 
        fullModel = formulaFullModel, BPPARAM=BPPARAM)

EvNEdxd = estimateExonFoldChanges(EvNEdxd, fitExpToVar = "condition", BPPARAM=BPPARAM)

EvNEdxr1 = DEXSeqResults(EvNEdxd)
EvNEdxr1
mcols(EvNEdxr1)$description

table(EvNEdxr1$padj<0.05)
table(tapply(EvNEdxr1$padj<0.05,EvNEdxr1$groupID,any))


#########################
# plots, results tables #
#########################

# plot TCF4
pdf(file.path(outdir, 'TCF4_EvNE_cov.pdf', fsep = '/'))
plotDEXSeq(EvNEdxr1, "ENSG00000196628.20", expression=FALSE, norCounts=TRUE,
   legend=TRUE)
dev.off()

#Visualising full results in HTML
DEXSeqHTML(EvNEdxr1, FDR=0.05, color=c("#FF000080", "#0000FF80"), path = path = file.path(outdir, 'DEXSeqReport_EvNE_cov', fsep = '/'))

#making result tables
#adding columns with gene names
EvNEdxr1.df <- as.data.frame(EvNEdxr1)
columns(org.Hs.eg.db)
ens.str <- substr(rownames(EvNEdxr1.df), 1, 15)
EvNEdxr1.df$symbol <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
EvNEdxr1.df$genename <- mapIds(org.Hs.eg.db,
                               keys=ens.str,
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first")
EvNEdxr1.df$entrez <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
names(EvNEdxr1.df)[1] <- "Gene"                   
head(EvNEdxr1.df)
setDT(EvNEdxr1.df)

#exporting all results
EvNEdxr1.df.rmNA <- replace_na(EvNEdxr1.df, list(padj=1))
head(EvNEdxr1.df)
EvNEdxr1.df.rmNA.pord <- EvNEdxr1.df.rmNA[order(EvNEdxr1.df.rmNA$padj), ]
head(EvNEdxr1.df.rmNA.pord)
fwrite(EvNEdxr1.df.rmNA.pord, file.path(outdir, "DEXSeq_EvNE_all.csv", fsep = '/'))

#exporting just sig events
EvNErdxr1.df.rmNA.pord <- EvNEdxr1.df.rmNA.pord %>% filter(padj<=0.05)
head(EvNErdxr1.df.rmNA.pord)
fwrite(EvNErdxr1.df.rmNA.pord, file.path(outdir, "DEXSeq_EvNE_sig.csv", fsep = '/'))
