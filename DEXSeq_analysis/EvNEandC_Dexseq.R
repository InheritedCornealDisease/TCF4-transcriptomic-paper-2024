library(data.table) 
library(DEXSeq)
library(BiocParallel)
library(org.Hs.eg.db)
library(tidyr)
library(dplyr)


#################
# read files in #
#################

inDir1 <- "../DEXseq_fromNAS/counts"
inDir2 <- "../DEXseq_fromNAS/counts/NE vs C"
outdir <- '../outputs/dexseq_covariates/EvNEandC'
countFiles1 = list.files(inDir1, pattern= "clean.txt$", full.names=TRUE)
countFiles2 = list.files(inDir2, pattern= "clean.txt$", full.names=TRUE)
countFilesAll <- append(countFiles1, countFiles2)
basename(countFilesAll)
countFiles <- countFilesAll[-c(8:11)] # remove duplicates
countFiles

gfDir <- "../ReferenceGenomes"
flattenedFile <- list.files(gfDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)

allSamples <- read.csv("../sample_info_w_cov.csv")
allSamples$condition[allSamples$condition != 'E_FECD'] <- 'not_expanded' # anything that is not Exp+ FECD, rename to 'not expanded'

dxd=DEXSeqDataSetFromHTSeq(
  countFiles, sampleData = allSamples, design= ~ sample + exon + sex:exon + ethnicity:exon + condition:exon, flattenedfile = flattenedFile)

dxd = estimateSizeFactors(dxd)


###################################
# specify model, run DEU analysis #
###################################

formulaFullModel    =  ~ sample + exon + sex:exon + ethnicity:exon + condition:exon
formulaReducedModel =  ~ sample + exon + sex:exon + ethnicity:exon 

BPPARAM = MulticoreParam(4)

dxd = estimateDispersions(dxd, formula = formulaFullModel, BPPARAM=BPPARAM)

dxd = testForDEU(dxd, 
    reducedModel = formulaReducedModel, 
        fullModel = formulaFullModel, BPPARAM=BPPARAM)

dxd = estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM=BPPARAM)

dxr1 = DEXSeqResults(dxd)
dxr1
mcols(dxr1)$description

table(dxr1$padj<0.05)
table(tapply(dxr1$padj<0.05,dxr1$groupID,any))


#########################
# plots, results tables #
#########################

pdf(file.path(outdir, 'TCF4_EvNEandC_cov.pdf', fsep = '/'))
plotDEXSeq( dxr1, "ENSG00000196628.20", expression=FALSE, norCounts=TRUE,
   legend=TRUE)
dev.off() 

#Visualising full results in HTML
DEXSeqHTML(dxr1, FDR=0.05, color=c("#FF000080", "#0000FF80"), path = file.path(outdir, 'DEXSeqReport_EvNEandC_cov', fsep = '/'))

#making result tables 
#adding columns with gene names
dxr1.df <- as.data.frame(dxr1)
columns(org.Hs.eg.db)
ens.str <- substr(rownames(dxr1.df), 1, 15)
dxr1.df$symbol <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="SYMBOL",
                             keytype="ENSEMBL",
                             multiVals="first")
dxr1.df$genename <- mapIds(org.Hs.eg.db,
                               keys=ens.str,
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first")
dxr1.df$entrez <- mapIds(org.Hs.eg.db,
                             keys=ens.str,
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
names(dxr1.df)[1] <- "Gene"                   
head(dxr1.df)
setDT(dxr1.df)

#exporting all results
dxr1.df.rmNA <- replace_na(dxr1.df, list(padj=1))
head(dxr1.df)
dxr1.df.rmNA.pord <- dxr1.df.rmNA[order(dxr1.df.rmNA$padj), ]
head(dxr1.df.rmNA.pord)
fwrite(dxr1.df.rmNA.pord, file.path(outdir, "DEXSeq_EvNEandC_all.csv", fsep = '/'))

#exporting just sig events
rdxr1.df.rmNA.pord <- dxr1.df.rmNA.pord %>% filter(padj<=0.05)
head(rdxr1.df.rmNA.pord)
fwrite(rdxr1.df.rmNA.pord, file.path(outdir, "DEXSeq_EvNEandC_sig.csv", fsep = '/'))