# Splicing analysis from short-read Bulk RNAseq with rMATS [^1] and maser [^2]
- Run rMATS analysis with the `fastq_run1.sh`. Input are FASTQ files.
- Get summary statistics from the rMATS results using `maser_summary.R` (more specifically, get number of significant splicing events in total and also by splice type)

## Versions
- Ensembl GRCh 38 release 105 for genome & transcriptome 
- STAR 2.7.10a
- rMATS-turbo 4.1.2
- maser 1.12.1
- R 4.1.2

[^1]: Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. Proc Natl Acad Sci U S A. 2014 Dec 23;111(51):E5593-601. doi: 10.1073/pnas.1419161111. Epub 2014 Dec 5. PMID: 25480548; PMCID: PMC4280593.

[^2]: F.T. Veiga D (2023). maser: Mapping Alternative Splicing Events to pRoteins. doi:10.18129/B9.bioc.maser, R package version 1.20.0, https://bioconductor.org/packages/maser. 
