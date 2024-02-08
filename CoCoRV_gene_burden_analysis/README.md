# Gene burden analysis using CoCoRV [^1]

## Analysis steps
1. Calculate read depths of **case** data, starting from BAM files
```
./mosdepth.sh

# the files we're interested in are the *.per-base.bed.gz files
```
2. Make a bed file containing only positions where >90% of **case** samples have read depth >10, and then intersect with gnomAD's [coverage table](https://gnomad.broadinstitute.org/downloads#v2-coverage) (regions >90% of **control** samples with read depth >10)
```
./case_coverage.sh

# to get gnomAD regions
zcat ${gnomad} | tail -n+2 | awk '$7>0.9 {print $1"\t"($2-1)"\t"$2}' | bedtools merge -i stdin > ${out}/gnomad.dp10.bed
```
3. Do joint variant calling using Genome Analysis Toolkit (GATK)'s [recommended workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode) for cohorts of samples. Requires `docker`. Includes a `bcftools norm` normalisation step at the end.
```
./gatk.sh
```
4. After we have our **case** vcf, we can do quality control for both **case** and **control** files. 1) Remove snps that didn't PASS in **case** (likewise for **control**), 2) Remove any snps in **case** that didn't PASS in **control**, 3) Subset both **case** & **control** vcfs for regions with good coverage (from step 2).
```
./vcf_qc.sh
```
5. Annotate both **case** & **control** files with Variant Effect Predictor (VEP). For consistency, we re-annotate the **control** file as well. Requires `docker`.
```
./vep.sh    # depending on how your system is set up, you might have to write the outputs to STDOUT and then pipe to the directory where the file is going to be stored

# I have found that it is necessary to extract the resulting CSQ annotations into their own lines after VEP annotation
bcftools +split-vep -c SYMBOL,CADD_PHRED:Int,gnomADg_AF:Float,gnomADe_AF:Float -d ${dir}/${vcfAnnotated} -Oz -o ${dir}/case_split_old-header.vcf.gz

# split-vep isn't able to change the gnomADg_AF from character to float for some reason, so I had to use bcftools reheader as well
bcftools head ${dir}/case_split_old-header.vcf.gz > ${dir}/case-header.hr  # manually change the type to Float
bcftools reheader -h ${dir}/case-header.hr ${dir}/case_split_old-header.vcf.gz > ${dir}/case_snps.vcf.gz
```
6. Convert both vcf to GDS format (see CoCoRV's bitbucket [documentation](https://bitbucket.org/Wenan/cocorv/src/master/))
7. Since annotation was done with VEP instead of ANNOVAR, use custom function to define variants of interest
```
cocorv_custom_filter_fxn.R
```
8. Prepare other files according to CoCoRV's [instructions](https://bitbucket.org/Wenan/cocorv/src/master/), then run the wrapper script
```
./cocorv_test.sh
```



### Versions
- CoCoRV 1.1
- R 4.2.1 (Funny-Looking Kid)
- GATK 4.4.0.0
- mosdepth 0.3.3
- VEP release 110 (GRCh 37)
- CADD 1.6 (GRCh 37)
- bedtools 2.31.0
- gnomAD 2.1.1 (GRCh 37)
- bcftools/htslib 1.18
- Fasta for joint genotyping: human_g1k_v37.fasta (Md5 sum: ce84c872fc0072a885926823dcd0338)

[^1]: Chen, W., Wang, S., Tithi, S.S. et al. A rare variant analysis framework using public genotype summary counts to prioritize disease-predisposition genes. Nat Commun 13, 2592 (2022). https://doi.org/10.1038/s41467-022-30248-0
