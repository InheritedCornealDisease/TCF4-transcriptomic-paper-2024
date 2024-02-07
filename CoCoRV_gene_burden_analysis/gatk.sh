# following steps in https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode

### 1. Variant Calling with HaplotypeCaller

# first half of the FECD cases, where the file name has 6 fields

CSV=../data/NE_FECD_samples_FINAL1.csv 
my_dir=../NiuzhengChai

while read p; do
    ID=$(echo ${p} | cut -f 3 -d "," | cut -f 6 -d "_")
    echo ${ID}
    LOC=$( echo "${p}" | cut -f 4 -d ",")
    echo ${LOC}
    BAM=$(find ${LOC} -name "*${ID}*final.bam" -o -name "*${ID}*sorted_unique.bam")
    BAM=${BAM//"/media"/}
    echo $BAM 
    
   docker run --rm -v ${my_dir}:/gatk/my_data broadinstitute/gatk:latest ./gatk --java-options "-Xmx8g" HaplotypeCaller  \
   -R my_data/${my_dir}/human_g1k_v37.fasta \
   -I my_data/${BAM} \
   -O my_data/${my_dir}/raw_variants_${ID}.g.vcf \
   -ERC GVCF \
   --native-pair-hmm-threads 8 \
   --standard-min-confidence-threshold-for-calling 30.0

done <${CSV}

# second half of the FECD cases, where the file name has 5 fields

CSV=../data/NE_FECD_samples_FINAL2.csv 

while read p; do
    ID=$(echo ${p} | cut -f 3 -d "," | cut -f 5 -d "_" )
    echo ${ID}
    LOC=$( echo "${p}" | cut -f 4 -d ",")
    echo ${LOC}
    BAM=$(find ${LOC} -name "*${ID}*final.bam" -o -name "*${ID}_*sorted_unique.bam" -o -name "*${ID}*recal_reads.bam")
    BAM=${BAM//"/media"/}
    echo $BAM 
    
   docker run --rm -v ${my_dir}:/gatk/my_data broadinstitute/gatk:latest ./gatk --java-options "-Xmx8g" HaplotypeCaller  \
   -R my_data/${my_dir}/human_g1k_v37.fasta \
   -I my_data/${BAM} \
   -O my_data/${my_dir}/raw_variants_${ID}.g.vcf \
   -ERC GVCF \
   --native-pair-hmm-threads 8 \
   --standard-min-confidence-threshold-for-calling 30.0

done <${CSV}


### 2. Data Aggregation with GenomicsDBImport

docker run --name nz_gatk --rm -v ${my_dir}:/gatk/my_data broadinstitute/gatk:latest ./gatk --java-options "-Xms128g -Xmx128g" GenomicsDBImport \
       --genomicsdb-workspace-path my_data/dbi/dbi_workspace \
       -L my_data/intervals.bed \
       -ip 100 \
       --sample-name-map my_data/cohort.sample_map \
       --tmp-dir my_data/dbi/dbi_tmp \
       --merge-input-intervals \
       --batch-size 50 \
       --max-num-intervals-to-import-in-parallel 16 \
       --genomicsdb-shared-posixfs-optimizations


### 3. Joint Genotyping with GenotypeVCFs

docker run --name nz_gatk --rm -v ${my_dir}:/gatk/my_data broadinstitute/gatk:latest ./gatk --java-options "-Xms128g -Xmx128g" GenotypeGVCFs \
   -R my_data/general_tools/reference_sequences/human_g1k_v37.fasta \
   -V gendb://my_data/cocorv/data/gatk/dbi/dbi_workspace \
   -O my_data/cocorv/data/gatk/raw_case_vcf.gz \
   --tmp-dir my_data/cocorv/data/gatk/genotype_tmp


### 4a) Variant Recalibration step 1: VariantRecalibrator

# snps
docker run --name nz_gatk --rm -v ${my_dir}:/gatk/my_data broadinstitute/gatk:latest ./gatk --java-options "-Xms128g -Xmx128g" VariantRecalibrator \
   -R my_data/general_tools/reference_sequences/human_g1k_v37.fasta \
   -V my_data/cocorv/data/gatk/raw_case.vcf.gz \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 my_data/cocorv/data/gatk/variantRecal/resources/hapmap_3.3.b37.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 my_data/cocorv/data/gatk/variantRecal/resources/1000G_omni2.5.b37.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 my_data/cocorv/data/gatk/variantRecal/resources/1000G_phase1.snps.high_confidence.b37.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 my_data/cocorv/data/gatk/variantRecal/resources/Homo_sapiens_assembly19.dbsnp138.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -tranche 89.9 -tranche 89.0 -tranche -88 \
   -mode SNP \
   -O my_data/cocorv/data/gatk/variantRecal/snps_output.recal \
   --tranches-file my_data/cocorv/data/gatk/variantRecal/snps_output.tranches \
   --rscript-file my_data/cocorv/data/gatk/variantRecal/snps_output.plots.R

# indels
docker run --name nz_gatk --rm -v ${my_dir}:/gatk/my_data broadinstitute/gatk:latest ./gatk --java-options "-Xms128g -Xmx128g" VariantRecalibrator \
   -R my_data/general_tools/reference_sequences/human_g1k_v37.fasta \
   -V my_data/cocorv/data/gatk/raw_case.vcf.gz \
   --resource:mills,known=false,training=true,truth=true,prior=12.0 my_data/cocorv/data/gatk/variantRecal/resources/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
   --resource:axiom,known=false,training=true,truth=false,prior=10.0 my_data/cocorv/data/gatk/variantRecal/resources/Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 my_data/cocorv/data/gatk/variantRecal/resources/Homo_sapiens_assembly19.dbsnp138.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -tranche 89.9 -tranche 89.0 -tranche -88 \
   -mode INDEL \
   -O my_data/cocorv/data/gatk/variantRecal/indel_output.recal \
   --tranches-file my_data/cocorv/data/gatk/variantRecal/indel_output.tranches \
   --rscript-file my_data/cocorv/data/gatk/variantRecal/indel_output.plots.R


### 4b) Variant Recalibration step 2: ApplyVQSR 

# snp
docker run --name nz_gatk --rm -v ${my_dir}:/gatk/my_data broadinstitute/gatk:latest ./gatk --java-options "-Xms128g -Xmx128g" ApplyVQSR \
   -V my_data/raw_case.vcf.gz \
   -O my_data/cases_vqsr_snps.vcf.gz \
   --truth-sensitivity-filter-level 90.0 \
   --tranches-file my_data/variantRecal/snps_output.tranches \
   --create-output-variant-index true \
   --recal-file my_data/variantRecal/snps_output.recal \
   -mode SNP


# indel
docker run --name nz_gatk --rm -v ${my_dir}:/gatk/my_data broadinstitute/gatk:latest ./gatk --java-options "-Xms128g -Xmx128g" ApplyVQSR \
   -V my_data/cases_vqsr_snps.vcf.gz \
   -O my_data/cases_vqsr_snps_indels.vcf.gz \
   --truth-sensitivity-filter-level 90.0 \
   --tranches-file my_data/variantRecal/indel_output.tranches \
   --create-output-variant-index true \
   --recal-file my_data/variantRecal/indel_output.recal \
   -mode INDEL


### 5) normalise case vcf file
bcftools norm -m -any -f ${my_dir}/Homo_sapiens_assembly19.fasta ${my_dir}/cases_norm.vcf | bgzip > ${my_dir}/cases_norm.vcf.gz
