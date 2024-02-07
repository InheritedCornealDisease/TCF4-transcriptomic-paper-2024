out=../outputs/outdir
tmp=../outputs/tmp
case_vcf=../data/cases_norm.vcf.gz
ctrl_vcf=../data/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz
bed=../data/combined.dp10.bed

# (1a) remove CTRL snps that didn't PASS
echo 'Starting gnomAD PASS snps'
bcftools view -f PASS ${ctrl_vcf} | bcftools sort -o ${out}/gnomad_pass_snps.vcf

bgzip -@ 8 ${out}/gnomad_pass_snps.vcf
tabix -p vcf ${out}/gnomad_pass_snps.vcf.gz
echo 'gnomAD PASS snps finished'

# (1b) remove CASE snps that didn't PASS
echo 'Starting case PASS snps'
bcftools view -f PASS ${case_vcf} | bcftools sort -o ${tmp}/cases_pass_snps1.vcf

bgzip -@ 8  ${tmp}/cases_pass_snps1.vcf
tabix -p vcf  ${tmp}/cases_pass_snps1.vcf.gz
echo 'case PASS snps finished'

# (2) remove any CASE snps that didn't PASS in the CTRL vcf
echo 'Creating list of gnomAD FAIL snps'
bcftools view -e "FILTER='PASS'"  ${ctrl_vcf} | bcftools sort -o ${tmp}/failed_gnomad_snps.vcf

bgzip -@ 8 ${tmp}/failed_gnomad_snps.vcf
tabix -p vcf ${tmp}/failed_gnomad_snps.vcf.gz
echo 'Done creating list of gnomAD FAIL snps'

echo 'removing failed gnomAD snps from case files'
bcftools isec -c all -C ${tmp}/cases_pass_snps1.vcf.gz /path/to/tmp/failed_gnomad_snps.vcf.gz -p ${tmp}
mv ${tmp}/0000.vcf ${tmp}/cases_pass_snps2.vcf
bcftools sort ${tmp}/cases_pass_snps2.vcf -o ${out}/cases_pass_snps.vcf

bgzip -@ 8 ${out}/cases_pass_snps.vcf
tabix -p vcf ${out}/cases_pass_snps.vcf.gz
echo 'cases PASS/FAIL filtering complete'

# (3a) get the 10X snps in CTRLS
echo 'read harmonisation for CTRLS'
bcftools view -R ${bed} ${ctrl_vcf} | bcftools sort -o ${out}/gnomad_pass_dp10.vcf

bgzip -@ 8 ${out}/gnomad_pass_dp10.vcf
tabix -p vcf ${out}/gnomad_pass_dp10.vcf.gz
echo 'read harmonisation for CTRLS done'

# (3b) get the 10X snps in CASES
echo 'read harmonisation for CASES'
bcftools view -R ${bed} ${out}/cases_pass_snps.vcf.gz | bcftools sort -o ${out}/case_pass_dp10.vcf

bgzip -@ 8 ${out}/case_pass_dp10.vcf
tabix -p vcf ${out}/case_pass_dp10.vcf.gz
echo 'read harmonisation for CASES done'

# (4) VEP annotation
separate script called vep.sh

## edit header before proceeding (gnomADg_AF type = Float) ### bcftools reheader -h header.hr input.vcf > output.vcf
