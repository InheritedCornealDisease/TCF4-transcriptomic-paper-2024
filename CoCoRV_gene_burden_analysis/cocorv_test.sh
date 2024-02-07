controlCount=../data/ctrl_snps.vcf.gz.gds
caseCount=../data/case_qc-and-normalised.vcf.gz.gds
sampleList=../data/samples_no_sas.txt
intersectBed=../data/combined.dp10_noXY.bed.gz
variantExclude=../cocorv/example/1KG/gnomAD.exclude.allow.segdup.lcr.v3.txt.gz # this is cocorv's 'blacklisted' variants
AFMax=1
maxAFPopmax=1
variantMissing=0.1
ACANConfig=../cocorv/example/1KG/stratified_config_gnomad_2.txt # modified cocorv's exmple to put our own sample info in
pLDControl=0.05
ancestryFile=../data/ethnicity_no_sas.txt
highLDVariantFile=./example/1KG/full_vs_gnomAD.p0.05.OR1.ignoreEthnicityInLD.rds # this is cocorv's precalculated linkage disequilibrium data, based off GRCh37
variantGroupCustom=../cocorv/example/1KG/cocorv_custom_filter_fxn.R # my custom fxn 
outputPrefix=../outputs/stratified_no-SAS_no-XY/run_21_0.005af_15cadd
annotationList=../data/anno.txt 

Rscript utilities/CoCoRV_wrapper.R \
  --sampleList ${sampleList} \
  --outputPrefix ${outputPrefix} \
  --maxAFPopmax ${maxAFPopmax} \
  --AFMax ${AFMax} \
  --bed ${intersectBed} \
  --variantMissing ${variantMissing} \
  --variantGroup "case_snps" \
  --variantGroupCustom ${variantGroupCustom} \
  --removeStar \
  --ACANConfig ${ACANConfig} \
  --caseGroup ${ancestryFile} \
  --groupColumn "SYMBOL" \
  --minREVEL 0 \
  --variantExcludeFile ${variantExclude}  \
  --checkHighLDInControl \
  --pLDControl ${pLDControl} \
  --highLDVariantFile ${highLDVariantFile} \
  --annotationUsed ${annotationList} \
  --fullCaseGenotype \
  --gnomADVersion 'v2genome' \
  ${controlCount} \
  ${caseCount}

echo "QQ plot of the dominant model, and estimated FDR"
  outputFile=${outputPrefix}.association.tsv
  Rscript utilities/QQPlotAndFDR.R ${outputFile} \
    ${outputFile}.dominant.nRep1000_quantile \
    --setID gene \
    --outColumns gene \
    --n 1000 \
    --pattern "case.*Mutation.*_DOM$|control.*Mutation.*_DOM$" \
    --nullBoxplot \
    --FDR \
    --ncore 8 \

