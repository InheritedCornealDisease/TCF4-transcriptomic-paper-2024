# docker pull ensemblorg/ensembl-vep


gnomadg=../data/gnomad.genomes.r2.1.1.sites.vcf.bgz
snv=../data/cadd/whole_genome_SNVs.tsv.gz
indel=../data/cadd/InDels.tsv.gz
my_dir=../NiuzhengChai

# download cache
# docker run -t -i -v ${my_dir}/general_tools/ens-vep:/data ensemblorg/ensembl-vep INSTALL.pl -a p --PLUGINS CADD

# add annotations
docker run --rm -v ${my_dir}:/data ensemblorg/ensembl-vep \
    vep     --cache \
            --dir_cache ${my_dir}/general_tools/ens-vep \
            --offline \
            --input_file ${my_dir}/gnomad_pass_dp10.vcf.gz \  # either case or gnomad vcf path
            --output_file STDOUT \  # pipe to file on NAS
            --vcf \
            --warning_file STDERR \
            --custom ${my_dir}/${gnomadg},gnomADg,vcf,exact,0,AF,AC,AN \
            --no_stats \
            --fork 36 \
            --af_gnomade \
            --canonical \
            --sift b \
            --polyphen b \
            --plugin CADD,${my_dir}/${snv},${my_dir}/${indel}
