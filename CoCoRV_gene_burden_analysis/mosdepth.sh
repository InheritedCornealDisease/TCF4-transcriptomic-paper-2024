# first half of the FECD cases, where the external id has 6 parts

CSV=../data/NE_FECD_samples_FINAL1.csv 
DIR=../data/exome_wide/mosdepth

while read p; do
    ID=$(echo ${p} | cut -f 3 -d "," | cut -f 6 -d "_")
    echo ${ID}
    LOC=$( echo "${p}" | cut -f 4 -d ",")
    echo ${LOC}
    BAM=$(find ${LOC} -name "*${ID}*final.bam" -o -name "*${ID}*sorted_unique.bam")
    echo $BAM
    mkdir -p ${DIR}/${ID}
    mosdepth -Q 20 -t 4 ${DIR}/${ID}/out_${ID} ${BAM}
done <${CSV}

gunzip ${DIR}/*/*.per-base.bed.gz

# second half of the FECD cases, where the external id has 5 parts

CSV=../data/exome_wide/NE_FECD_samples_FINAL2.csv 
DIR=../data/exome_wide/mosdepth

while read p; do
    ID=$(echo ${p} | cut -f 3 -d "," | cut -f 5 -d "_" )
    echo ${ID}
    LOC=$( echo "${p}" | cut -f 4 -d ",")
    echo ${LOC}
    BAM=$(find ${LOC} -name "*${ID}*final.bam" -o -name "*${ID}_*sorted_unique.bam" -o -name "*${ID}*recal_reads.bam")
    echo $BAM
    mkdir -p ${DIR}/${ID}
    if [ ! -f ${DIR}/${ID}/out_${ID}.per-base.bed ]; then
        echo "running mosdepth"
        mosdepth -Q 20 -t 4 ${DIR}/${ID}/out_${ID} ${BAM}
    fi
done <${CSV}

gunzip ${DIR}/*/*.per-base.bed.gz