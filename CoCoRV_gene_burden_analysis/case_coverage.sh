DIR=../data/exome_wide/mosdepth
TMP=../data/exome_wide/dp10/tmp
MERGED=../data/exome_wide/dp10

for bed in ${DIR}/*/*.per-base.bed
    do 
        ID=$(basename ${bed} | cut -d "." -f 1 | cut -d "_" -f 2)
        echo ${ID}
        awk '$4>10' ${bed} > ${TMP}/${ID}_dp10.bed
    done

# list bed files with space separation
BEDS=$( ls -1 ${TMP}/*.bed | tr '\n' ' ')

echo ${BEDS}

multiIntersectBed -i ${BEDS} > ${MERGED}/merged_cases_dp10.bed

# sorting based on chr and position
sort -k 1,1 -k2,2n ${MERGED}/merged_cases_dp10.bed > ${MERGED}/sorted_merged_cases_dp10.bed

# 40,440,246 rows

awk '$4>126 {print $1"\t"$2"\t"$3}' ${MERGED}/sorted_merged_cases_dp10.bed > ${MERGED}/filtered_sorted_merged_cases_dp10.bed #126 is 90% of 141

bedtools merge -i ${MERGED}/filtered_sorted_merged_cases_dp10.bed > ${MERGED}/final_cases_dp10.bed

bedtools intersect -a ../data/exome_wide/gnomad_data/gnomad.dp10.bed -b ${MERGED}/final_cases_dp10.bed | sort -k1,1n -k2,2n | bedtools merge -i stdin > ${MERGED}/combined.dp10.bed
