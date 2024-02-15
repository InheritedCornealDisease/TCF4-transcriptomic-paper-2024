#!/bin/bash

# control vs expansion+
rmats.py --s1 ../data/fastq/C.txt --s2 ../data/fastq/E.txt --gtf ../data/fastq/references/Homo_sapiens.GRCh38.105.gtf --bi ../data/fastq/star_index -t paired --variable-read-length --readLength 100 --nthread 128 --od ../output/fastq/rmats/run1_CvE --tmp ../output/fastq/rmats/run1_CvE/tmp

# control vs expansion-
rmats.py --s1 ../data/fastq/C.txt --s2 ../data/fastq/NE.txt --gtf ../data/fastq/references/Homo_sapiens.GRCh38.105.gtf --bi ../data/fastq/star_index -t paired --variable-read-length --readLength 100 --nthread 128 --od ../output/fastq_start/rmats/run1_CvNE --tmp ../output/fastq_start/rmats/run1_CvNE/tmp

# expansion+ vs expansion-
rmats.py --s1 ../data/fastq/E.txt --s2 ../data/fastq/NE.txt --gtf ../data/fastq/references/Homo_sapiens.GRCh38.105.gtf --bi ../data/fastq/star_index -t paired --variable-read-length --readLength 100 --nthread 128 --od ../output/fastq_start/rmats/run1_EvNE --tmp ../output/fastq_start/rmats/run1_EvNE/tmp
