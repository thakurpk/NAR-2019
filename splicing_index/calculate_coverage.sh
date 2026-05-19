#!/bin/bash
#SBATCH --job-name=splicing_coverage

FEATURE_DIR="splicing_features"
OUTPUT_DIR="splicing_coverage"

mkdir -p ${OUTPUT_DIR}

# Calculate coverage for spliced reads (exons)
bedtools coverage -s -counts -split \
    -a ${FEATURE_DIR}/exon_pc.bed \
    -b ${FEATURE_DIR}/gapped.bam \
    > ${OUTPUT_DIR}/gapped_pc.coverage

bedtools coverage -s -counts -split \
    -a ${FEATURE_DIR}/exon_lnc.bed \
    -b ${FEATURE_DIR}/gapped.bam \
    > ${OUTPUT_DIR}/gapped_lnc.coverage

# Calculate coverage for unspliced reads (3'SS regions)
bedtools coverage -s -counts -split \
    -a ${FEATURE_DIR}/intron_3ss_pc.bed \
    -b ${FEATURE_DIR}/ungapped.bam \
    > ${OUTPUT_DIR}/ungapped_pc.coverage

bedtools coverage -s -counts -split \
    -a ${FEATURE_DIR}/intron_3ss_lnc.bed \
    -b ${FEATURE_DIR}/ungapped.bam \
    > ${OUTPUT_DIR}/ungapped_lnc.coverage

# Prepare files for R analysis
awk '{print $4,$7}' ${OUTPUT_DIR}/gapped_pc.coverage > ${OUTPUT_DIR}/included_pc.csv
awk '{print $4,$7}' ${OUTPUT_DIR}/gapped_lnc.coverage > ${OUTPUT_DIR}/included_lnc.csv
awk '{print $4,$7}' ${OUTPUT_DIR}/ungapped_pc.coverage > ${OUTPUT_DIR}/excluded_pc.csv
awk '{print $4,$7}' ${OUTPUT_DIR}/ungapped_lnc.coverage > ${OUTPUT_DIR}/excluded_lnc.csv

echo "Coverage calculation complete. Output in ${OUTPUT_DIR}/"
