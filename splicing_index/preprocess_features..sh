#!/bin/bash

BAM="ENCFF000FEH.bam"
GENCODE_GTF="gencode.v27.annotation.gtf"
OUTPUT_DIR="splicing_features"

mkdir -p ${OUTPUT_DIR}

# Split reads into gapped (spliced) and ungapped (unspliced)
samtools view -h ${BAM} | \
awk '{if($0 ~ /^@/ || $6 ~ /N/) print $0}' | \
samtools view -Sb - > ${OUTPUT_DIR}/gapped.bam

samtools view -h ${BAM} | \
awk '{if($0 ~ /^@/ || $6 !~ /N/) print $0}' | \
samtools view -Sb - > ${OUTPUT_DIR}/ungapped.bam

samtools index ${OUTPUT_DIR}/gapped.bam
samtools index ${OUTPUT_DIR}/ungapped.bam

# Extract exons of protein-coding genes
awk '$3=="exon" && $0~/protein_coding/ {print $0}' ${GENCODE_GTF} | \
grep -w "exon" | \
awk -F'\t' 'BEGIN{OFS="\t"} {
    split($9,a,";");
    for(i in a){if(a[i]~/gene_name/){split(a[i],b,"\""); gene=b[2]}}
    print $1,$4-1,$5,gene,".",$7
}' | sort -k1,1 -k2,2n > ${OUTPUT_DIR}/exon_pc.bed

# Extract exons of lncRNAs
awk '$3=="exon" && $0~/lncRNA/ {print $0}' ${GENCODE_GTF} | \
grep -w "exon" | \
awk -F'\t' 'BEGIN{OFS="\t"} {
    split($9,a,";");
    for(i in a){if(a[i]~/gene_name/){split(a[i],b,"\""); gene=b[2]}}
    print $1,$4-1,$5,gene,".",$7
}' | sort -k1,1 -k2,2n > ${OUTPUT_DIR}/exon_lnc.bed

# Define 3' splice site regions from introns: last intronic base + 25 bp downstream exonic region
bedtools subtract -a <(awk '$3=="gene"' ${GENCODE_GTF} | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$4-1,$5,"gene",".",$7}') \
    -b ${OUTPUT_DIR}/exon_pc.bed -s | \
awk 'BEGIN{OFS="\t"} {
    if($6=="+"){
        print $1,$3-1,$3+25,$4,".",$6
    } else {
        print $1,$2-25,$2+1,$4,".",$6
    }
}' > ${OUTPUT_DIR}/intron_3ss_pc.bed

bedtools subtract -a <(awk '$3=="gene"' ${GENCODE_GTF} | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$4-1,$5,"gene",".",$7}') \
    -b ${OUTPUT_DIR}/exon_lnc.bed -s | \
awk 'BEGIN{OFS="\t"} {
    if($6=="+"){
        print $1,$3-1,$3+25,$4,".",$6
    } else {
        print $1,$2-25,$2+1,$4,".",$6
    }
}' > ${OUTPUT_DIR}/intron_3ss_lnc.bed

echo "Preprocessing complete. Output in ${OUTPUT_DIR}/"
