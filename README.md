Splicing Index Analysis Pipeline

This pipeline calculates the Splicing Index (SI) to assess exon inclusion levels in protein-coding and 
lncRNA genes using RNA-seq BAM files and gene annotations. It measures how frequently specific exons are
 included in mature transcripts based on RNA-seq read coverage.
 
 Directory Structure
 
├── preprocess_features.sh          # Step 1: Prepare exon and intron features
├── calculate_coverage.sh           # Step 2: Calculate read coverage using Bedtools
├── splicing_index.R                # Step 3: Compute splicing index and generate plots
├── splicing_features/              # Contains exon/intron BEDs and split BAMs
├── splicing_coverage/              # Contains coverage results
├── splicing_index_results.csv      # Final SI values for all genes
├── splicing_index_comparison.pdf   # Boxplot: protein-coding vs lncRNA
├── splicing_index_density.pdf      # Density plot of SI values
└── statistical_results.txt         # Wilcoxon test results and summary statistics



Requirements Software

samtools
bedtools
awk
R and Rscript
R Packages
Install required R packages by running:
R
install.packages(c("ggplot2", "dplyr"))


How to Run the Pipeline
Run the three steps sequentially in a Linux terminal:

Step 1: Preprocess Features

sh preprocess_features.sh

This script:

Splits the input BAM into:
Gapped reads (spliced, containing introns)
Ungapped reads (unspliced)
Extracts exon annotations for:
Protein-coding genes
lncRNA genes

Defines 3' splice site (3'SS) regions:
Last intronic base + 25 bp downstream

Input files:

ENCFF000FEH.bam – RNA-seq BAM file
gencode.v27.annotation.gtf – Gene annotation (GTF)

Output (in splicing_features/):
gapped.bam, ungapped.bam – Split BAMs
exon_pc.bed, exon_lnc.bed – Exon annotations
intron_3ss_pc.bed, intron_3ss_lnc.bed – 3' splice site regions

Step 2: Calculate Coverage

sh calculate_coverage.sh

This script:

Uses bedtools coverage to compute:

Inclusion counts: gapped reads overlapping exons

Exclusion counts: ungapped reads overlapping 3'SS introns

Extracts and formats gene-level coverage data

Output (in splicing_coverage/):

Raw coverage:

gapped_pc.coverage, gapped_lnc.coverage

ungapped_pc.coverage, ungapped_lnc.coverage

Cleaned count tables:

included_pc.csv, included_lnc.csv

excluded_pc.csv, excluded_lnc.csv

Step 3: Splicing Index Analysis in R

Rscript splicing_index.R

This R script:

Merges coverage data

Filters genes with sufficient inclusion and exclusion support

Calculates the Splicing Index (SI) using:


SI = Inclusion+Exclusion /Inclusion

 
Visualizes SI:

Boxplot comparing protein-coding vs lncRNA
Density plot of SI distribution
Performs Wilcoxon rank sum test to assess differences

Output:

splicing_index_results.csv – Per-gene SI values
splicing_index_comparison.pdf – Boxplot
splicing_index_density.pdf – Density plot
statistical_results.txt – Test statistics and summaries

Reference : 

Krchnáková Z, Thakur PK, Krausová M, Bieberstein N, Haberman N, Müller-McNicoll M, Stanek D. Splicing of long non-coding RNAs primarily depends on polypyrimidine tract and 5' splice-site sequences due to weak interactions with SR proteins. Nucleic Acids Res. 2019 Jan 25;47(2):911-928. doi: 10.1093/nar/gky1147. PMID: 30445574; PMCID: PMC6344860.


