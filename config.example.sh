#!/usr/bin/env bash
# Copy this file to config.sh and edit the paths for one species/dataset.

PROJECT_ID="example_species"
PROJECT_ROOT="/path/to/project"
REFERENCE="/path/to/reference/genome.fa"
SAMPLES_TSV="/path/to/repository/samples.tsv"

# Two alternative sample-comparison structures used by detect_mutations.pl.
# GROUP_B_FILE may be left empty when only one comparison structure is needed.
GROUP_A_NAME="comparison_set_a"
GROUP_A_FILE="/path/to/groups_a.txt"
GROUP_B_NAME="comparison_set_b"
GROUP_B_FILE="/path/to/groups_b.txt"

FASTQ_DIR="${PROJECT_ROOT}/fastq"
BAM_DIR="${PROJECT_ROOT}/bam"
HC_DIR="${PROJECT_ROOT}/vcf_hc"
UG_DIR="${PROJECT_ROOT}/vcf_ug"
READCOUNT_DIR="${PROJECT_ROOT}/readcounts"
RESULT_DIR="${PROJECT_ROOT}/mutation_results"
QC_DIR="${PROJECT_ROOT}/qc"

THREADS=4
RUN_REALIGNMENT=true

# Legacy callers used in the original analysis.
GATK3_JAR="/path/to/GenomeAnalysisTK-3.7.jar"
VARSCAN_JAR="/path/to/VarScan.v2.3.6.jar"

# Custom utilities used by the original mutation-screening workflow.
# Put them in scripts/vendor/ or provide absolute paths here.
FILL_VCF_DEPTH="/path/to/fillVcfDepth.pl"
DETECT_MUTATIONS="/path/to/detect_mutations.pl"
VCF_PROCESS="/path/to/vcf_process.pl"
MASK_VCF="/path/to/maskVCF.pl"

# Filtering parameters retained from the original scripts.
MIN_MAPQ=20
MIN_BASEQ=20
MIN_SUPPORT_DEPTH=5
MAX_COMPARATOR_DEPTH=2
MAX_COMPARATOR_TOTAL=5
MAX_COMPARATOR_MISSING=5
MIN_SUPPORT_PLUS=1
MIN_SUPPORT_MINUS=1
MIN_SITE_DEPTH=3
MAX_SITE_DEPTH=150
