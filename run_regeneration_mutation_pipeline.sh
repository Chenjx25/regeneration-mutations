#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash scripts/run_regeneration_mutation_pipeline.sh CONFIG STAGE

Stages:
  check             validate configuration, inputs, and software
  qc                FastQC and MultiQC
  map               BWA-MEM mapping, duplicate marking, optional GATK3 realignment
  bam-qc            samtools flagstat/stats for final BAMs
  call-hc           GATK3 HaplotypeCaller GVCFs and joint genotyping
  call-ug           GATK3 UnifiedGenotyper per sample and VCF merge
  candidate-sites   union SNV sites from both callers
  readcounts        three pileup/readcount sets: MQ20.AR, MQ0.AR, MQ20.NAR
  screen-snv        refill allele depths and screen candidates by sample groups
  combine-snv       combine readcount modes, callers, and group structures
  all               run all stages above in order
EOF
}

[[ $# -eq 2 ]] || { usage >&2; exit 2; }
CONFIG=$1
STAGE=$2
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }
need_file() { [[ -s "$1" ]] || die "missing or empty file: $1"; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "command not found: $1"; }
need_tool() {
    if [[ "$1" == */* ]]; then need_file "$1"; else need_cmd "$1"; fi
}

required_vars=(
    PROJECT_ID PROJECT_ROOT REFERENCE SAMPLES_TSV GROUP_A_NAME GROUP_A_FILE
    BAM_DIR HC_DIR UG_DIR READCOUNT_DIR RESULT_DIR QC_DIR THREADS
    GATK3_JAR VARSCAN_JAR FILL_VCF_DEPTH DETECT_MUTATIONS VCF_PROCESS MASK_VCF
    MIN_MAPQ MIN_BASEQ MIN_SUPPORT_DEPTH MAX_COMPARATOR_DEPTH
    MAX_COMPARATOR_TOTAL MAX_COMPARATOR_MISSING MIN_SUPPORT_PLUS
    MIN_SUPPORT_MINUS MIN_SITE_DEPTH MAX_SITE_DEPTH
)

check_config() {
    local name
    for name in "${required_vars[@]}"; do
        [[ -n "${!name:-}" ]] || die "configuration variable is unset: $name"
    done
    need_file "$REFERENCE"
    need_file "$SAMPLES_TSV"
    need_file "$GROUP_A_FILE"
    [[ -z "${GROUP_B_FILE:-}" ]] || need_file "$GROUP_B_FILE"
    need_file "$GATK3_JAR"
    need_file "$VARSCAN_JAR"
    for name in bwa samtools bcftools bgzip tabix bedtools vcf-annotate fastqc multiqc \
        java perl awk sort sed head tail tr cut uniq cp; do
        need_cmd "$name"
    done
    for name in "$FILL_VCF_DEPTH" "$DETECT_MUTATIONS" "$VCF_PROCESS" "$MASK_VCF"; do
        need_tool "$name"
    done

    local header
    header=$(head -n 1 "$SAMPLES_TSV" | tr -d '\r')
    [[ "$header" == $'sample_id\tread1\tread2' ]] || \
        die "samples TSV header must be: sample_id<TAB>read1<TAB>read2"
    [[ "$GROUP_A_NAME" =~ ^[A-Za-z0-9_.-]+$ ]] || die "unsafe GROUP_A_NAME"
    if [[ -n "${GROUP_B_FILE:-}" ]]; then
        [[ "${GROUP_B_NAME:-}" =~ ^[A-Za-z0-9_.-]+$ ]] || die "unsafe GROUP_B_NAME"
    fi

    local duplicate_sample
    duplicate_sample=$(tail -n +2 "$SAMPLES_TSV" | sed '/^[[:space:]]*$/d; /^[[:space:]]*#/d' | \
        cut -f1 | sort | uniq -d | head -n 1)
    [[ -z "$duplicate_sample" ]] || die "duplicate sample ID: $duplicate_sample"

    local sample read1 read2 extra
    while IFS=$'\t' read -r sample read1 read2 extra; do
        [[ -n "$sample" ]] || continue
        [[ "$sample" != *[[:space:]]* ]] || die "sample ID contains whitespace: $sample"
        need_file "$read1"
        need_file "$read2"
    done < <(tail -n +2 "$SAMPLES_TSV" | sed '/^[[:space:]]*$/d; /^[[:space:]]*#/d')
    log "configuration and dependencies passed"
}

prepare_dirs() {
    mkdir -p "$BAM_DIR" "$HC_DIR" "$UG_DIR" "$READCOUNT_DIR" "$RESULT_DIR" "$QC_DIR"
}

sample_rows() {
    tail -n +2 "$SAMPLES_TSV" | sed '/^[[:space:]]*$/d; /^[[:space:]]*#/d'
}

final_bam() {
    local sample=$1
    if [[ "${RUN_REALIGNMENT:-true}" == true ]]; then
        printf '%s/%s.dedup.realn.bam\n' "$BAM_DIR" "$sample"
    else
        printf '%s/%s.dedup.bam\n' "$BAM_DIR" "$sample"
    fi
}

run_qc() {
    prepare_dirs
    local fastqs=()
    local sample read1 read2 extra
    while IFS=$'\t' read -r sample read1 read2 extra; do
        fastqs+=("$read1" "$read2")
    done < <(sample_rows)
    fastqc --noextract --threads "$THREADS" --outdir "$QC_DIR" "${fastqs[@]}"
    multiqc --force --outdir "$QC_DIR" "$QC_DIR"
}

run_map() {
    prepare_dirs
    local sample read1 read2 extra sorted dedup target final log_file
    while IFS=$'\t' read -r sample read1 read2 extra; do
        sorted="$BAM_DIR/${sample}.sort.bam"
        dedup="$BAM_DIR/${sample}.dedup.bam"
        log_file="$BAM_DIR/${sample}.bwa.log"
        log "mapping $sample"
        bwa mem -t "$THREADS" -M \
            -R "@RG\\tID:${sample}\\tLB:${sample}\\tPL:ILLUMINA\\tPU:${sample}\\tSM:${sample}" \
            "$REFERENCE" "$read1" "$read2" 2>"$log_file" | \
            samtools fixmate -@ "$THREADS" -m -O bam - - | \
            samtools sort -@ "$THREADS" -o "$sorted" -
        samtools markdup -@ "$THREADS" -O BAM -f "$BAM_DIR/${sample}.markdup.stats" \
            "$sorted" "$dedup"
        samtools index "$dedup"

        if [[ "${RUN_REALIGNMENT:-true}" == true ]]; then
            target="$BAM_DIR/${sample}.realn.intervals"
            final="$BAM_DIR/${sample}.dedup.realn.bam"
            java -jar "$GATK3_JAR" -R "$REFERENCE" -T RealignerTargetCreator \
                -nt "$THREADS" -I "$dedup" -o "$target" >>"$log_file" 2>&1
            java -jar "$GATK3_JAR" -R "$REFERENCE" -T IndelRealigner \
                -targetIntervals "$target" -I "$dedup" -o "$final" >>"$log_file" 2>&1
            samtools index "$final"
        fi
    done < <(sample_rows)
}

run_bam_qc() {
    prepare_dirs
    local sample read1 read2 extra bam
    while IFS=$'\t' read -r sample read1 read2 extra; do
        bam=$(final_bam "$sample")
        need_file "$bam"
        samtools flagstat -@ "$THREADS" "$bam" >"$QC_DIR/${sample}.flagstat.txt"
        samtools stats -@ "$THREADS" "$bam" >"$QC_DIR/${sample}.samtools.stats.txt"
    done < <(sample_rows)
}

run_hc() {
    prepare_dirs
    local sample read1 read2 extra bam gvcf
    local variants=()
    while IFS=$'\t' read -r sample read1 read2 extra; do
        bam=$(final_bam "$sample")
        gvcf="$HC_DIR/${sample}.hc.g.vcf"
        java -jar "$GATK3_JAR" -R "$REFERENCE" -T HaplotypeCaller \
            --emitRefConfidence GVCF --variant_index_type LINEAR \
            --variant_index_parameter 128000 -dt NONE -I "$bam" -o "$gvcf" \
            >"$HC_DIR/${sample}.hc.log" 2>&1
        bgzip -f "$gvcf"
        tabix -f -p vcf "${gvcf}.gz"
        variants+=( -V "${gvcf}.gz" )
    done < <(sample_rows)
    java -jar "$GATK3_JAR" -R "$REFERENCE" -T GenotypeGVCFs -nt 1 \
        -stand_call_conf 30.0 "${variants[@]}" -o "$HC_DIR/${PROJECT_ID}.hc.vcf" \
        >"$HC_DIR/${PROJECT_ID}.hc.joint.log" 2>&1
    bgzip -f "$HC_DIR/${PROJECT_ID}.hc.vcf"
    tabix -f -p vcf "$HC_DIR/${PROJECT_ID}.hc.vcf.gz"
}

run_ug() {
    prepare_dirs
    local sample read1 read2 extra bam vcf
    local vcfs=()
    while IFS=$'\t' read -r sample read1 read2 extra; do
        bam=$(final_bam "$sample")
        vcf="$UG_DIR/${sample}.ug.vcf"
        java -jar "$GATK3_JAR" -R "$REFERENCE" -T UnifiedGenotyper \
            -glm BOTH -nt "$THREADS" -stand_call_conf 30.0 \
            -rf MappingQuality -mmq "$MIN_MAPQ" -dt NONE -I "$bam" -o "$vcf" \
            >"$UG_DIR/${sample}.ug.log" 2>&1
        bgzip -f "$vcf"
        tabix -f -p vcf "${vcf}.gz"
        vcfs+=("${vcf}.gz")
    done < <(sample_rows)
    bcftools merge --missing-to-ref -Oz -o "$UG_DIR/${PROJECT_ID}.ug.vcf.gz" "${vcfs[@]}"
    tabix -f -p vcf "$UG_DIR/${PROJECT_ID}.ug.vcf.gz"
}

run_candidate_sites() {
    prepare_dirs
    local hc="$HC_DIR/${PROJECT_ID}.hc.vcf.gz"
    local ug="$UG_DIR/${PROJECT_ID}.ug.vcf.gz"
    need_file "$hc"; need_file "$ug"
    {
        bcftools view -v snps "$hc" -Ou | bcftools query -f '%CHROM\t%POS0\t%END\n'
        bcftools view -v snps "$ug" -Ou | bcftools query -f '%CHROM\t%POS0\t%END\n'
    } | sort -k1,1 -k2,2n | bedtools merge -i - >"$RESULT_DIR/${PROJECT_ID}.snv_candidates.bed"
    {
        bcftools view -v indels "$hc" -Ou | bcftools query -f '%CHROM\t%POS0\t%END\n'
        bcftools view -v indels "$ug" -Ou | bcftools query -f '%CHROM\t%POS0\t%END\n'
    } | sort -k1,1 -k2,2n | bedtools merge -i - >"$RESULT_DIR/${PROJECT_ID}.indel_candidates.bed"
}

make_readcounts() {
    local sample=$1 bam=$2 mode=$3 mapq=$4 anomalous=$5
    local pileup="$READCOUNT_DIR/${sample}.${mode}.mpileup"
    local output="$READCOUNT_DIR/${sample}.${mode}.readcounts"
    local opts=(-d 100000 -q "$mapq" -f "$REFERENCE" \
        -l "$RESULT_DIR/${PROJECT_ID}.snv_candidates.bed")
    [[ "$anomalous" == true ]] && opts=(-A "${opts[@]}")
    samtools mpileup "${opts[@]}" "$bam" | awk '$4 > 0' >"$pileup"
    java -jar "$VARSCAN_JAR" readcounts "$pileup" --min-base-qual "$MIN_BASEQ" \
        --min-coverage 1 --output-file "$output"
}

run_readcounts() {
    prepare_dirs
    local sample read1 read2 extra bam
    while IFS=$'\t' read -r sample read1 read2 extra; do
        bam=$(final_bam "$sample")
        make_readcounts "$sample" "$bam" MQ20.AR "$MIN_MAPQ" true
        make_readcounts "$sample" "$bam" MQ0.AR 0 true
        make_readcounts "$sample" "$bam" MQ20.NAR "$MIN_MAPQ" false
    done < <(sample_rows)

    local mode list
    for mode in MQ20.AR MQ0.AR MQ20.NAR; do
        list="$READCOUNT_DIR/${PROJECT_ID}.${mode}.list"
        : >"$list"
        while IFS=$'\t' read -r sample read1 read2 extra; do
            printf '%s\t%s\t%s\n' "$sample" "$sample" \
                "$READCOUNT_DIR/${sample}.${mode}.readcounts" >>"$list"
        done < <(sample_rows)
    done
}

refill_one() {
    local caller=$1 mode=$2 joint output list
    if [[ "$caller" == hc ]]; then joint="$HC_DIR/${PROJECT_ID}.hc.vcf.gz";
    else joint="$UG_DIR/${PROJECT_ID}.ug.vcf.gz"; fi
    list="$READCOUNT_DIR/${PROJECT_ID}.${mode}.list"
    output="$RESULT_DIR/${PROJECT_ID}.${caller}.snv.${mode}.vcf.gz"
    "$FILL_VCF_DEPTH" --vcf "$joint" --list "$list" --minimum-vcf --update-AD | \
        bgzip -c >"$output"
    tabix -f -p vcf "$output"
}

screen_one() {
    local caller=$1 mode=$2 group_name=$3 group_file=$4
    local input="$RESULT_DIR/${PROJECT_ID}.${caller}.snv.${mode}.vcf.gz"
    local output="$RESULT_DIR/${PROJECT_ID}.${caller}.snv.${mode}.${group_name}.mut.vcf"
    "$DETECT_MUTATIONS" -v "$input" \
        --max-cmp-depth "$MAX_COMPARATOR_DEPTH" \
        --max-cmp-total "$MAX_COMPARATOR_TOTAL" \
        --min-supp-depth "$MIN_SUPPORT_DEPTH" \
        --max-cmp-miss "$MAX_COMPARATOR_MISSING" \
        --min-supp-plus "$MIN_SUPPORT_PLUS" \
        --min-supp-minus "$MIN_SUPPORT_MINUS" -g "$group_file" | \
        vcf-annotate -f "c=${MIN_SITE_DEPTH},${MAX_SITE_DEPTH}" --fill-type | \
        "$MASK_VCF" --input - --seq "$REFERENCE" --add-pass --output "$output"
}

run_screen_snv() {
    prepare_dirs
    local caller mode
    for caller in hc ug; do
        for mode in MQ20.AR MQ0.AR MQ20.NAR; do
            refill_one "$caller" "$mode"
            screen_one "$caller" "$mode" "$GROUP_A_NAME" "$GROUP_A_FILE"
            if [[ -n "${GROUP_B_FILE:-}" ]]; then
                screen_one "$caller" "$mode" "$GROUP_B_NAME" "$GROUP_B_FILE"
            fi
        done
    done
}

combine_modes() {
    local caller=$1 group=$2 prefix="$RESULT_DIR/${PROJECT_ID}.${caller}.snv"
    "$VCF_PROCESS" --vcf "${prefix}.MQ20.NAR.${group}.mut.vcf" \
        --secondary-vcf "${prefix}.MQ20.AR.${group}.mut.vcf" \
        --primary-tag NAR --secondary-tag AR --intersect-tag PF | \
        "$VCF_PROCESS" --vcf - \
        --secondary-vcf "${prefix}.MQ0.AR.${group}.mut.vcf" \
        --primary-tag MQ20 --secondary-tag MQ0 --intersect-tag MQPASS \
        >"$RESULT_DIR/${PROJECT_ID}.${caller}.snv.${group}.combined.vcf"
}

combine_callers() {
    local group=$1
    "$VCF_PROCESS" --vcf "$RESULT_DIR/${PROJECT_ID}.hc.snv.${group}.combined.vcf" \
        --secondary-vcf "$RESULT_DIR/${PROJECT_ID}.ug.snv.${group}.combined.vcf" \
        --combine-rows 0 1 --compare-rows 2 3 4 \
        --primary-tag HC_GVCF --secondary-tag UG_Single \
        --intersect-tag 'UG_Single+HC_GVCF' \
        >"$RESULT_DIR/${PROJECT_ID}.snv.${group}.callers_combined.vcf"
}

run_combine_snv() {
    local caller
    for caller in hc ug; do combine_modes "$caller" "$GROUP_A_NAME"; done
    combine_callers "$GROUP_A_NAME"

    if [[ -n "${GROUP_B_FILE:-}" ]]; then
        for caller in hc ug; do combine_modes "$caller" "$GROUP_B_NAME"; done
        combine_callers "$GROUP_B_NAME"
        "$VCF_PROCESS" \
            --vcf "$RESULT_DIR/${PROJECT_ID}.snv.${GROUP_B_NAME}.callers_combined.vcf" \
            --secondary-vcf "$RESULT_DIR/${PROJECT_ID}.snv.${GROUP_A_NAME}.callers_combined.vcf" \
            --combine-rows 0 1 --compare-rows 2 3 4 \
            --primary-tag "SET_${GROUP_B_NAME}" --secondary-tag "SET_${GROUP_A_NAME}" \
            --intersect-tag SET_BOTH \
            >"$RESULT_DIR/${PROJECT_ID}.snv.mutations.combined.vcf"
    else
        cp "$RESULT_DIR/${PROJECT_ID}.snv.${GROUP_A_NAME}.callers_combined.vcf" \
            "$RESULT_DIR/${PROJECT_ID}.snv.mutations.combined.vcf"
    fi
}

run_stage() {
    case "$1" in
        check) check_config ;;
        qc) check_config; run_qc ;;
        map) check_config; run_map ;;
        bam-qc) check_config; run_bam_qc ;;
        call-hc) check_config; run_hc ;;
        call-ug) check_config; run_ug ;;
        candidate-sites) check_config; run_candidate_sites ;;
        readcounts) check_config; run_readcounts ;;
        screen-snv) check_config; run_screen_snv ;;
        combine-snv) check_config; run_combine_snv ;;
        all)
            check_config
            run_qc; run_map; run_bam_qc; run_hc; run_ug
            run_candidate_sites; run_readcounts; run_screen_snv; run_combine_snv
            ;;
        *) usage >&2; die "unknown stage: $1" ;;
    esac
}

run_stage "$STAGE"
log "stage completed: $STAGE"
