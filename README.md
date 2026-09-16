# Regeneration-associated mutation detection pipeline

This repository contains a configurable WGS workflow for detecting candidate de novo mutations in plant regenerants. The public script was reorganized from the same analysis framework used for *Arabidopsis thaliana*, rice, and *Liriodendron* datasets. Species-specific paths, sample names, reference genomes, and comparison structures are supplied through configuration files rather than being embedded in the code.

## Scope

The workflow performs:

1. FASTQ quality control with FastQC and MultiQC.
2. BWA-MEM alignment, duplicate marking, and optional GATK3 local realignment.
3. Variant discovery with two legacy GATK3 calling strategies:
   - HaplotypeCaller in GVCF mode followed by joint genotyping;
   - single-sample UnifiedGenotyper followed by multi-sample merging.
4. Construction of the union of candidate SNV and indel loci.
5. Allele-depth recounting under three read-selection conditions:
   - `MQ20.AR`: mapping quality at least 20, including anomalous read pairs;
   - `MQ0.AR`: no mapping-quality threshold, including anomalous read pairs;
   - `MQ20.NAR`: mapping quality at least 20, proper-pair behavior retained.
6. Candidate mutation screening against user-provided biological comparison groups.
7. Integration across read-count conditions, variant callers, and alternative group structures.

The current public driver produces the combined candidate SNV VCF. The original downstream evidence assessment and final table filtering rely on the external utilities listed below. Their exact tested versions must be identifiable through stable citations or repository links before claiming complete end-to-end reproducibility.

## Repository layout

```text
.
├── README.md
├── UPLOAD_CHECKLIST.md
├── config.example.sh
├── samples.example.tsv
├── .gitignore
├── .gitattributes
├── metadata/
│   ├── README.md
│   ├── sample_id_mapping.example.tsv
│   └── software_versions.example.tsv
└── scripts/
    ├── run_regeneration_mutation_pipeline.sh
    ├── rename_sample_ids.py
    └── vendor/                         # add the original custom utilities here
```

## Required inputs

### Reference files

- reference genome FASTA;
- BWA index files;
- FASTA index (`.fai`);
- GATK sequence dictionary (`.dict`);
- optional tandem-repeat BED files for downstream false-positive assessment.

Large reference files should normally be obtained from their original repositories rather than committed to GitHub. Record the assembly name, release, download URL, and checksum in the repository.

### Sample manifest

Copy `samples.example.tsv` to `samples.tsv`. It must be tab-delimited and contain:

```text
sample_id    read1    read2
```

Use the exact sample identifiers present in the BAM read groups and group-definition files. FASTQ paths may be absolute or relative to the directory from which the pipeline is run.

### Biological comparison groups

The original mutation detector uses group files to distinguish focal regenerants from comparator or parental samples. The exact group files used for each species are part of the analysis metadata and should be deposited even though they are small.

This repository does not invent a replacement format because the parser belongs to `detect_mutations.pl`. Record the tested version and its `--help` output, cite or link the stable code release, and provide example group files with non-sensitive sample IDs.

For the three plant systems, use separate configuration files, for example:

```text
config/arabidopsis.sh
config/rice.sh
config/liriodendron.sh
metadata/arabidopsis_samples.tsv
metadata/rice_samples.tsv
metadata/liriodendron_samples.tsv
metadata/arabidopsis_groups_*.txt
metadata/rice_groups_*.txt
metadata/liriodendron_groups_*.txt
```

## Configuration

Copy and edit the example:

```bash
cp config.example.sh config.sh
```

Make the driver and deposited helper tools executable on Linux:

```bash
chmod +x scripts/run_regeneration_mutation_pipeline.sh scripts/vendor/*
```

Important settings include the project identifier, output directories, reference FASTA, sample manifest, group files, GATK3 and VarScan JAR files, custom utility paths, and filtering thresholds.

The defaults reproduce the principal thresholds in the original scripts:

| Parameter | Default | Meaning |
|---|---:|---|
| `MIN_MAPQ` | 20 | minimum mapping quality for the primary read counts |
| `MIN_BASEQ` | 20 | minimum base quality |
| `MIN_SUPPORT_DEPTH` | 5 | minimum supporting depth in a focal sample |
| `MAX_COMPARATOR_DEPTH` | 2 | maximum alternative depth in an individual comparator |
| `MAX_COMPARATOR_TOTAL` | 5 | maximum total alternative support among comparators |
| `MAX_COMPARATOR_MISSING` | 5 | maximum missing comparator genotypes |
| `MIN_SUPPORT_PLUS` | 1 | minimum forward-strand support for SNVs |
| `MIN_SUPPORT_MINUS` | 1 | minimum reverse-strand support for SNVs |
| `MIN_SITE_DEPTH` | 3 | lower site-depth filter used after screening |
| `MAX_SITE_DEPTH` | 150 | upper site-depth filter used after screening |

Do not change these values silently. If a species-specific analysis used different thresholds, place them in that species configuration and report them in the Methods.

## Software

The legacy analysis framework used the following software families:

- Bash, Java, Perl, Python;
- BWA;
- SAMtools, BCFtools, bgzip, and tabix;
- BEDTools;
- FastQC and MultiQC;
- GATK 3.x, including HaplotypeCaller, GenotypeGVCFs, UnifiedGenotyper, RealignerTargetCreator, and IndelRealigner;
- VarScan 2.x `readcounts`;
- VCFtools `vcf-annotate`;
- optional Picard and tandem-repeat tools for downstream assessment.

For reproducibility, add an `environment.yml`, container recipe, or exact version table based on the computing environment actually used. GATK3 and UnifiedGenotyper are retained here to document the original analysis; they should not be silently replaced with current GATK commands because that would define a different workflow.

## External custom utilities and provenance

The core custom mutation-detection utilities used by this workflow were previously reported in **Wang et al. (2019), Ren et al. (2021), and Wang et al. (2023)**. These publications should be cited in the Methods and Code availability sections. Before public release, replace these author-year citations with complete references and persistent code links or DOIs.

If the exact, unmodified versions used here are already available from a stable public repository, they do not need to be duplicated in this repository. Record the repository URL, release or commit, checksum, and license. If the versions are not publicly obtainable, deposit them when redistribution is permitted or describe the access restriction explicitly.

The following non-standard commands are called by the original scripts but were not present in the local `reg` directory inspected during repository preparation:

### Core mutation screening

- `fillVcfDepth.pl`
- `detect_mutations.pl`
- `vcf_process.pl`
- `maskVCF.pl`
- `flt_vcf_AD` (needed for the original indel branch)

Without the first four utilities, the public driver can perform QC, mapping, calling, candidate-site generation, and read counting, but it cannot reproduce the original mutation-screening results.

### False-positive evidence assessment

- `check_bam_supports.pl`
- `annotate_mut_vcf.pl`
- `map_pos2intervals.pl`

### Convenience commands used only for summaries

- `body`
- `multijoin`

The last two are not required by the reorganized driver and may be replaced with documented standard `awk`, Python, or R code. Every deposited custom utility should include its license or permission statement, author/provenance, version or commit identifier, and a minimal test input/output pair.

On the original Linux server, locate and checksum the commands before copying them:

```bash
for tool in fillVcfDepth.pl detect_mutations.pl vcf_process.pl maskVCF.pl \
  flt_vcf_AD check_bam_supports.pl annotate_mut_vcf.pl map_pos2intervals.pl; do
  command -v "$tool" || true
done

sha256sum scripts/vendor/* > metadata/custom_script_sha256.tsv
```

## Running the workflow

First validate the configuration and all inputs:

```bash
bash scripts/run_regeneration_mutation_pipeline.sh config.sh check
```

Run individual stages:

```bash
bash scripts/run_regeneration_mutation_pipeline.sh config.sh qc
bash scripts/run_regeneration_mutation_pipeline.sh config.sh map
bash scripts/run_regeneration_mutation_pipeline.sh config.sh bam-qc
bash scripts/run_regeneration_mutation_pipeline.sh config.sh call-hc
bash scripts/run_regeneration_mutation_pipeline.sh config.sh call-ug
bash scripts/run_regeneration_mutation_pipeline.sh config.sh candidate-sites
bash scripts/run_regeneration_mutation_pipeline.sh config.sh readcounts
bash scripts/run_regeneration_mutation_pipeline.sh config.sh screen-snv
bash scripts/run_regeneration_mutation_pipeline.sh config.sh combine-snv
```

Or run all currently implemented stages:

```bash
bash scripts/run_regeneration_mutation_pipeline.sh config.sh all
```

The principal result is:

```text
mutation_results/<PROJECT_ID>.snv.mutations.combined.vcf
```

## Sample-ID normalization

The original Arabidopsis script contained many dataset-specific `sed -i` substitutions. These should not be embedded in a cross-species pipeline. Store mappings in a two-column TSV:

```text
old_id    new_id
old_name_1    public_name_1
```

Then run:

```bash
python scripts/rename_sample_ids.py \
  --input input.vcf \
  --output renamed.vcf \
  --mapping metadata/sample_id_mapping.tsv \
  --column ID
```

## Reproducibility and validation

Before release, document the following for each species:

- number of FASTQ pairs, BAMs, and samples in the joint VCF;
- reference assembly release;
- exact sample/group definition files;
- software versions and custom-script checksums;
- candidate counts after each major filtering stage;
- final SNV and indel counts per sample;
- a checksum or summary comparison against the files used in the manuscript.


## Data availability



## Citation and license


