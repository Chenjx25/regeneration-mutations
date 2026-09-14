# GitHub upload checklist

## Upload now

- `README.md`
- `.gitignore`
- `.gitattributes`
- `scripts/run_regeneration_mutation_pipeline.sh`
- `scripts/rename_sample_ids.py`
- `config.example.sh`
- `samples.example.tsv`
- `metadata/sample_id_mapping.example.tsv`
- `metadata/software_versions.example.tsv`
- species-specific sample manifests with public IDs
- species-specific group-definition files
- sample-ID mapping tables
- reference/annotation version and checksum table
- filtering-stage count summary used to verify manuscript results

## Previously reported external utilities

The core utilities below were reported in Wang et al. (2019), Ren et al. (2021), and Wang et al. (2023). If the exact versions are already publicly accessible, cite and link the stable release instead of uploading duplicate copies. Record the version, commit/checksum, license, and whether the scripts were modified.

Core scripts required for the mutation calls:

- `fillVcfDepth.pl`
- `detect_mutations.pl`
- `vcf_process.pl`
- `maskVCF.pl`
- `flt_vcf_AD`

Scripts required for the reported false-positive assessment:

- `check_bam_supports.pl`
- `annotate_mut_vcf.pl`
- `map_pos2intervals.pl`

Optional convenience tools:

- `body`
- `multijoin`

For every custom script, record provenance, redistribution permission, version/checksum, command-line help, dependencies, and a minimal test case.

## Add before public release

- exact software versions or an environment/container definition
- `LICENSE`
- `CITATION.cff`
- manuscript citation or preprint DOI
- data repository accessions
- a small test dataset or toy VCF with expected output
- a regression-test summary for Arabidopsis, rice, and Liriodendron

## Do not upload to ordinary Git history

- raw FASTQ files
- BAM/CRAM files and indexes
- full reference genomes
- large intermediate VCF/GVCF files
- server logs and temporary pileups
- private server paths, credentials, tokens, or personal metadata
