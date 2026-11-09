# Workflow Development Language scripts

Modified from a WDL previously developed and used by the Fungal Group of the Broad Institute:

* `DownloadFromSRA_Paired.wdl`: Download reads from an SRA accession.
* `paired-fastq-to-unmapped-bam.wdl` and `paired-fastq-to-unmapped-bam.smsheet.wdl`: Convert paired FASTQ reads into an unmapped BAM for compatibility with variant calling workflow
* `fungal-variant-call-gatk4_Part1.wdl`: Per-sample operations of variant calling pipeline: Map reads, format alignments, and genotype sample (HaplotypeCaller)
* `fungal-variant-call-gatk4_Part2-GenomicsDB.wdl`: Pooled operations of variant calling pipeline: Split variants by chromosome, jointly genotype, and filter variants.
* `fungal-variant-call-gatk4_Part1and2_smsheet.wdl`: Run both per-sample and pooled operations for a sample sheet input of samples and BAMs

Developed for this project:

* `magicblast_ubam_smsheet.wdl`: Query SRA accessions from a sample sheet for reads matching a BLAST database.
Usage note: BLAST database must be built separately and added into a Docker configured for Magic-BLAST.
* `refilter_annotate_vcf.wdl`: Filter variants with bcftools, censor genotypes, and annotate with SnpEff.
The genotype censoring script ran within this workflow was previously developed and used by the Fungal Group of the Broad Institute. It is available in the Docker container and accessible here: https://github.com/broadinstitute/broad-fungalgroup/blob/master/scripts/SNPs/filterGatkGenotypes.py

# Python scripts

* `compute_pn_ps_dn_ds.py`: Compute per-site polymorphism (`Pn`, `Ps`, `Pn_indel`, `Pn_indel_fs`) and optionally
per-site fixed-difference (`Dn`, `Ds`, `Dn_indel`, `Dn_indel_fs`) counts from a
SnpEff-annotated VCF and a sample-to-clade population map.
* `nei_gojobori_sites.py`: Nei-Gojobori (1986) synonymous/nonsynonymous site counts per CDS.
Supplemental documentation available in `../docs`

# R scripts

* `plot_pn_ps_nummuts_normalized.R`: Used to parse pN/pS data and recurrent mutation analysis, and visualize lollipop plots
* `identify_gene_variants.R`: Identify amino acid variant haplotype of each sample in a VCF.
Supplemental documentation available in `../docs`
