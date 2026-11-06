# Workflow Development Language scripts

* `DownloadFromSRA_Paired.wdl`: Download reads from an SRA accession
* `magicblast_ubam_smsheet.wdl`: Query SRA accessions from a sample sheet for reads matching a BLAST database.
Usage note: BLAST database must be built separately and added into a Docker configured for Magic-BLAST.
* `paired-fastq-to-unmapped-bam.wdl` and `paired-fastq-to-unmapped-bam.smsheet.wdl`: Convert paired FASTQ reads into an unmapped BAM for compatibility with variant calling workflow
* `fungal-variant-call-gatk4_Part1.wdl`: Per-sample operations of variant calling pipeline: Map reads, format alignments, and genotype sample (HaplotypeCaller)
* `fungal-variant-call-gatk4_Part2-GenomicsDB.wdl`: Pooled operations of variant calling pipeline: Split variants by chromosome, jointly genotype, and filter variants.
* `fungal-variant-call-gatk4_Part1and2_smsheet.wdl`: Run both per-sample and pooled operations for a sample sheet input of samples and BAMs
* `refilter_annotate_vcf.wdl`: Filter variants with bcftools, censor genotypes, and annotate with SnpEff.
