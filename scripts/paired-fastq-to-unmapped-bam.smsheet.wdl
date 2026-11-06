version 1.0

# Enables use of sample sheet input for fastq to unmapped BAM workflow

import "paired-fastq-to-unmapped-bam.wdl" as fq2ubam

workflow Fq2UbamFromSampleSheet {
  input {
    File sample_sheet_tsv

    String? gcs_output_dir
    String platform_name = "ILLUMINA"

    String gatk_docker = "broadinstitute/fungal_gatk4:gatk4.1.8.1"
    String gatk_path = "gatk"

    Float disk_multiplier = 2.5
  }

  Array[Array[String]] rows = read_tsv(sample_sheet_tsv)

  scatter (row in rows) {
    String isolate = row[0]
    File fastq1 = row[1]
    File fastq2 = row[2]

    call fq2ubam.ConvertPairedFastQsToUnmappedBamWf as fq2ubam_one {
      input:
        sample_name = isolate,
        library_name = isolate,
        fastqs_1 = [fastq1],
        fastqs_2 = [fastq2],

        platform_name = platform_name,
        gcs_output_dir = gcs_output_dir,

        gatk_docker = gatk_docker,
        gatk_path = gatk_path,

        disk_multiplier = disk_multiplier
    }
  }

  output {
    Array[File?] output_unmapped_bam_merged = fq2ubam_one.output_unmapped_bam_merged
    Array[File?] output_unmapped_bam_unmerged = fq2ubam_one.output_unmapped_bam_unmerged
    Array[File?] unmapped_bam_list_merged = fq2ubam_one.unmapped_bam_list_merged
    Array[File?] unmapped_bam_list_unmerged = fq2ubam_one.unmapped_bam_list_unmerged
  }
}
