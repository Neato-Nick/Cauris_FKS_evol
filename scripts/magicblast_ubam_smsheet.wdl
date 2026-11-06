version 1.0

# MagicBlast -> RevertSam WDL using read_tsv(sample_sheet) and scatter
# Defaults chosen to match your environment; runtime singularity_image defaults to your local .sif path.
#
# Usage:
# - Provide sample_sheet (TSV two columns: sample_name<TAB>SRR) as input.
# - The workflow will read_tsv(sample_sheet) into sm_entries and scatter over rows.

workflow MagicBlast_From_SampleSheet {
  input {
    # Docker image name (used if singularity_image not set or not found on node)
    String docker = "ncauldronbroad/gatk_ncbi:latest"

    # Path to gatk inside the container (if 'gatk' is on PATH in the image, use "gatk")
    String gatk_path = "gatk"

    # TSV sample sheet (two columns: sample_name<TAB>SRR), no header assumed.
    File sample_sheet

    # Parse the TSV into an array of rows; each row is an Array[String] of columns
    Array[Array[String]] sm_entries = read_tsv(sample_sheet)

    # MagicBlast DB prefix passed to -db (no default — user should provide)
    String ref_db

    # Output prefix used in filenames
    String out_prefix = "magicblast_ubam"

    # Resource tuning
    Int num_threads = 1
    Int mem_gb = 8
  }

  # Scatter across parsed sample entries (sm_entries[i][0] = sample_name, sm_entries[i][1] = SRR)
  scatter (i in range(length(sm_entries))) {
    call RunMagicBlast {
      input:
        sample_name = sm_entries[i][0],
        SRR = sm_entries[i][1],
        ref_db = ref_db,
        out_prefix = out_prefix,
        num_threads = num_threads,
        mem_gb = mem_gb,
        docker = docker,
        gatk_path = gatk_path
    }
  }

  output {
    Array[File] unmapped_bams = RunMagicBlast.unmapped_bam
    Array[String] task_stdouts = RunMagicBlast.task_stdout
  }
}

task RunMagicBlast {
  input {
    String SRR
    String ref_db
    String sample_name
    String out_prefix
    Int num_threads = 1
    Int mem_gb = 8
    String docker
    String gatk_path = "gatk"
  }

  command {
    set -euxo pipefail

    # Sanitize possible CR characters from read_tsv
    sample_name="$(echo -n '${sample_name}' | tr -d '\r')"
    SRR="$(echo -n '${SRR}' | tr -d '\r')"

    mkdir -p unmapped_bams

    echo "Running magicblast for ${SRR} against db ${ref_db}"
    CMD="magicblast \
      -no_unaligned \
      -num_threads ${num_threads} \
      -splice F \
      -sra ${SRR} \
      -db ${ref_db} > unmapped_bams/${out_prefix}.${sample_name}.paired.sam"
    echo $CMD
    eval $CMD
    date

    echo "Condensing sam to bam"
    CMD="samtools view -f 0x2 \
      -o unmapped_bams/${out_prefix}.${sample_name}.paired.bam unmapped_bams/${out_prefix}.${sample_name}.paired.sam"
    echo $CMD
    eval $CMD
    date

    echo
    echo "magicblast -> samtools completed, produced: unmapped_bams/${out_prefix}.${sample_name}.paired.bam"

    echo "Running ${gatk_path} RevertSam to create unmapped BAM"
    ${gatk_path} RevertSam \
      --SANITIZE true \
      -ALIAS ${sample_name} \
      -I unmapped_bams/${out_prefix}.${sample_name}.paired.bam \
      -O unmapped_bams/${out_prefix}.${sample_name}.unmapped.bam

    date

    # If magicblast produced real output, remove the intermediate paired.bam
    if [[ -s "unmapped_bams/${out_prefix}.${sample_name}.paired.bam" ]]; then
      rm -f unmapped_bams/${out_prefix}.${sample_name}.paired.bam
      rm -f unmapped_bams/${out_prefix}.${sample_name}.paired.sam
    fi

    echo
    echo "Done cleaning up files. Final output bam of unaligned-reads:"
    ls -lh unmapped_bams/${out_prefix}.${sample_name}.unmapped.bam || true
  }

  runtime {
    # These runtime attributes will be passed through to Cromwell and used by the submit-docker wrapper.
    docker: docker

    # resource hints; Cromwell will convert these into the slurm request via the submit config
    cpu: num_threads
    memory: "${mem_gb} GB"
  }

  output {
    File unmapped_bam = "unmapped_bams/${out_prefix}.${sample_name}.unmapped.bam"
    String task_stdout = read_string(stdout())
  }
}
