version 1.0

# Variant calling pipeline for fungal haploid genomes: Part I Single BAM
# Developed by Xiao Li (xiaoli@broadinstitute.org) and Aina Martinez Zurita (amartine@broadinstitute.org)
# Fungal Genomics Group, Infectious Disease and Microbiome Program.
# The Broad Institute of MIT and Harvard


workflow GATK4_BAM_VariantCalling_Workflow {
  input {
    # Config Parameters
    String gatk_docker
    String gatk_path
    Int general_mem_size_gb
    Int cpus_bwa = 4
  
    # Reference files
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa

    # Input Files
    String sample_name
    File input_bam

    # Optional parameters
    String? HC_extra_args
    String? gcs_output_dir
    File? gatk_intervals_file
  }

  #Document Tool
  call MarkIlluminaAdapters {
      input:
        input_bam=input_bam,
        input_sample_name=sample_name,
        gatk_path=gatk_path,
        docker=gatk_docker,
        mem_size_gb=general_mem_size_gb
  }

  #Document tool
  call SamToFastqAllignMerge {
    input:
      input_bam=MarkIlluminaAdapters.bam,
      input_sample_name=sample_name,
      input_unmapped_bam=MarkIlluminaAdapters.qname_bam,
      gatk_path=gatk_path,
      docker=gatk_docker,
      n_cpus=cpus_bwa,

      ref=ref,
      dict=dict,
      amb=amb,
      ann=ann,
      bwt=bwt,
      fai=fai,
      pac=pac,
      sa=sa
  }

  #Document Tool
  call MarkDuplicates{
    input:
      input_sorted_bam=SamToFastqAllignMerge.bam,
      input_sample_name=sample_name,
      mem_size_gb=general_mem_size_gb,
      docker=gatk_docker,
      gatk_path=gatk_path

  }

  call ReorderBam {
    input:
      input_bam=MarkDuplicates.bam,
      input_sample_name=sample_name,

      ref=ref,
      dict=dict,

      mem_size_gb=general_mem_size_gb,
      docker=gatk_docker,
      gatk_path=gatk_path,
      gcs_output_dir = gcs_output_dir
  }

  call HaplotypeCaller {
    input:
      input_bam=ReorderBam.bam,
      input_bam_index=ReorderBam.bai,
      input_sample_name=sample_name,

      ref_dict=dict,
      ref=ref,
      ref_index=fai,

      HC_extra_args=HC_extra_args,
      gatk_intervals_file=gatk_intervals_file,
      gcs_output_dir = gcs_output_dir,

      docker=gatk_docker,
      gatk_path=gatk_path
  }

  # Variant info for sample, for QC and for informative stdout
  call GvcfStats {
    input:
      input_gvcf=HaplotypeCaller.output_gvcf
  }

  output {
    File processed_bam=ReorderBam.bam
    File processed_bai=ReorderBam.bai

    File gvcf=HaplotypeCaller.output_gvcf
    File gvcf_index=HaplotypeCaller.output_gvcf_index
    String gvcf_summary=GvcfStats.bcf_stats_sn
  }

}

### Task Definitions ###

task MarkIlluminaAdapters {
  input {
    File input_bam
    String input_sample_name

    String gatk_path
    String docker
    Int mem_size_gb
  }


  command {
    # e - exit when a command fails
    # u - exit when script tries to use undeclared variables
    # x - trace what gets executed
    # o - exit status of the last command that threw a non-zero exit code is returned
    set -euxo pipefail
    # Queryname sort according to GATK, which is different than samtools querysort
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" SortSam -I ${input_bam} \
    -O ${input_sample_name}_SortSam.bam -SORT_ORDER queryname

    # Mark Illumina adapters
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MarkIlluminaAdapters -I ${input_sample_name}_SortSam.bam \
    -O ${input_sample_name}_markilluminaadapters.bam -M ${input_sample_name}_markilluminaadapters_metrics.txt
  }

  runtime {
    docker: docker
  }

  output {
    File bam = "${input_sample_name}_markilluminaadapters.bam"
    File metrics_adapters = "${input_sample_name}_markilluminaadapters_metrics.txt"
    File qname_bam = "${input_sample_name}_SortSam.bam"
  }
}


task SamToFastqAllignMerge {
  input {
    File input_bam
    String input_sample_name
    File input_unmapped_bam

    String docker
    String gatk_path
    Int mem_size_gb=15
    Int n_cpus=4

    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
  }

  String read_group = "'@RG\\tID:FLOWCELL_${input_sample_name}\\tSM:${input_sample_name}\\tPL:ILLUMINA\\tLB:LIB_${input_sample_name}'"

  command {
    # e - exit when a command fails
    # u - exit when script tries to use undeclared variables
    # x - trace what gets executed
    # o - exit status of the last command that threw a non-zero exit code is returned
    set -euxo pipefail

    #We need to assign RG to all the unmapped reads
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" AddOrReplaceReadGroups \
    -I ${input_unmapped_bam} -O ${input_sample_name}.readgroups_unmapped.bam -ID FLOWCELL_${input_sample_name} \
    -LB LIB_${input_sample_name} -PL ILLUMINA -SM ${input_sample_name} \
    -PU unit1

    #Piped command containing SamToFastq | bwa mem | samtools view | MergeBamAlignment
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" SamToFastq -I ${input_bam} \
    --FASTQ /dev/stdout  -CLIP_ATTR XT  -CLIP_ACT 2  -INTER true -NON_PF true | \
    bwa mem -t ${n_cpus} -M -R ${read_group} -p ${ref} /dev/stdin | samtools view -1 | \
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MergeBamAlignment \
    -ALIGNED /dev/stdin -UNMAPPED ${input_sample_name}.readgroups_unmapped.bam \
    -O ${input_sample_name}.sorted.bam -R ${ref}

  }

  runtime {
    docker: docker
    cpu: n_cpus
  }

  output {
    File bam="${input_sample_name}.sorted.bam"
  }
}


task MarkDuplicates {
  input {
    File input_sorted_bam
    String input_sample_name

    Int mem_size_gb
    String docker
    String gatk_path
  }

  command {
    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" MarkDuplicates \
    -I ${input_sorted_bam} -O ${input_sample_name}.marked_duplicates.bam \
    -M ${input_sample_name}.marked_duplicates.metrics
  }

  output {
    File bam = "${input_sample_name}.marked_duplicates.bam"
    File metrics_duplicates = "${input_sample_name}.marked_duplicates.metrics"
  }

  runtime {
    docker: docker
  }
}


task ReorderBam {
  input {
    File input_bam
    String input_sample_name

    File ref
    File dict

    String? gcs_output_dir

    Int mem_size_gb
    String docker
    String gatk_path
  }

  Int cmd_mem_size_gb = mem_size_gb - 1

  command {
    # reorder bam
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" ReorderSam \
    -I ${input_bam} -O ${input_sample_name}.reordered.bam -SD ${ref}

    # then index
    ${gatk_path} --java-options "-Xmx${cmd_mem_size_gb}G" BuildBamIndex \
    -I ${input_sample_name}.reordered.bam

    gcs_dir="~{default="" gcs_output_dir}"
    if [[ -n "$gcs_dir" ]]; then
      echo "Copying reordered.bam and index to google cloud storage bucket"
      date
      gsutil cp "${input_sample_name}.reordered.bam" "$gcs_dir/bams_imported/" || true
      gsutil cp "${input_sample_name}.reordered.bai" "$gcs_dir/bams_imported/" || true
      echo "Finished copying"
      date
    fi
  }

  output {
    File bam = "${input_sample_name}.reordered.bam"
    File bai = "${input_sample_name}.reordered.bai"
  }

  runtime {
    docker: docker
    maxRetries: 2
  }
}

task HaplotypeCaller {
  input {
    File input_bam
    File input_bam_index
    String input_sample_name
    String gvcf_name = "${input_sample_name}.g.vcf.gz"
    String gvcf_index = "${input_sample_name}.g.vcf.gz.tbi"

    File ref_dict
    File ref
    File ref_index

    String? HC_extra_args
    String? gcs_output_dir
    File? gatk_intervals_file

    Int mem_size_gb=30
    String docker
    String gatk_path
  }

  command {
    gatk_intervals="~{default="" gatk_intervals_file}"
    if [[ -n "$gatk_intervals" ]]; then
      HC_extra_args="${HC_extra_args} -L $gatk_intervals"
    fi

    ${gatk_path} --java-options "-Xmx${mem_size_gb}G" HaplotypeCaller \
    ${HC_extra_args} \
    -R ${ref} -I ${input_bam} -O ${gvcf_name} -ERC GVCF -ploidy 1


    gcs_dir="~{default="" gcs_output_dir}"
    if [[ -n "$gcs_dir" ]]; then
      echo "Copying reordered.bam and index to google cloud storage bucket"
      date
      gsutil cp "${gvcf_name}" "$gcs_dir/gvcfs_imported/" || true
      gsutil cp "${gvcf_index}" "$gcs_dir/gvcfs_imported/" || true
      echo "Finished copying"
      date
    fi

  }

  output {
    File output_gvcf = "${gvcf_name}"
    File output_gvcf_index = "${gvcf_index}"
  }

  runtime {
    docker: docker
  }

}

task GvcfStats {
  input { 
    File input_gvcf
  }

  command {
    module load bcftools && \
    bcftools stats ${input_gvcf} | grep "^SN" 
  }

  output {
    String bcf_stats_sn = read_string(stdout())
  }
}
