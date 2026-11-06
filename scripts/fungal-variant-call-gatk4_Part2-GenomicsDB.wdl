version 1.0

workflow GATK4_JointGenotyping_Workflow {
  input {
    String gatk_docker
    String python2_docker
    String gatk_path
    Int general_mem_size_gb
    Boolean validate_gvcfs = false

    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa

    # Optional parameters
    File? gatk_intervals_file

    # Run config files and strings
    String run_name
    Array[File] gvcf_files
    Array[File] gvcf_index_files
    #Array[String] sample_names
    #File samplesheet
    #Array[Array[String]] sm_entries = read_tsv(samplesheet)
    # Extracting columns into their respective arrays (explicit indexing)
    #Array[String] sample_names = [for row in sm_entries: row[0]]  # First column
    #Array[File] gvcf_files = [for row in sm_entries: row[1]]      # Second column
    #Array[File] gvcf_index_files = [for row in sm_entries: row[2]] # Third column
    
  }

  # Extract chromosome names from the reference .fai file (column 1)
  call ExtractChromosomes {
    input:
      fai = fai
  }

  Array[String] chromosomes = read_lines(ExtractChromosomes.chrom_list)

  # Validate GVCFs (optional, not scattered by chromosome)
  if(validate_gvcfs) {
    scatter (i in range(length(gvcf_files))) {
      call ValidateGVCF {
        input:
          gvcf_file = gvcf_files[i],
          gvcf_index_file = gvcf_index_files[i],
          gatk_path = gatk_path,
          docker = gatk_docker,
          mem_size_gb = general_mem_size_gb,
          ref_fasta = ref,
          ref_fai = fai,
          ref_dict = dict
      }
    }
  }

  # Prepare sample map for GenomicsDBImport (required format: Sample[TAB]File)
  # need gcloud bucket location (prefix gs://) not the mounted location from cromwell (prefix //mnt)
  #Array[String] gvcf_files_gs = gvcf_files
  #call WriteSampleMap {
  #  input:
  #    sample_names = sample_names,
  #    gvcf_files_gs = gvcf_files_gs
  #}

  # Per-chromosome joint calling
  scatter (chrom in chromosomes) {
    call GenomicsDBImport {
      input:
    #    sample_map = WriteSampleMap.sample_map,
    #    sample_names = sample_names,
        gvcf_files = gvcf_files,
        gvcf_index_files = gvcf_index_files,
        interval = chrom,
        workspace_dir_name = "genomicsdb_" + chrom,
        disk_size = 100,
        gatk_path = gatk_path,
        docker = gatk_docker,
        gatk_intervals_file=gatk_intervals_file
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = GenomicsDBImport.output_genomicsdb,
        interval = chrom,
        output_vcf_filename = "genotyped." + chrom + ".vcf.gz",
        ref_fasta = ref,
        ref_fasta_index = fai,
        ref_dict = dict,
        gatk_path = gatk_path,
        gatk_docker = gatk_docker,
        disk_size = 30,
        gatk_intervals_file=gatk_intervals_file
    }
  }

  # Gather all per-chromosome VCFs
  call GatherVcfs {
    input:
      #input_vcfs_fofn = write_lines(GenotypeGVCFs.output_vcf),

      # DF modification 20260507: Commented out input_vcfs_fofn and made the Array[File] type input_vcfs. write_lines was not working correctly.
      input_vcfs = GenotypeGVCFs.output_vcf,
      output_vcf_name = run_name + ".joint.vcf.gz",
      gatk_path = gatk_path,
      docker = gatk_docker,
      disk_size = 30
  }

  # Hard filter on the gathered VCF
  call HardFiltration {
    input:
      ref = ref,
      ref_dict = dict,
      ref_index = fai,
      vcf = GatherVcfs.output_vcf,
      vcf_index = GatherVcfs.output_vcf_index,
      output_filename = run_name + ".hard_filtered.vcf.gz",
      docker = gatk_docker,
      gatk_path = gatk_path,
      mem_size_gb = general_mem_size_gb
  }

  call MakeReadableVCF {
    input:
      vcf_file = HardFiltration.all_filtered_variants,
      output_filename = run_name,
      docker = gatk_docker,
      mem_size_gb = general_mem_size_gb
  }

  call CostumeVCFFilter {
    input:
      vcf_file = MakeReadableVCF.vcf_output,
      output_filename = run_name,
      docker = python2_docker,
      mem_size_gb = general_mem_size_gb
  }

  output {
    # File sample_map = WriteSampleMap.sample_map
    File raw_vcf = GatherVcfs.output_vcf
    File filtered_vcf = HardFiltration.all_filtered_variants
    File filtered_snps = HardFiltration.snps
    File filtered_indels = HardFiltration.indels
    File costume_filtered_vcf = CostumeVCFFilter.vcf_filtered
    File costume_filtered_stats = CostumeVCFFilter.vcf_filter_stats
    File multialelic_stats = CostumeVCFFilter.multialelic_stats
  }
}

task ExtractChromosomes {
  input {
    File fai
  }
  command <<<
    cut -f1 ~{fai} | uniq > chrom.list
  >>>
  output {
    File chrom_list = "chrom.list"
  }
  runtime {
    memory: "512 MB"
    docker: "ubuntu:latest"
  }
}

task WriteSampleMap {
  input {
    Array[String] sample_names
    Array[String] gvcf_files_gs
  }

  command <<<
    paste <(printf "%s\n" ~{sep=' ' sample_names}) <(printf "%s\n" ~{sep=' ' gvcf_files_gs}) > sample_map.txt
  >>>
  output {
    File sample_map = "sample_map.txt"
  }
  runtime {
    docker: "ubuntu:latest"
    memory: "512 MB"
  }
}

task ValidateGVCF {
  input {
    File gvcf_file
    File gvcf_index_file
    String gatk_path
    String docker
    Int mem_size_gb
    File ref_fasta
    File ref_fai
    File ref_dict
  }
  Int cmd_mem_size_gb = mem_size_gb - 1

  command <<<
    echo "Validating GVCF: ~{gvcf_file} with its index ~{gvcf_index_file}"
    ~{gatk_path} --java-options "-Xmx~{cmd_mem_size_gb}G" ValidateVariants \
      -V ~{gvcf_file} -gvcf -R ~{ref_fasta}
    if [ $? -ne 0 ]; then
      echo "GVCF validation returned a non-zero exit status."
      exit 1
    else
      touch "validated.ok"
    fi
  >>>
  output {
    File gvcf_validation_markers = "validated.ok"
  }
  runtime {
    max_retries: 2
    docker: docker
    memory: mem_size_gb + " GB"
  }
}

task GenomicsDBImport {
  input {
    #File sample_map
    String interval
    String workspace_dir_name
    Int disk_size
    String gatk_path
    String docker
    File? gatk_intervals_file
    #Array[String] sample_names
    Array[File] gvcf_files
    Array[File] gvcf_index_files
    # --sample-name-map ~{sample_map}
  }

  command <<<
    set -euxo pipefail
    gatk_intervals="~{default="" gatk_intervals_file}"
    if [[ -n "$gatk_intervals" ]]; then
      intervals="intervals_in_chr.bed"
      grep ~{interval} $gatk_intervals > $intervals
    else
      intervals="interval_chr.list"
      echo ~{interval} > $intervals
    fi

    ~{gatk_path} --java-options "-Xmx8g" GenomicsDBImport \
      --genomicsdb-workspace-path ~{workspace_dir_name} \
      --batch-size 50 \
      -V ~{sep=" -V " gvcf_files} \
      -L "$intervals"
    tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}
  >>>
  output {
    File output_genomicsdb = "~{workspace_dir_name}.tar"
  }
  runtime {
    memory: "32 GB"
    disks: "local-disk ~{disk_size} HDD"
    docker: docker
  }
}

task GenotypeGVCFs {
  input {
    File workspace_tar
    String interval
    String output_vcf_filename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String gatk_path
    String gatk_docker
    Int disk_size
    File? gatk_intervals_file
  }
  command <<<
    set -euxo pipefail
    tar -xf ~{workspace_tar}
    WORKSPACE=$( basename ~{workspace_tar} .tar)
    gatk_intervals="~{default="" gatk_intervals_file}"
    if [[ -n "$gatk_intervals" ]]; then
      grep ~{interval} $gatk_intervals > intervals_in_chr.bed
      intervals="intervals_in_chr.bed"
    else
      echo ~{interval} > "interval_chr.list"
      intervals="interval_chr.list"
    fi

    ~{gatk_path} --java-options "-Xmx8g" GenotypeGVCFs \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename} \
      -V gendb://$WORKSPACE \
      -L "$intervals"
  >>>
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
  runtime {
    memory: "32 GB"
    disks: "local-disk ~{disk_size} HDD"
    docker: gatk_docker
  }
}

task GatherVcfs {
  input {
    #File input_vcfs_fofn

    # DF modification 20260507 - Commented out intput_vcfs_fofn and made the Array[File] type input_vcfs
    Array[File] input_vcfs
    String output_vcf_name
    String gatk_path
    String docker
    Int disk_size
  }
  command <<<
    # mv ~ {input_vcfs_fofn} inputs.list

    # DF modification 20260507 - Commented out the mv command above and added a space between ~ and {. Added in the echo command to make inputs.list from input_vcfs
    echo -e "~{sep='\n' input_vcfs}" > inputs.list

    ~{gatk_path} --java-options "-Xmx6g" GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input inputs.list \
      --output ~{output_vcf_name}

    tabix ~{output_vcf_name}
  >>>
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
  runtime {
    memory: "7 GB"
    disks: "local-disk ~{disk_size} HDD"
    docker: docker
  }
}

task HardFiltration {
  input {
    File ref
    File vcf
    File vcf_index
    File ref_dict
    File ref_index
    String output_filename
    String docker
    String gatk_path
    Int mem_size_gb
  }
  Int cmd_mem_size_gb = mem_size_gb - 1

  command <<<
    # Select SNPs
    ~{gatk_path} --java-options "-Xmx~{cmd_mem_size_gb}G" SelectVariants \
      -V ~{vcf} -R ~{ref} -select-type SNP -O snps.vcf.gz

    # Select Indels
    ~{gatk_path} --java-options "-Xmx~{cmd_mem_size_gb}G" SelectVariants \
      -V ~{vcf} -R ~{ref} -select-type INDEL -O indels.vcf.gz

    # Filter SNPs
    ~{gatk_path} --java-options "-Xmx~{cmd_mem_size_gb}G" VariantFiltration \
      -V snps.vcf.gz \
      -R ~{ref} \
      -filter "QD < 20.0" --filter-name "QD20" \
      -filter "QUAL < 30.0" --filter-name "QUAL30" \
      -filter "SOR > 3.0" --filter-name "SOR3" \
      -filter "FS > 60.0" --filter-name "FS60" \
      -filter "MQ < 40.0" --filter-name "MQ40" \
      -O snps_filtered.vcf.gz

    # Filter Indels
    ~{gatk_path} --java-options "-Xmx~{cmd_mem_size_gb}G" VariantFiltration \
      -V indels.vcf.gz \
      -R ~{ref} \
      -filter "QD < 20.0" --filter-name "QD20" \
      -filter "QUAL < 30.0" --filter-name "QUAL30" \
      -filter "FS > 200.0" --filter-name "FS200" \
      -O indels_filtered.vcf.gz

    # Merge filtered VCFs
    ~{gatk_path} --java-options "-Xmx~{cmd_mem_size_gb}G" MergeVcfs \
      -R ~{ref} \
      -I snps_filtered.vcf.gz \
      -I indels_filtered.vcf.gz \
      -O ~{output_filename}
  >>>
  output {
    File indels = "indels_filtered.vcf.gz"
    File snps = "snps_filtered.vcf.gz"
    File all_filtered_variants = output_filename
  }
  runtime {
    docker: docker
    memory: mem_size_gb + " GB"
  }
}

task MakeReadableVCF {
  input {
    File vcf_file
    String output_filename
    String docker
    Int mem_size_gb
  }
  command <<<
    bcftools view ~{vcf_file} > ~{output_filename}.nonbinary.vcf
  >>>
  output {
    File vcf_output = output_filename + ".nonbinary.vcf"
  }
  runtime {
    docker: docker
    memory: mem_size_gb + " GB"
  }
}

task CostumeVCFFilter {
  input {
    File vcf_file
    String output_filename
    String docker
    Int mem_size_gb
  }
  command <<<
    ln -s ~{vcf_file} out.vcf
    python /broad-fungalgroup/scripts/SNPs/filterGatkGenotypes.py --min_GQ 50 \
      --min_percent_alt_in_AD 0.8 --min_total_DP 10 out.vcf \
      > ~{output_filename}.filtered_SNPs_GQ50_AD08_DP10.vcf 2> ~{output_filename}.variant_qc_genotype_filter.tsv

    mv out.vcf.multiallelic.txt ~{output_filename}.multiallelic.txt
  >>>
  output {
    File vcf_filtered = output_filename + ".filtered_SNPs_GQ50_AD08_DP10.vcf"
    File vcf_filter_stats = output_filename + ".variant_qc_genotype_filter.tsv"
    File multialelic_stats = output_filename + ".multiallelic.txt"
  }
  runtime {
    docker: docker
    memory: mem_size_gb + " GB"
  }
}
