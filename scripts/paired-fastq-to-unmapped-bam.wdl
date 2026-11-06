version 1.0
##Copyright Broad Institute, 2018
## 
## This WDL converts paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (one file per orientation)
## - The following metada descriptors per sample:
##  - readgroup
##  - sample_name
##  - library_name
##  - platform_unit
##  - run_date
##  - platform_name
##  - sequecing_center
##
## Outputs :
## - Set of unmapped BAMs, one per read group
## - File of a list of the generated unmapped BAMs
##
## Cromwell version support 
## - Successfully tested on v47
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow ConvertPairedFastQsToUnmappedBamWf {
  input {
    String sample_name 
    String? SRA_ID
    Array[File] fastqs_1 
    Array[File] fastqs_2
    String readgroup_name = sample_name
    String? library_name = SRA_ID
    String platform_name = "ILLUMINA"

    #String gatk_docker = "broadinstitute/fungal_gatk4:latest"
    #String gatk_path = "/gatk/gatk"
    String gatk_docker = "broadinstitute/fungal_gatk4:gatk4.1.8.1"
    String gatk_path = "gatk"
    #String gatk_docker = "broadinstitute/fungal_gatk4:gatk4.2.4.1"
    #String gatk_path = "/gatk/gatk"
    #String gatk_path = "gatk"
    String? gcs_output_dir
    Array[File]? SRR_succeeds
    
    # Sometimes the output is larger than the input, or a task can spill to disk.
    # In these cases we need to account for the input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
    Float disk_multiplier = 2.5
  }
    String ubam_list_name = sample_name

#  String platform_unit = if infer_runinfo then GetSRAInfo.platform_unit else "unknown"
#  String platform_model = if infer_runinfo then GetSRAInfo.platform_model else "unknown"
#  String run_date = if infer_runinfo then GetSRAInfo.platform_model else "unknown"
#  String sequencing_center = if infer_runinfo then GetSRAInfo.platform_model else "unknown"

  # Convert pairs of FASTQs to uBAMs
  Int n_parts = length(fastqs_1)
  Array[Int] i_parts = range(n_parts)
  if (n_parts > 1) {
    scatter (fastq_pair in zip(zip(fastqs_1, fastqs_2), i_parts)) {
      call MultiPairedFastQsToUnmappedBAM {
        input:
          sample_name = sample_name,
          fastq_1 = fastq_pair.left.left,
          fastq_2 = fastq_pair.left.right,
          part_num = fastq_pair.right,
          readgroup_name = readgroup_name,
          library_name = library_name,
          platform_name = platform_name,
          gatk_path = gatk_path,
          docker = gatk_docker,
          disk_multiplier = disk_multiplier
      }
    }
    Array[File] merge_ubams = MultiPairedFastQsToUnmappedBAM.output_unmapped_bam
    Int num_ubams_to_merge = length(merge_ubams)
    call MergeUnmappedBAMs {
      input:
        sample_name = sample_name,
        ubams = merge_ubams,
        n_ubams = num_ubams_to_merge,
        docker = gatk_docker
    }
    File merged_ubam = MergeUnmappedBAMs.merged_ubam
    # Create a file with the generated ubam
    call CreateFoFN as CreateFoFNMerged {
        input:
          ubam = merged_ubam,
          fofn_name = ubam_list_name + ".ubam",
          gcs_output_dir = gcs_output_dir
    }
  }
  
  if (n_parts == 1) {
      call SinglePairedFastQsToUnmappedBAM {
        input:
          sample_name = sample_name,
          fastq_1 = select_first(fastqs_1),
          fastq_2 = select_first(fastqs_2),
          part_num = 0,
          readgroup_name = readgroup_name,
          library_name = library_name,
          platform_name = platform_name,
          gatk_path = gatk_path,
          docker = gatk_docker,
          disk_multiplier = disk_multiplier
      }

    File unmerged_ubam = SinglePairedFastQsToUnmappedBAM.output_unmerged_ubam
    # Create a file with the generated ubam
    call CreateFoFN as CreateFoFNUnmerged {
        input:
          ubam = unmerged_ubam,
          fofn_name = ubam_list_name + ".ubam",
          gcs_output_dir = gcs_output_dir
    }
  }
  
  # Outputs that will be retained when execution is complete
  output {
    File? output_unmapped_bam_merged = merged_ubam
    File? output_unmapped_bam_unmerged = unmerged_ubam
    File? unmapped_bam_list_merged = CreateFoFNMerged.fofn_list
    File? unmapped_bam_list_unmerged = CreateFoFNUnmerged.fofn_list
  }
}

# TASK DEFINITIONS

# Get SRA submission info
task GetSRAInfo {
  input {
    String? SRA_ID
  }

  command {
    # PM Platform Model
    PM=$(esearch -db sra -query ~{SRA_ID} | efetch -format xml | xtract -pattern PLATFORM -element INSTRUMENT_MODEL)
    # Platform unit too challenging to infer - skip
    PU=$(echo ~{SRA_ID})

    grep --help
    printf "$PM" > pm.txt
    printf "$PU" > pu.txt
  }

  runtime {
    maxRetries: 2
    docker: "biocontainers/ncbi-entrez-direct:v10.9.20190219ds-1b10-deb_cv1"
  }

  output {
    String platform_unit = read_string("pu.txt")
    String platform_model = read_string("pm.txt")
  }
}

# Convert a pair of FASTQs to uBAM
task MultiPairedFastQsToUnmappedBAM {
  input {
    # Command parameters
    String sample_name
    File? fastq_1
    File? fastq_2
    String readgroup_name
    String? library_name
    String platform_name
    String gatk_path
    Int part_num
    
    # Runtime parameters
    Int addtional_disk_space_gb = 10
    Int machine_mem_gb = 12
    String docker
    Float disk_multiplier
  }
    Int command_mem_gb = machine_mem_gb - 1
    Float fastq_size = size(fastq_1, "GB") + size(fastq_2, "GB")
    Int disk_space_gb = ceil(fastq_size + (fastq_size * disk_multiplier ) + addtional_disk_space_gb)
  command {
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    FastqToSam \
    --FASTQ ~{fastq_1} \
    --FASTQ2 ~{fastq_2} \
    --OUTPUT ~{readgroup_name}.part~{part_num}.unmapped.bam \
    --READ_GROUP_NAME ~{readgroup_name} \
    --SAMPLE_NAME ~{sample_name} \
    --LIBRARY_NAME ~{library_name} \
    --PLATFORM ~{platform_name}
  }
  runtime {
    maxRetries: 3
    docker: docker
  }
  output {
    File output_unmapped_bam = "~{readgroup_name}.part~{part_num}.unmapped.bam"
  }
}

# Convert a pair of FASTQs to uBAM
task SinglePairedFastQsToUnmappedBAM {
  input {
    # Command parameters
    String sample_name
    File? fastq_1
    File? fastq_2
    String readgroup_name
    String? library_name
    String platform_name
    String gatk_path
    Int part_num
    
    # Runtime parameters
    Int addtional_disk_space_gb = 10
    Int machine_mem_gb = 12
    String docker
    Float disk_multiplier
  }
    Int command_mem_gb = machine_mem_gb - 1
    Float fastq_size = size(fastq_1, "GB") + size(fastq_2, "GB")
    Int disk_space_gb = ceil(fastq_size + (fastq_size * disk_multiplier ) + addtional_disk_space_gb)
  command {
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    FastqToSam \
    --FASTQ ~{fastq_1} \
    --FASTQ2 ~{fastq_2} \
    --OUTPUT ~{readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ~{readgroup_name} \
    --SAMPLE_NAME ~{sample_name} \
    --LIBRARY_NAME ~{library_name} \
    --PLATFORM ~{platform_name}
  }
  runtime {
    maxRetries: 3
    docker: docker
  }
  output {
    File output_unmerged_ubam = "~{readgroup_name}.unmapped.bam"
  }
}
# Creats a file of file names of the uBAM, which is a text file with each row having the path to the file.
# In this case there will only be one file path in the txt file but this format is used by 
# the pre-processing for variant discovery workflow. 
task CreateFoFN {
  input {
    # Command parameters
    File ubam
    String fofn_name
    String? gcs_output_dir
  }
  command {
    echo ~{ubam} > ~{fofn_name}.list
    gcs_dir="~{default="" gcs_output_dir}"
    if [[ -n "$gcs_dir" ]]; then
      gsutil cp ~{ubam} ~{gcs_output_dir}/ || true
    fi  
  }
  output {
    File fofn_list = "~{fofn_name}.list"
  }
  runtime {
    maxRetries: 2
    docker: "quay.io/broad-long-read-pipelines/lr-cloud-downloader:0.2.3"
  }
}

# Merges multiple unmapped bams from each sample into just one unmapped bam
task MergeUnmappedBAMs {
  input {
    # Command parameters
    String sample_name
    Array[File] ubams
    Int n_ubams
    String docker
  }
  
  # Merge unmapped bams if there are at least 2;
  # otherwise just link the one
  command {
    if [ ~{n_ubams} -gt 1 ]; then
        printf "Merging multiple unmapped bam files (n=~{n_ubams})\n"
        samtools merge ~{sample_name}.unmapped.bam ~{sep=" " ubams} # for samtools 1.10, -o for output file causes error
    else
        printf "Something went wrong; expected multiple unmapped bams but got only 1 (n=~{n_ubams})\n"
    fi
  }
  output {
    File merged_ubam = "~{sample_name}.unmapped.bam"
  }
  runtime {
    docker: docker
  }
}
