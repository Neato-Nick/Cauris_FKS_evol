version 1.0

# Variant calling pipeline using sample sheet input
# Developed by Nick Cauldron (nicholas_cauldron@brown.edu) and Paul Cao (paul_cao@brown.edu)
# Molecular Microbiology and Immunology; Center for Computation and Visualization
# Brown University

# Implements sample-sheet input for variant calling pipeline
import "fungal-variant-call-gatk4_Part1.wdl" as ubam2gvcf
#import "fungal-variant-call-gatk4_Part2.wdl" as gvcf2vcf
import "fungal-variant-call-gatk4_Part2-GenomicsDB.wdl" as gvcf2vcf

workflow fungal_multisamples {
  input {
    # Config Parameters
    String gatk_docker = "broadinstitute/fungal_gatk4:latest"
    String python2_docker = "broadinstitute/python_fungal:latest"
    String gatk_path = "gatk"
    Int general_mem_size_gb = 30

    # Reference files
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa

    # Input files Part 1
    File samplesheet
    Array[Array[String]] sm_entries = read_tsv(samplesheet) 
    # Input files Part 2
    String run_name

    # Optional parameters
    String? HC_extra_args
    String? gcs_output_dir
    File? gatk_intervals_file
  }

    # Parse sample sheet and scatter tasks
    scatter (sm_entry in sm_entries) {
        call ubam2gvcf.GATK4_BAM_VariantCalling_Workflow {
            input: sample_name = sm_entry[0], input_bam = sm_entry[1],
                   gatk_docker = gatk_docker, gatk_path = gatk_path, general_mem_size_gb = general_mem_size_gb,
                   ref = ref, dict = dict, amb = amb, ann = ann, bwt = bwt, fai = fai, pac = pac, sa = sa,
                   HC_extra_args = HC_extra_args, gcs_output_dir = gcs_output_dir, gatk_intervals_file = gatk_intervals_file
        }
    }

    # Genotype all gvcfs into one vcf
    call gvcf2vcf.GATK4_JointGenotyping_Workflow {
        input: gvcf_files = GATK4_BAM_VariantCalling_Workflow.gvcf, gvcf_index_files = GATK4_BAM_VariantCalling_Workflow.gvcf_index,
               gatk_docker = gatk_docker, gatk_path = gatk_path, general_mem_size_gb = general_mem_size_gb, python2_docker = python2_docker,
               ref = ref, dict = dict, amb = amb, ann = ann, bwt = bwt, fai = fai, pac = pac, sa = sa, run_name = run_name
    }

    # All files expected in output directory
    output {
        # Part 1
        Array[String] gvcf_summary = GATK4_BAM_VariantCalling_Workflow.gvcf_summary
        Array[File] gvcfs = GATK4_BAM_VariantCalling_Workflow.gvcf
        Array[File] gvcf_indexes = GATK4_BAM_VariantCalling_Workflow.gvcf_index
        Array[File] bams = GATK4_BAM_VariantCalling_Workflow.processed_bam
        Array[File] bais = GATK4_BAM_VariantCalling_Workflow.processed_bai
        # Part 2
        File raw_vcf = GATK4_JointGenotyping_Workflow.raw_vcf
        File filtered_vcf = GATK4_JointGenotyping_Workflow.filtered_vcf
        File filtered_snps = GATK4_JointGenotyping_Workflow.filtered_snps
        File filered_indels = GATK4_JointGenotyping_Workflow.filtered_indels
        File costume_filtered_vcf = GATK4_JointGenotyping_Workflow.costume_filtered_vcf
        File costume_filtered_stats = GATK4_JointGenotyping_Workflow.costume_filtered_stats
        File multialelic_stats = GATK4_JointGenotyping_Workflow.multialelic_stats
    }
}
