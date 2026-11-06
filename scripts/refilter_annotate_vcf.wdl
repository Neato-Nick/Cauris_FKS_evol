version 1.0

workflow variant_processing {

  input {
    File input_vcf           # Input VCF (raw, unfiltered, e.g. genotyped_gvcfs.vcf.gz)
    String output_basename   # Basename for input VCF file (e.g., "genotyped_gvcfs")
    String scripts_docker = "ncauldronbroad/refilter_vcf:latest"
    String snpeff_config = "/app/local/snpEff-5.2a/snpEff.config"
    File transcript_list = "/app/local/snpEff-5.2a/data/Cauris_B11205_mt/Cand_auris_B11205.add_mt.mRNA.list"
    String snpeff_ref = "Cauris_B11205_mt"
    Array[Int] min_total_DP_array = [10, 4, 6, 8]
    Int min_GQ=50
    Float min_AD=0.8
    Int mem_gb = 32
    String snpeff_jar="/app/local/snpEff-5.2a/snpEff.jar"
    String snpsift_jar="/app/local/snpEff-5.2a/SnpSift.jar"
    String snpeff_data_dir="/app/local/snpEff-5.2a/data"
  }

  call SNPFiltrationBcftools {
    input:
      vcf=input_vcf,
      vcf_basename=output_basename,
      docker=scripts_docker
  }

  call INDELFiltrationBcftools {
    input:
      vcf=input_vcf,
      vcf_basename=output_basename,
      docker=scripts_docker
  }

  call MergeCleanBcftools {
    input:
      snps_vcf=SNPFiltrationBcftools.snps_pass_vcf,
      snps_vcf_index=SNPFiltrationBcftools.snps_pass_vcf_index,
      indels_vcf=INDELFiltrationBcftools.indels_pass_vcf,
      indels_vcf_index=INDELFiltrationBcftools.indels_pass_vcf_index,
      vcf_basename=output_basename,
      docker=scripts_docker
  }

  call SnpEffAnnotate {
    input:
      java_path="/usr/bin/java",
      snpeff_jar=snpeff_jar,
      snpeff_config = snpeff_config,
      snpeff_data_dir=snpeff_data_dir,
      ref=snpeff_ref,
      vcf=MergeCleanBcftools.cleaned_vcf,
      vcf_index=MergeCleanBcftools.cleaned_vcf_index,
      vcf_name=output_basename,
      transcript_list=transcript_list,
      docker=scripts_docker,
      mem_gb=mem_gb
  }

  scatter (min_total_DP in min_total_DP_array) {
    String name_stem_out = "~{output_basename}.filtered_SNPs_GQ~{min_GQ}_AD~{sub(min_AD, '\\.', '')}_DP~{min_total_DP}"
    call RefilterVCF {
      input:
        filter_script="/app/local/filterGatkGenotypes_stdin.py",
        vcf_f_in=SnpEffAnnotate.vcf_snpeff_gz,
        vcf_f_in_tbi=SnpEffAnnotate.vcf_snpeff_gz_tbi,
        min_GQ=min_GQ,
        min_AD=min_AD,
        min_tot_DP=min_total_DP,
        name_stem_out=name_stem_out,
        vcf_name=output_basename,
        python_path="python3",
        docker=scripts_docker,
        mem_gb = mem_gb,
        java_path="/usr/bin/java",
        snpeff_jar=snpeff_jar,
        snpsift_jar=snpsift_jar,
        snpeff_config = snpeff_config,
        snpeff_data_dir=snpeff_data_dir,
        ref=snpeff_ref,
        transcript_list=transcript_list
    }
  }

  output {
    File filtered_annot_vcf = SnpEffAnnotate.vcf_snpeff_gz
    File snpeff_summary = SnpEffAnnotate.snpeff_summary
    Array[File] censored_vcfs = RefilterVCF.filtered_vcf
    Array[File] censored_stats = RefilterVCF.filter_tsv
  }
}

task SNPFiltrationBcftools {
  input {
    File vcf
    String vcf_basename
    String docker
  }
  command <<<
    set -euo pipefail
    tabix -f ~{vcf}
    # Task 1: SNP hard filter
    bcftools norm -Ou -m- --atom-overlaps '*' --atomize ~{vcf} | \
      bcftools view -Ou -a | \
      bcftools view -Ou -i 'TYPE="SNP"' | \
      bcftools filter -s 'QD20' -e 'INFO/QD < 20.0' -Ou | \
      bcftools filter -s 'QUAL30' -e 'QUAL < 30.0' -m+ -Ou | \
      bcftools filter -s 'SOR3' -e 'INFO/SOR > 3.0' -m+ -Ou | \
      bcftools filter -s 'FS60' -e 'INFO/FS > 60.0' -m+ -Ou | \
      bcftools filter -s 'MQ40' -e 'INFO/MQ < 40.0' -m+ -Oz -o snps.bcftools_filter.vcf.gz
    tabix -f snps.bcftools_filter.vcf.gz
    bcftools view -i "FILTER='PASS'" -Ou snps.bcftools_filter.vcf.gz | \
      bcftools norm -m+ --atom-overlaps '*' -Oz -o snps.bcftools_filter.pass.vcf.gz
    tabix -f snps.bcftools_filter.pass.vcf.gz
  >>>
  output {
    File snps_pass_vcf = "snps.bcftools_filter.pass.vcf.gz"
    File snps_pass_vcf_index = "snps.bcftools_filter.pass.vcf.gz.tbi"
  }
  runtime {
    docker: docker
    memory: "32G"
    cpu: 4
  }
}

task INDELFiltrationBcftools {
  input {
    File vcf
    String vcf_basename
    String docker
  }
  command <<<
    set -euo pipefail
    tabix -f ~{vcf}
    # Task 2: INDEL hard filter
    bcftools norm -Ou -m- --atom-overlaps '*' --atomize ~{vcf} | \
      bcftools view -Ou -a | \
      bcftools view -Ou -i 'TYPE!="SNP"' | \
      bcftools filter -s 'QD20' -e 'INFO/QD < 20.0' -Ou | \
      bcftools filter -s 'QUAL30' -e 'QUAL < 30.0' -m+ -Ou | \
      bcftools filter -s 'FS200' -e 'INFO/FS > 200.0' -m+ -Oz -o indels.bcftools_filter.vcf.gz
    tabix -f indels.bcftools_filter.vcf.gz
    bcftools view -Ou -i "FILTER='PASS'" indels.bcftools_filter.vcf.gz | \
      bcftools norm -m+ --atom-overlaps '*' -Oz -o indels.bcftools_filter.pass.vcf.gz
    tabix -f indels.bcftools_filter.pass.vcf.gz
  >>>
  output {
    File indels_pass_vcf = "indels.bcftools_filter.pass.vcf.gz"
    File indels_pass_vcf_index = "indels.bcftools_filter.pass.vcf.gz.tbi"
  }
  runtime {
    docker: docker
    memory: "32G"
    cpu: 4
  }
}

task MergeCleanBcftools {
  input {
    File snps_vcf
    File snps_vcf_index
    File indels_vcf
    File indels_vcf_index
    String vcf_basename
    String docker
  }
  command <<<
    
    CMD="ls /app/*"
    echo $CMD
    eval $CMD
    
    set -euo pipefail

    # Task 3: Merge, remove spanning indels, clean
    bcftools concat --allow-overlaps -Ou ~{snps_vcf} ~{indels_vcf} | \
      bcftools norm -m+any --atom-overlaps '*' -Ou | \
      bcftools view -Oz -o ~{vcf_basename}.bcft_hard_filtered.pass.vcf.gz -e "N_ALT=1 & ALT='*'"
    tabix -f ~{vcf_basename}.bcft_hard_filtered.pass.vcf.gz

    bcftools norm -Oz -o ~{vcf_basename}.bcft_hard_filtered.pass.norm.vcf.gz -m- --atom-overlaps '*' --atomize ~{vcf_basename}.bcft_hard_filtered.pass.vcf.gz
    tabix ~{vcf_basename}.bcft_hard_filtered.pass.norm.vcf.gz

    python3 /app/local/filter_spanning_indels_v1.py -v ~{vcf_basename}.bcft_hard_filtered.pass.norm.vcf.gz -o ~{vcf_basename}.bcft_hard_filtered.pass.rm_stars.vcf.gz
    tabix ~{vcf_basename}.bcft_hard_filtered.pass.rm_stars.vcf.gz

    bcftools view -Ou -a ~{vcf_basename}.bcft_hard_filtered.pass.rm_stars.vcf.gz | \
      bcftools view -Ou --min-alleles 2 | \
      bcftools norm -m+any --atom-overlaps '*' -Ou | \
      bcftools view -e 'N_ALT=1 & ALT="*"' -Ou | \
      bcftools +fill-tags -Oz -o ~{vcf_basename}.bcft_hard_filtered.pass.rm_stars.clean.vcf.gz -- -t AN,AC,AF,MAF,F_MISSING
    mv ~{vcf_basename}.bcft_hard_filtered.pass.rm_stars.clean.vcf.gz ~{vcf_basename}.vcf.gz
    tabix -f ~{vcf_basename}.vcf.gz
  >>>
  output {
    File cleaned_vcf = "~{vcf_basename}.vcf.gz"
    File cleaned_vcf_index = "~{vcf_basename}.vcf.gz.tbi"
  }
  runtime {
    docker: docker
    memory: "32G"
    cpu: 4
  }
}

task SnpEffAnnotate {
  input {
    String java_path
    String snpeff_jar
    String snpeff_config
    String snpeff_data_dir
    String ref
    File vcf
    String vcf_name
    File vcf_index
    File transcript_list
    String docker
    Int mem_gb
  }
  command <<<
    set -euo pipefail
    export SNPEFF_HOME=$(dirname ~{snpeff_jar})
    mkdir stats
    ~{java_path} -Xmx~{mem_gb}G -jar ~{snpeff_jar} -i vcf -v \
      -s stats/~{vcf_name}.snpEff_summary.html \
      -c $SNPEFF_HOME/snpEff.config \
      -onlyTr ~{transcript_list} \
      -nodownload \
      -datadir ~{snpeff_data_dir} \
      ~{ref} ~{vcf} > ~{vcf_name}.snpEff.vcf
    bgzip ~{vcf_name}.snpEff.vcf
    tabix -f ~{vcf_name}.snpEff.vcf.gz
  >>>
  output {
    File vcf_snpeff_gz = "~{vcf_name}.snpEff.vcf.gz"
    File vcf_snpeff_gz_tbi = "~{vcf_name}.snpEff.vcf.gz.tbi"
    File snpeff_summary = "stats/~{vcf_name}.snpEff_summary.html"
  }
  runtime {
    docker: docker
    memory: mem_gb + " G"
    cpu: 2
  }
}

task RefilterVCF {
  input {
    String filter_script
    File vcf_f_in
    File vcf_f_in_tbi
    Int min_GQ
    Float min_AD
    Int min_tot_DP
    String vcf_name
    String python_path
    String docker
    Int mem_gb
    String name_stem_out
    String java_path
    String snpsift_jar
    String snpeff_jar
    String snpeff_config
    String snpeff_data_dir
    String ref
    File transcript_list
  }
  command <<<
    set -euo pipefail
    
    # Decompress VCF (filtering py script does not work on compressed VCFs)
    bcftools view ~{vcf_f_in} > ~{vcf_name}.nonbinary.vcf
        
    # Filter VCF, remove unused alleles (since there likely will be many), update INFO
    ~{python_path} ~{filter_script} \
      --min_GQ ~{min_GQ} \
      --min_percent_alt_in_AD ~{min_AD} \
      --min_total_DP ~{min_tot_DP} \
      ~{vcf_name}.nonbinary.vcf 2> ~{name_stem_out}.variant_qc_genotype_filter.tsv | \
      bcftools view -a -O u | \
      bcftools view --min-alleles 2 -O u | \
      bcftools view -O u -e 'COUNT(GT="ref")+COUNT(GT="mis")=N_SAMPLES' | \
      bcftools +fill-tags -O z -o ~{name_stem_out}.allele_cleanup.vcf.gz -- -t AN,AC,AF,MAF,F_MISSING

    # Index VCF
    tabix -f ~{name_stem_out}.allele_cleanup.vcf.gz

    # Remove unused alleles from snpeff annotations
    # (To do this, snpeff documentation says to remove all info and then reannotate,
    # though a python script to match Alleles in ANN vs alleles in ALT would be useful)
    # First clear snpeff annotations
    ~{java_path} -Xmx~{mem_gb}G -jar ~{snpsift_jar} \
      rmInfo ~{name_stem_out}.allele_cleanup.vcf.gz ANN > ~{name_stem_out}.allele_cleanup.rmANN.vcf

    # Now reannotate
    export SNPEFF_HOME=$(dirname ~{snpeff_jar})
    mkdir stats
    ~{java_path} -Xmx~{mem_gb}G -jar ~{snpeff_jar} -i vcf -v \
      -s stats/~{name_stem_out}.allele_cleanup.snpEff_summary.html \
      -c $SNPEFF_HOME/snpEff.config \
      -onlyTr ~{transcript_list} \
      -nodownload \
      -datadir ~{snpeff_data_dir} \
      ~{ref} ~{name_stem_out}.allele_cleanup.rmANN.vcf > ~{name_stem_out}.vcf
    bgzip ~{name_stem_out}.vcf
    tabix -f ~{name_stem_out}.vcf.gz
    

    # Clean up files
    if [[ $? -ne 0 ]]; then
      echo "Filtering failed"
      exit 1
    else
      rm ~{vcf_name}.nonbinary.vcf
      rm ~{name_stem_out}.allele_cleanup.vcf.gz
      rm ~{name_stem_out}.allele_cleanup.vcf.gz.tbi
      rm ~{name_stem_out}.allele_cleanup.rmANN.vcf
    fi
  >>>
  output {
    File filtered_vcf = "~{name_stem_out}.vcf.gz"
    File filter_tsv = "~{name_stem_out}.variant_qc_genotype_filter.tsv"
  }
  runtime {
    docker: docker
    memory: mem_gb + " G"
    cpu: 2
  }
}
