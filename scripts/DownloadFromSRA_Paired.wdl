version 1.0

# Copyright Broad Institute, 2019
#
# About:
#   This WDL pipeline downloads data from SRA in parallel and stores the results in the
#   specified GCS dir.  This pipeline is essentially a Cromwell/GCP reimagining of the
#   Nextflow/AWS downloading pipeline from @alaincoletta (see: http://broad.io/aws_dl).
#
# Description of inputs:
#   Required:
#       Array[String] SRA_IDs                 - The SRA ID(s) of the datasets to download.
#
# Licensing:
#   This script is released under the WDL source code license (BSD-3) (see LICENSE in
#   https://github.com/broadinstitute/wdl). Note however that the programs it calls may
#   be subject to different licenses. Users are responsible for checking that they are
#   authorized to run all programs before running this script.

workflow DownloadFromSRA {
    input {
        Array[String] SRA_IDs
        String NCBI_KEY="" # For security concerns, personal API key has been omitted
        Int max_disk_space_gb = 8
    }

    scatter (SRA_ID in SRA_IDs) {
        call GetSRAIDs {
            input:
                SRA_ID                         = SRA_ID,
                NCBI_KEY                       = NCBI_KEY
        }

        scatter (row in read_tsv(GetSRAIDs.sra_table)) {
            call printRow {
                input: length = length(row), row = row
            }

            # ... But, only run fasterq-dump if the library layout is PAIRED (pacbio, nanopore, and 454 are SINGLE-end reads)
            # Do not dump if either single-end OR layout not recognized
            String library_layout = row[2]
            Boolean paired_end = if (library_layout == "PAIRED") then true else false
            if (paired_end) {
              call FastqDump {
                  input:
                      SRR_ID                     = row[0],
                      uncompressed_disk_space_mb = row[1],
                      max_disk_space_gb = max_disk_space_gb
              }
            }
            if (!paired_end) {
              call SkipFastqDump {
                  input:
                      SRR_ID     = row[0],
                      lib_layout = library_layout
              }
            } 
        }
        Array[File] reads_fwd_no_none = select_all(FastqDump.reads_fwd)
        Array[File] reads_rev_no_none = select_all(FastqDump.reads_rev)
        Array[File] upload_list_no_none = select_all(FastqDump.upload_list)
        Array[File] succeed_map_no_none = select_all(FastqDump.succeed_map)
        Array[File] fail_list_no_none = select_all(SkipFastqDump.fail_list)
    }
    # DO not actually need to do a glob-based handoff. New flattening strategy works great
    output {
        Array[Boolean] is_paired = flatten(select_all(paired_end))
        Array[File] SRR_table = GetSRAIDs.sra_table
        Array[File] fail_list = flatten(fail_list_no_none)
        Array[Array[File?]] reads_fwd_nulls = FastqDump.reads_fwd
        Array[File] reads_fwd = flatten(reads_fwd_no_none)
        Array[File] reads_rev = flatten(reads_rev_no_none)
        Array[File] upload_list = flatten(upload_list_no_none)
        Array[File] succeed_map = flatten(succeed_map_no_none)
    }
}

# This task takes an SRA ID of some sort (e.g. study ID, run ID, experiment ID, etc.), finds the SRR/ERR run IDs
# associated with the provided SRA ID, and writes their accession numbers and file sizes to a tsv table.
#
# In inspecting the filesizes returned (column 8), it appears empirically to be around 5 times larger than the
# gzipped fastq files themselves.  Thus I assume this number is the uncompressed size of the file, but I cannot
# find any documentation that explicitly states this.  I've assumed a lesser compression level of 2-fold below.

task GetSRAIDs {
    input {
        String SRA_ID
        String NCBI_KEY
    }

    command <<<
        export NCBI_API_KEY=~{NCBI_KEY}
        e_fails=0
        max_e_retries=10
        run_esearch() {
            ESEARCH="esearch -db sra -query ~{SRA_ID} | efetch --format runinfo | grep '[SED]RR' | cut -f 1,8,16 -d',' | sed 's/,/\t/g' > sra.txt"
            echo $ESEARCH
            eval $ESEARCH

            # Even with API key, we are limited to 10 queries/sec. This limit can be hit with large populations so retry if output is empty
            if [ -s sra.txt ]; then
                echo "sra.txt not empty, implying Eutils completed successfully"
            else
                e_fails=${1}
                let e_fails++
                if [[ $e_fails -gt $max_e_retries ]]; then echo "Last attempt hit max tries: ${e_fails}/${max_e_retries}" && exit 2; fi
                sleep_secs=$(echo $((RANDOM % 60)))
                printf "Eutils failed on attempt:${e_fails}. Sleeping ${sleep_secs} and re-running.\n"
                sleep "${sleep_secs}"
                RERUN="run_esearch ${e_fails}"
                echo $RERUN && eval $RERUN
            fi
        }
        
        RUN="run_esearch ${e_fails}"
        echo $RUN && eval $RUN
    >>>

    output {
        File sra_table = "sra.txt"
    }

    runtime {
        maxRetries: 1
        docker: "biocontainers/ncbi-entrez-direct:v10.9.20190219ds-1b10-deb_cv1"
    }
}

# This task checks to see if a to-be-downloaded file exists at the specified GCS filepath, and if not, initiates a
# parallel download process.  The downloaded file is then immediately uploaded to the specified GCS filepath.  It
# is then deleted from local storage here by Cromwell, preventing redundant storage of the data.
task FastqDump {
    input {
        String SRR_ID
        Int uncompressed_disk_space_mb
        Int max_disk_space_gb
    }

    Int mb_to_gb_divisor = 1000
    Int compression_factor = 2
    Int disk_space_divisor = mb_to_gb_divisor * compression_factor
    Int minimum_disk_size_gb = 10
    Int compressed_disk_space_gb = ceil(uncompressed_disk_space_mb / disk_space_divisor)

    # We're going to request a LOCAL disk which comes in minimum increments of 375 GB and ignores our disk space
    # request.  But we're going to compute it anyway should we decide to change our request to an HDD/SSD later.
    Int final_disk_space_gb = if (compressed_disk_space_gb < minimum_disk_size_gb) then minimum_disk_size_gb else compressed_disk_space_gb
    
    # Fasterq-dump should be quicker than fastq-dump,
    # and split-3 resolves issues where the number of reads mismatched among the forward and reverse libraries
    command <<<
          fasterq-dump --split-3 ~{SRR_ID}
          # If failed, exit workflow so cromwell doesnt think it worked
          fq_dump_rc="$?"
          if [ $fq_dump_rc -ne 0 ]; then
              echo "Error: fasterq-dump did not succeed, exit code $fq_dump_rc"
              exit $fq_dump_rc
          else
              echo "Successful fasterq-dump, exit code $fq_dump_rc"
          fi
          rm ~{SRR_ID}.fastq
          gzip ~{SRR_ID}*
          ln -s ~{SRR_ID}* ./
          ls -lh *_1.fastq.gz > upload_list.txt
          ls -lh *_2.fastq.gz >> upload_list.txt
          
          ls *_*.fastq.gz | paste - - > ~{SRR_ID}.fastqs.tsv
          echo ~{SRR_ID} | cat - ~{SRR_ID}.fastqs.tsv | paste - - >> SRR_succeed.tsv
    >>>
  
    output {
        File? upload_list = "upload_list.txt"
        File? reads_fwd = "${SRR_ID}_1.fastq.gz"
        File? reads_rev = "${SRR_ID}_2.fastq.gz"
        File? succeed_map = "SRR_succeed.tsv"
    }

    runtime {
        maxRetries: 3
        docker: "quay.io/biocontainers/sra-tools:3.2.1--h4304569_1"
    }
}

task printRow {
    input {
        Int length
        Array[String] row
    }

    command <<<
       echo "Parsed array length" > array_parsing.txt
       echo ~{length} >> array_parsing.txt
       echo "Array values with sep=," >> array_parsing.txt
       echo ~{sep="," row} >> array_parsing.txt
    >>>

    output {
        File array_parsing_info = "array_parsing.txt"
    }
}

task SkipFastqDump {
    input {
        String SRR_ID
        String lib_layout
    }

    command <<<
       echo "Skipping fastq dump. Observed sequence library layout: ~{lib_layout}"
       echo ~{SRR_ID} >> SRR_fail.list
    >>>

    output {
       File? fail_list = "SRR_fail.list"
    }
}

