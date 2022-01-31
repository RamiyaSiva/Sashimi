##
## This workflow utilizes the deepTools "bamCoverage" command to
## generate coverage tracks (bigwig files) of each region (chr1-22,X,Y)
## of the WASP filtered BAM file.
##
## Inputs:
## - Per-sample:
##   - Sample Name
##   - Filtered BAM file
##   - BAM index file
##
##
## - VM configuration
##   - Docker image url
##   - VM disk size
##   - VM memory
##   - Runtime zones, ex: "us-central1-a us-central1-b"
##   - Num CPU cores
##   - Number of times to try the workflow with a preemptible VM before
##     falling back to a full-price VM.
##
## Outputs:
##   - Array of Files: <sample_name>.<region>.bw
##
workflow create_bigwig_file {
    File filtered_bam
    String sample_name
    String deeptools_docker
    Int vm_disk_size_gb
    String vm_memory
    String runtime_zones
    Int num_cpu_cores
    Int preemptible_tries
    Float timeout_hours
    File index_file
	File black_list


	call bigwig {
    	input:
    		filtered_bam=filtered_bam,
    		sample_name=sample_name,
    		deeptools_docker=deeptools_docker,
  			vm_disk_size_gb=vm_disk_size_gb,
  			vm_memory=vm_memory,
  			runtime_zones=runtime_zones,
  			num_cpu_cores=num_cpu_cores,
  			preemptible_tries=preemptible_tries,
            timeout_hours=timeout_hours,
			index_file=index_file,
            black_list=black_list


    	}

    output {
    	Array[File] bigwig_file= bigwig.bigwig_file
       }
     }

task bigwig {
    File filtered_bam
    String sample_name

    String deeptools_docker
    Int vm_disk_size_gb
    String vm_memory
    String runtime_zones
    Int num_cpu_cores
    Int preemptible_tries
    Float timeout_hours
    File index_file
	File black_list

    command <<<
        set -o errexit
    	set -o nounset
    	set -o pipefail
        
        echo "$(date "+%Y-%m-%d %H:%M:%S") Copying Index File to WD"
        
        cp ${index_file} ./

		echo "$(date "+%Y-%m-%d %H:%M:%S")Creating BigWig File"
        
        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr1.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr1 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20
        
        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr2.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr2 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr3.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr3 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr4.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr4 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr5.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr5 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr6.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr6 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr7.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr7 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr8.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr8 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr9.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr9 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr10.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr10 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr11.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr11 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr12.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr12 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr13.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr13 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr14.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr14 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr15.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr15 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr16.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr16 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr17.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr17 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr18.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr18 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr19.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr19 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr20.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr20 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr21.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr21 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20
        
        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chr22.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chr22 --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chrX.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chrX --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

        bamCoverage -b ${filtered_bam} --outFileName ${sample_name}.chrY.bw --outFileFormat bigwig --verbose --numberOfProcessors max --region chrY --blackListFileName ${black_list} --ignoreDuplicates --minMappingQuality 20

		echo "$(date "+%Y-%m-%d %H:%M:%S") Complete."
    >>>

    runtime {
        docker: deeptools_docker
        memory: vm_memory
        disks: "local-disk " + vm_disk_size_gb + " HDD"
        cpu: num_cpu_cores
        preemptible: preemptible_tries
       }

    output {
        Array[File] bigwig_file = glob("*.bw")
    }
  }
