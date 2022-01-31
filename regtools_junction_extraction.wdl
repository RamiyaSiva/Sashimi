##
## This workflow utilizes the regtools "junctions extract" command to
## extract the exon-exon junctions from the WASP filtered BAM files.
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
##   - <sample_name>.junc
##
 workflow junction_extract_regtools {
    File filtered_bam
    String sample_name
    String regtools_docker
    Int vm_disk_size_gb
    String vm_memory
    String runtime_zones
    Int num_cpu_cores
    Int preemptible_tries
    Float timeout_hours
    File index_file



	call junction_extract {
    	input:
    		filtered_bam=filtered_bam,
    		sample_name=sample_name,
    		regtools_docker=regtools_docker,
  			vm_disk_size_gb=vm_disk_size_gb,
  			vm_memory=vm_memory,
  			runtime_zones=runtime_zones,
  			num_cpu_cores=num_cpu_cores,
  			preemptible_tries=preemptible_tries,
            timeout_hours=timeout_hours,
            index_file=index_file
    	}

    output {
    	File regtools_junction_file= junction_extract.regtools_junction_file
       }
     }

task junction_extract {
	File index_file
    File filtered_bam
    String sample_name

    String regtools_docker
    Int vm_disk_size_gb
    String vm_memory
    String runtime_zones
    Int num_cpu_cores
    Int preemptible_tries
    Float timeout_hours


    command <<<
        set -o errexit
    	set -o nounset
    	set -o pipefail

		echo "$(date "+%Y-%m-%d %H:%M:%S") Copying Index File to WD"
        
        cp ${index_file} ./
        
        echo "$(date "+%Y-%m-%d %H:%M:%S") Extracting Junction File"
        
        regtools junctions extract -a 8 -m 50 -M 500000 -s 1 ${filtered_bam} > ${sample_name}.junc
		
		echo "$(date "+%Y-%m-%d %H:%M:%S") Done"
    >>>

    runtime {
        docker: regtools_docker
        memory: vm_memory
        disks: "local-disk " + vm_disk_size_gb + " HDD"
        cpu: num_cpu_cores
        preemptible: preemptible_tries
        checkpointFile: "my_checkpoint"
       }

    output {
        File regtools_junction_file= "${sample_name}.junc"
        Array[String] out = read_lines("my_checkpoint")
    }
  }
