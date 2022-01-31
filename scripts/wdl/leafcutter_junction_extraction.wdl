##
## This workflow utilizes a bash script, developed as part of the leafcutter tool implementation
## to extract the exon-exon junctions from the WASP filtered BAM files.
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
##   - <sample_name>.leafcutter.junc.gz
##


 workflow leafcutter_bam2junc {
    File wasp_filtered_bam
    String sample_name
    
    String leafcutter_docker
    Int vm_disk_size_gb
    String vm_memory
    String runtime_zones
    Int num_cpu_cores
    Int preemptible_tries
    Float timeout_hours


	call leafcutter_bam_to_junc {
    	input:
    		wasp_filtered_bam=wasp_filtered_bam,
    		sample_name=sample_name,
    		leafcutter_docker=leafcutter_docker,
  			vm_disk_size_gb=vm_disk_size_gb,
  			vm_memory=vm_memory,
  			runtime_zones=runtime_zones,
  			num_cpu_cores=num_cpu_cores,
  			preemptible_tries=preemptible_tries,
            timeout_hours=timeout_hours
    	}
        
    output {
    	File leafcutter_junction_file= leafcutter_bam_to_junc.leafcutter_junction_file
       }
     }

task leafcutter_bam_to_junc {

    File wasp_filtered_bam
    String sample_name
    
    String leafcutter_docker
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

		echo "$(date "+%Y-%m-%d %H:%M:%S") Starting Script"


		bed_file=${wasp_filtered_bam}.bed
        leafcutter_junction_file=${sample_name}.leafcutter.junc

		echo "$(date "+%Y-%m-%d %H:%M:%S") Bam to Junc"
		timeout "${timeout_hours}"h samtools view ${wasp_filtered_bam} | /usr/bin/filter_cs.py | /usr/bin/sam2bed.pl --use-RNA-strand - $bed_file
        timeout "${timeout_hours}"h /usr/bin/bed2junc.pl $bed_file $leafcutter_junction_file

		echo "$(date "+%Y-%m-%d %H:%M:%S") Zipping Files"
        gzip $leafcutter_junction_file

		echo "$(date "+%Y-%m-%d %H:%M:%S") Done"
    >>>

    runtime {
        docker: leafcutter_docker
        memory: vm_memory
        disks: "local-disk " + vm_disk_size_gb + " HDD"
        cpu: num_cpu_cores
        preemptible: preemptible_tries
       }

    output {
        File leafcutter_junction_file="${sample_name}.leafcutter.junc.gz"
    }
  }
