##
## This workflow sorts and filters the BAM file from variant aware STAR alignment
## using samtools. Reads with passing WASP tags are retained. The new filtered
## bam files are indexed.
##
## Inputs:
## - Per-sample:
##   - Sample Name
##   - BAM file
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
##   - <sample_name>.waspFiltered.bam

workflow remove_wasp_tags {
  String sample_name
  File bamFile

  String samtools_docker
  Int vm_disk_size_gb
  String vm_memory
  String runtime_zones
  Int num_cpu_cores
  Int preemptible_tries
  Float timeout_hours

  call wasp_tag_filtering {
    input:
      sample_name=sample_name,
      bamFile=bamFile,
      samtools_docker=samtools_docker,
      vm_disk_size_gb=vm_disk_size_gb,
      vm_memory=vm_memory,
      runtime_zones=runtime_zones,
      num_cpu_cores=num_cpu_cores,
      preemptible_tries=preemptible_tries,
      timeout_hours=timeout_hours
   }

  output {

    File filtered_bam = wasp_tag_filtering.filtered_bam
    File bai_output_path = wasp_tag_filtering.bai_output_path
  }
}

task wasp_tag_filtering {

  String sample_name
  File bamFile

  String samtools_docker
  Int vm_disk_size_gb
  String vm_memory
  String runtime_zones
  Int num_cpu_cores
  Int preemptible_tries
  Float timeout_hours

  command<<<
    set -o errexit
    set -o nounset
    set -o pipefail
    
    samtools sort -o "${sample_name}.samtools.bam" "${bamFile}"
    
    (samtools view -H "${sample_name}.samtools.bam"; samtools view "${sample_name}.samtools.bam" | egrep -v "vW:i:2|vW:i:3|vW:i:4|vW:i:5|vW:i:6|vW:i:7") | samtools view -b -h -o "${sample_name}.waspFiltered.bam" -

    samtools index "${sample_name}.waspFiltered.bam"
    
  >>>

  runtime {
    docker: samtools_docker

    disks: "local-disk " + vm_disk_size_gb + " HDD"

    memory: vm_memory
    cpu: num_cpu_cores
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {

    File filtered_bam = "${sample_name}.waspFiltered.bam"
    File bai_output_path = "${sample_name}.waspFiltered.bam.bai"
  }
}
