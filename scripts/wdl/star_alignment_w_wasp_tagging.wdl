## This workflow runs variant aware STAR alignment on a pair of FASTQ file lists and one vcf file. The workflow produces an aligned, unsorted BAM with WASP tagging.
##
## Inputs:
## - Per-sample:
##   - Sample name
##   - R1 FASTQ file list
##   - R2 FASTQ file list
##   - VCF file
##
## - Reference:
##   - STAR index directory tar file (.tar.gz)
##
## - VM configuration
##   - docker image url
##   - VM disk size
##   - VM memory
##   - Runtime zones, ex: "us-central1-c us-central1-b"
##   - Number of times to try the workflow with a preemptible VM before
##     falling back to a full-price VM.
##
## Outputs:
## - Per-sample
##   - <sample-id>.Aligned.out.bam
##   - <sample-id>.Log.final.out
##   - <sample-id>.Log.out
##   - <sample-id>.Log.progress.out

workflow RNAAlignment {
  String sample_name
  Array[String] fastq_1
  Array[String] fastq_2
  File vcf

  # WDL does not support localizing a directory.
  # Workaround that by passing a TAR file.
  File star_index

  String star_docker
  Int star_vm_disk_size_gb
  String star_vm_memory
  String runtime_zones
  Int preemptible_tries
  Float star_timeout_hours

  call star_align {
    input:
      fastq_1 = fastq_1,
      fastq_2 = fastq_2,
      vcf = vcf,
      sample_name = sample_name,
      star_index = star_index,
      star_vm_disk_size_gb = star_vm_disk_size_gb,
      star_vm_memory = star_vm_memory,
      star_docker = star_docker,
      runtime_zones = runtime_zones,
      preemptible_tries = preemptible_tries,
      star_timeout_hours = star_timeout_hours,
   }

  output {
    Array[File] outputFiles = star_align.outputFiles
    File bamFile = star_align.bamFile

  }
}

task star_align {

  Array[File] fastq_1
  Array[File] fastq_2
  File vcf
  String sample_name
  File star_index

  Int star_vm_disk_size_gb
  String star_vm_memory
  String star_docker
  String runtime_zones
  Int preemptible_tries
  Float star_timeout_hours

  command<<<
    set -o errexit
    set -o nounset
    set -o pipefail
    set -o xtrace

    # Untar the STAR index
    mkdir -v star_index
    tar vxfz "${star_index}" --directory star_index --strip-components=2	
    
    mkdir output_dir
	
    # Run STAR --runMode alignReads
    # We run using timeout to prevent long-running workflows
    # from wasteful preemption loops.
    # The requirements for variant aware alignment is provided as an argument 
    timeout "${star_timeout_hours}"h STAR \
      --genomeDir star_index \
      --runMode alignReads \
      --readFilesCommand zcat \
      --readFilesIn ${sep=',' fastq_1} ${sep=',' fastq_2} \
      --runThreadN 16 \
      --varVCFfile ${vcf} \
      --waspOutputMode SAMtag \
      --outSAMtype BAM Unsorted \
      --outSAMattributes NH HI AS nM NM MD vA vG vW \
      --twopassMode Basic \
      --outFileNamePrefix output_dir/${sample_name}. 



	#Printing what is in the output_dir 
    # Move the resulting BAM to its own directory
     	mkdir -v output_dir/BAM
	    cp -f output_dir/*.bam output_dir/BAM
        ls output_dir/BAM
        rm output_dir/*.bam
  >>>

  runtime {
    docker: star_docker

    disks: "local-disk " + star_vm_disk_size_gb + " HDD"

    memory: star_vm_memory
    cpu: 16
    preemptible: preemptible_tries
    zones: runtime_zones
    checkpointFile: "my_checkpoint"
  }
  output {
    File bamFile = "output_dir/BAM/${sample_name}.Aligned.out.bam"
    Array[File] outputFiles = glob("output_dir/*")
    Array[String] out = read_lines("my_checkpoint")
  }
}
