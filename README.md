# AMP-PD Allele Aware RNA Analysis

This project can be divided into three phases: 
1. Running the Variant Aware Alignment Workspace on Terra 
2. Formatting the BigWig and Junction Files 
3. Visualizing the Output 

3 Studies were included in this analysis:

|   | Number of Samples Sucessfully Processed  |
|---|---|
| BioFIND | 165 |
| PDBP | 3296 |   
| PPMI | 4622 | 

Note: Phase 1 of the analysis of PPMI data was carried out locally.

## Phase 1: Running the Variant Aware Alignment Workspace on AMP-PD RNA Seq Data

### Definitions
1. Terra:
Terra is an application that allows for users to run bioinformatic analyses using cloud-based resources. The analyses can be defined in the form of a workspace. The defined workspaces can be shared with other researchers, allowing for a collaborative environment to be fostered. 

2. Docker:
Docker is an open source platform that enables users to containerize applications. In Terra, the tools and environment required for a workflow to run are packaged into a Docker image. 

### Overview
A workspace was created on Terra.bio to carry out this analysis. There are a total of 6 workflows in the workspace.

- **Step 0** - Parse File Listing - This workflow takes a string containing the paths to the FASTQ files for a sample and converts them into a file array. Use of this workspace is optional. Alternatively, paths to the FASTQ files may be entered using the Terra GUI.
- **Step 1** - STAR Alignment with WASP Tagging - This workflow runs variant-aware STAR alignment using a set of FASTQ files (R1/R2) and one vcf file. The workflow produces an aligned, unsorted BAM file with WASP tagging.
- **Step 2** - BAM Sorting and Filtering - The BAM file produced in step 1 is used as input. The file is sorted, and filtered for passing WASP tags. The new filtered BAM files are indexed.
- **Step 3a/3b** - Junction Extraction - The BAM files are converted into junction files. This can be done using one of two workflows in this workspace. Workflow 3 utilizes the Regtools “junction extract” command. Workflow 3b utilizes a bash script, developed as part of the Leafcutter tool.
- **Step4** - Creating Bigwig File - This workspace utilizes the deepTools “bamCoverage” command to generate coverage tracks (bigwig files) of each region (chr 1-22, X, Y) of the WASP filtered BAM file.

### Tools Utilized:
- **STAR Aligner** - Version 2.7.3a
- **Samtools**  - Version 1.7.0_cv4
- **Regtools** - Version 0.5.2
- **DeepTools** - Version 3.5.0

### Scripts:
- The scripts for each workflow were written in workflow description language (wdl)
- The scripts can be found in the following folder: [Link](https://github.com/RamiyaSiva/Sashimi/tree/7b615c42712411b3d41e2104d3e9a89b1af4ee95/scripts/wdl) 
  
### Input Files:
- Per Sample Input 
  - FASTQ files 
  - VCF files 

- Common Input for All Samples
  - Genome Index file
  - Blacklist BED file 

- The input files for 3486 samples were processed, of which 3461 (165 BioFIND, 3296 PDBP) were successful 

### Output Files Per Sample:
- Junction File
  - Columns: chromosome, start position, end position, period, score, strand

    <img src="https://res.cloudinary.com/dogejctwp/image/upload/v1643638108/amp/Screen_Shot_2022-01-31_at_7.33.18_PM_jxypva.png" alt="drawing" height="200"/>

- Bigwig files seperated by chromosome:
  - Bigwig files are compressed binary files. They can be converted to the human-readable wig format using the deepTools `bigWigToWig` command.
  - The files contain header lines for each section of the genome (e.g. `#bedGraph section chr10:55700-365900`)
  - Columns: chromosome, start position, end position, coverage
    
    <img src="https://res.cloudinary.com/dogejctwp/image/upload/v1643638108/amp/Screen_Shot_2022-01-31_at_7.36.37_PM_izbdjm.png" alt="drawing" height="200"/>

## Phase 2: Formatting the BigWig and Junction Files  

### Overview
The junction and bigwig files that were produced using the variant aware alignment workspace were used as inputs for this phase of the project. 

- **Step 1a** The junction files were grouped according to samples that shared a particular characteristic(Cases vs. Controls, Sex etc). This resulted in a file in which rows represented a junction position and columns represented the counts values for a sample. The number of columns corresponded with the number of samples in a particular group. 
- **Step 2a** The median value of counts for each junction position was calculated on the grouped files using `datamash`. 

- **Step 1b** Similar to the junction files, the coverage files were grouped according to samples that shared a particular characteristic(Cases vs. Controls, Sex etc). This resulted in one large file in which rows were a particular region in chromosome and columns were coverage per sample. 
- **Step 2b** The median value of coverage for each region was calculated on the grouped files using `datamash`. 

- **Step 3** The files were formatted to include category and slice information. 

  - The grouped files can be divided into 6 categories, with 22 slices or individual groups. 
    - ALL_SAMPLES
      - ALL_SAMPLES

    - GENETIC_CONDITION
      - GBA_NEG
      - GBA_POS
      - LRRK_NEG
      - LRRK_POS
      - PD_MUT_NEG
      - PD_MUT_POS
      - SNCA_NEG
      - SNCA_POS

    - SEX
      - MALE
      - FEMALE

    - STATUS
      - CASES
      - CONTROLS

    - STUDY
      - BIOFIND
      - PDBP
      - PPMI

    - TIMEPOINT
      - BLM0T1
      - SVM6T1
      - SVM12T1
      - SVM18T1
      - SVM24T1
      - SVM36T1

Formatted Table of Junction Counts: 

  <img src="https://res.cloudinary.com/dogejctwp/image/upload/v1643639754/amp/Screen_Shot_2022-01-31_at_8.05.21_PM_cdi1wi.png" alt="drawing" height="200"/>

Formatted Table of Coverage: 

  <img src="https://res.cloudinary.com/dogejctwp/image/upload/v1643639755/amp/Screen_Shot_2022-01-31_at_8.04.43_PM_gacs97.png" alt="drawing" height="200"/>


- **Step 4** The data tables were then imported into a mongoDB database. 
  - Database Name: AMP22A
  - Collections: `coverage` and `counts`

### Tools Utilized:
- deepTools bigwig to wig 
- datamash
- mongoDB 

## Phase 3: Visualizing the Output

### Overview 
- The coverage and junction files will be used as input data to create a sashimi plot. 
- This plot will be created using vega.js, an application in which you can define interactive vizualizations in JSON. These plot can be hosted on a web browser. 

#NoteToSelf- add an image of a sashimi plot, maybe use bioRender to create it 

### Script:
- The vega script for creating a sashimi plot can be found here: [Link](https://github.com/RamiyaSiva/Sashimi/tree/7b615c42712411b3d41e2104d3e9a89b1af4ee95/scripts/vega) 

### Tools Utilized:
- vega.js

