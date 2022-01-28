## Phase 1: Perform Variant Aware Alignment on AMP-PD RNA Seq Data

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
- The scripts can be found in the following folder: Link /insert hyperlink/ 
  
### Input Files:
- Per Sample Input 
  - FASTQ files 
  - VCF files 

- Common Input for All Samples
  - Genome Index file
  - Blacklist BED file 

- The input files for 3486 samples were processed, of which 3461 were successful 
| Study | Number of Successfully Processed Samples | 
| BioFIND | 165 | 
| PDBP | 3296 | 

### Output Files Per Sample:
- Junction File 
- Bigwig files per chromosome 


## Phase 2: Use the output of Variant Aware Alignment to create informative visualizations

### Overview
The junction and bigwig files that were produced using variant aware alignment were used as inputs for this phase of the project. 

- **Step 1a** The junction files were grouped according to samples that shared a particular characteristic(Cases vs. Controls, Sex etc). This resulted in a file in which rows were with a junction position and columns represented the counts values for a sample. The number of columns corresponded with the number of samples in a particular group. 
- **Step 2a** The median value of counts for each junction was calculated for each group. 

- **Step 1b** Similar to the junction files, the coverage files were grouped according to samples that shared a particular characteristic(Cases vs. Controls, Sex etc). This resulted in one large file in which rows were a particular region in chromosome and columns were coverage per sample. 
- **Step 2b** The median value of coverage for each region was calculated on the grouped files using `datamash`. 

### Tools Utilized:
- deepTools bigwig to wig 
- datamash 
  
#note to self - group might be confusing, change the way that line is worded 

The grouped files can be divided into 6 categories, with 22 slices or individual groups. 

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

- **Step 3** The data tables were then imported into a mongoDB database. 
  - Database Name: AMP22A
  - Collections: `coverage` and `counts`

## Phase 3: Visualization of the Data


####
Pending Topics:
- Terra
- Docker
- Vega

Scripts to Upload: 
- wdl 
- vega scripts
- scripts for formatting 

