Workspace Parameters: 

1. STAR Alignment with WASP Tagging: 
  Docker: "gcr.io/optimum-surface-304420/sha256:e0de6c4f957c37d788ac5117063931e14d8b3e00cf46dca5ac16c3c5d3dfdce4"
  Disk Size: 200GB
  Memory: 64GB
3. BAM Sorting and Filtering:
  Docker: "biocontainers/samtools:v1.7.0_cv4"
  Disk Size: 300GB
  CPU Cores: 2
  Memory: 3.75GB
5. Junction Extraction with Regtools: 
  Docker: "gcr.io/optimum-surface-304420/sha256:69354604f562c9571adb9d3a000a819b1776f98c9c74e86bc7820dd07ab7ca56"
  Disk Size: 100GB
  CPU Cores:2
  Memory: 5GB
7. Junction Extraction with Leafcutter: 
  Docker: "gcr.io/optimum-surface-304420/sha256:bac8b238c1d13275171877c883b2903760d4915a99cac760eb9927d2c89617cd"
  Disk Size: 100GB
  CPU Cores: 8
  Memory: 30GB
9. Creating Bigwig File:
  Docker: "gcr.io/optimum-surface-304420/sha256:c71d990c8c4959574b52f234965a6cd9e4af9fc479ce25a401af63fbf98370ca"
  Disk Size: 150GB
  CPU Cores: 4
  Memory: 16GB
