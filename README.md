# SACOMA - Spatially-Aware Clustering for Co-Methylation Analysis

This repository contains all the R codes for implementing SACOMA and reproducing the findings that can be found in our research article - A spatially-aware unsupervised pipeline to identify co-methylation regions in DNA methylation data.

## Implementation

For the implementation of SACOMA, the following R packages are needed to be installed beforehand.

For CRAN packages:

``` r
install.packages(c("data.table", "foreach", "doParallel", "ClustGeo", "readr", "dplyr", "dendextend", "peakRAM"))
```

For Bioconductor packages:

``` r
install.packages("BiocManager")
BiocManager::install("bumphunter")
```

The main function used in this workflow can be found: **[SACOMA Function](https://github.com/SChatLab/SACOMA/tree/main/SACOMA%20Implementation)** 

To run SACOMA, you will also need array manifest files for the Infinium 450k or EPIC platforms.  
These can be downloaded from following links:

- **[450k manifest](https://indiana-my.sharepoint.com/:t:/g/personal/suvchat_iu_edu/EUDB4zlBM7NPt6BLLVwiq8sBuT1pdhKPcsWlXuD_-jNc2A?e=UlUx8U)**
- **[EPIC manifest](https://indiana-my.sharepoint.com/:t:/g/personal/suvchat_iu_edu/EWOBny4uHHpDhef2_0OQcisBxvtm0kHPMBX9zzNK4Y-puQ?e=d4dJGZ)**  

Download the files and point the `manifest.dir` argument to the folder where they are stored when running the function.

## Arguments

`Sacoma` has nine main arguments which the user needs to supply/define.

### The table below details the required arguments:

| Parameter   | Default  | Description                                                                                                                                                                                                                                                         |
|:------------|:--------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| bval.dnam    |          | A dataframe of beta values with CpGs in the row and samples in the column.|
| data_type    |          | The type of array the data is sourced from. The user can choose between '450k' and 'EPIC'.|
| manifest.dir |       | Path to the directory where manifest files for 450k or EPIC array are stored.|
| minCpGs      | 3 | Minimum number of CpGs in a region.|
| minCpGs.thres.type      |     | Specifies whether regions must contain exactly minCpGs CpGs ("equal") or at least that many ("gr.equal").|
| maxGap      |  200   | maximum allowed genomic distance (in base pairs) between adjacent CpG sites for them to be grouped into the same region.|
| rthresold    |  0.5   | Threshold for minimum correlation among CpGs in a region.|
| method     |    | Chooses the method to calculate correlation among CpGs in a region.|
| ncores      |    1     | The number of cores on which the analysis will be serially performed.|

## Values

`Sacoma` outputs has 8 value arguments.

### The table below details the values returned by Sacoma:

| Value   | Description                                                                                       |
|:--------|:--------------------------------------------------------------------------------------------------|
| Probe_ID   | ID of the CpG probe included in the region.|
| CHR  | Chromosome where the CpG probe is located.|
| POS      | Genomic position of the CpG probe.|
| bin     | Bin assigned during region-building.|
| prediction     | The upper confidence interval for the genes that were analysed                                    |
| Region    | Classification of the region as signal or noise. If signal, its a co-methylated region.|
| Time.min | Total runtime (in minutes) for the Sacoma pipeline.|
| peak.memory.gb | 	Peak memory usage (in gigabytes) recorded during processing.|

## Working Example

### Importing libraries

``` r
packages <- c('foreach', 'doParallel', 'data.table', 'tidyverse',
              'readr', 'ClustGeo', 'dendextend', 'dplyr', 'peakRAM')
for (i in packages){
  print(i)
  print(packageVersion(i))
  suppressWarnings(
    suppressPackageStartupMessages(library(i, character.only = TRUE)))
}
```

### Loading Example Data

Loading a sample data

``` r
load("~/SACOMA Implementation/Data/sample_data.RData") 
```

### Snapshot of data

A typical DNA methylation beta value data frame looks something
like below. Note, the CpGs should be in the rows and the samples in
columns.

``` r
beta_vals[1:5, 1:5]
```

    ##                EUB172    EUB085      EUB051      EUB126      EUB136
    ## cg00004979 0.961173615 0.9548424 0.962959693 0.946701205 0.961559923
    ## cg00005543 0.006018365 0.0149589 0.008378839 0.004339915 0.009621588
    ## cg00007032 0.031160492 0.0351314 0.025661302 0.032440554 0.021019590
    ## cg00012194 0.876145321 0.8700672 0.900874663 0.882595185 0.893456708
    ## cg00016439 0.906904795 0.9439028 0.940404015 0.934304551 0.944407266

### Identifying co-methylated regions using SACOMA

To obtain co-methylated regions using SACOMA please use the code below.

``` r
res <- Sacoma(bval.dnam = beta_vals,
              data_type = "450k",
              manifest.dir = paste(dir, "Data", sep = '/'),
              minCpGs = 3,
              minCpGs.thres.type = "gr.equal",
              maxGap = 200,
              rthresold = 0.5, 
              method = "spearman",
              ncores = 1)
```

The results table from SACOMA should look something like the following:

``` r
res[1:5,]
```

    ##      Probe_ID   CHR       POS bin prediction Region Time.min peak.memory.gb
    ##  1 cg06550165 chr10   7618453  49      noise      1    0.092      0.4654629
    ##  2 cg00510087 chr10   7618556  49      noise      1    0.092      0.4654629
    ##  3 cg20023675 chr10   7618614  49      noise      1    0.092      0.4654629
    ##  4 cg01146320 chr11 120080945 637      noise      2    0.092      0.4654629
    ##  5 cg12526834 chr11 120081089 637      noise      2    0.092      0.4654629

Filter the results table by prediction = 'signal' to obtain co-methylated regions.

``` r
res[res$prediction == 'signal',]
```

    ##        Probe_ID  CHR      POS  bin prediction Region Time.min peak.memory.gb
    ##  174 cg02134976 chr6 42110740 6282     signal     55    0.092      0.4654629
    ##  175 cg16669650 chr6 42110833 6282     signal     55    0.092      0.4654629
    ##  176 cg02382532 chr6 42110867 6282     signal     55    0.092      0.4654629


## Analysis Scripts and Data Availability

All R scripts used for the simulation study and real-data analysis in the manuscript are available in this repository at the following locations:

- **Simulation analysis scripts:** [Simulations](https://github.com/SChatLab/SACOMA/tree/main/Simulations)
- **Real data analysis scripts:** [Real Data Analysis](https://github.com/SChatLab/SACOMA/tree/main/Real%20Data%20Analysis)

The datasets used in the paper were obtained from GEO under the following accession IDs:

- **450k dataset:** GSE281199
- **EPIC dataset:** GSE169338

After downloading, place the data files in a preferred directory and update the paths inside the analysis scripts accordingly.

# Citation

Meshram, S., Fadikar, A., Arunkumar, G., & Chatterjee, S. (2025). A spatially-aware unsupervised pipeline to identify co-methylation regions in DNA methylation data. bioRxiv, 2025-11.
