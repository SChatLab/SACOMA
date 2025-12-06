# SACOMA - Spatially-Aware Clustering for Co-Methylation Analysis

This repository contains all the R codes for implementing SACOMA and reproducing the findings that can be found in our research article - A spatially-aware pipeline to identify co-methylation regions in DNA methylation data.

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

- **450k manifest:** <link-here>  
- **EPIC manifest:** <link-here>  

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

# Citation

Meshram, S., Fadikar, A., Arunkumar, G., & Chatterjee, S. (2025). A spatially-aware unsupervised pipeline to identify co-methylation regions in DNA methylation data. bioRxiv, 2025-11.
