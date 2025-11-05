#### Load packages ####

packages <- c('pkgmaker', 'stringr', 'foreach', 'doParallel', 'data.table', 'comeback',
              'coMethDMR', 'tidyverse', 'readr', 'ClustGeo', 'dendextend', 'SpaCCr',
              'doSNOW', 'gtools', 'dplyr', 'ChIPpeakAnno', 'pROC', 'sesameData')
for (i in packages){
  print(i)
  print(packageVersion(i))
  suppressWarnings(
    suppressPackageStartupMessages(library(i, character.only = TRUE)))
}

#### Set Path ####

dir <- 'path_to_directory/Real Data Analysis'

#### Load required functions ####

pkgmaker::source_files(paste(dir, "Codes", "fit_functions", sep = '/'),'*.R')

#### Load data ####

manifest_name <- c('HM450.h19.manifest.txt', 'EPIC.hg38.manifest.txt')

if(datatype %in% '450k'){
  if(data == 'ROSMAP_450K'){
    load(paste(dir, 'Data/450K', paste(data, 'RData', sep = '.'), sep = '/'))
    bval_dat <- data$beta_vals
  }else{
    bval_dat <- data.frame(fread(paste(dir, "Data/450K", paste(data, 'csv', sep = '.'), sep = '/')), row.names = 1)
  }
  manifest <- data.frame(fread(paste(dir, "Data", manifest_name[1], sep = '/')))
}else{
  bval_dat <- data.frame(fread(paste(dir, "Data/EPIC", paste(data, 'csv', sep = '.'), sep = '/')), row.names = 1)
  manifest <- data.frame(fread(paste(dir, "Data", manifest_name[2], sep = '/')))
  manifest <- manifest[, c(1:3, ncol(manifest), 4:(ncol(manifest)-1))]
}

bval_dat[] <- lapply(bval_dat, function(x) as.numeric(x))


bval_dat <- na.omit(bval_dat)
bval_dat[bval_dat < 0] <- 0.00001
bval_dat[bval_dat > 1] <- 0.99999

if(nrow(manifest) != nrow(bval_dat)){
  manifest <- manifest[manifest$Probe_ID %in% rownames(bval_dat), ]
}

#### Get Regions with required number of CpGs ####

maxGap <- 200
minCpgs.regions <- get.Genoregs(manifest, 
                                maxGap = maxGap, 
                                minCpgs = minCpgs, 
                                ceiling = 'gr.equal', 
                                intergenic = FALSE)

#### Getting required input elements for models ####

cpgLocs <- do.call('rbind', minCpgs.regions)
rownames(cpgLocs) <- cpgLocs$Probe_ID
cpgs_ls <- get.cpgList(minCpgs.regions)

#### Running Models on Evaluation Data ####

if(Method %in% 'Sacoma'){
  
  #### Running Sacoma on Simulated Data ####
  
  fit <- Sacoma(bval.dnam = bval_dat,
                annotated.dnam =  cpgLocs,
                bins_valid = cpgs_ls,
                minCpGs = minCpgs, 
                rthresold = 0.5, 
                method = "spearman",
                ncores = 1)
}else if(Method %in% 'coMethDMR'){
  
  #### Running coMethDMR on Simulated Data ####
  
  fit <- trigger.coMethDMR(dnam = bval_dat, 
                           annotated.dnam = cpgLocs,
                           CpGs_ls = cpgs_ls,
                           minCpGs = minCpgs, 
                           rthresold = 0.4, 
                           method = "spearman", 
                           genome = "hg19",
                           arrayType = datatype, 
                           returnAllCpGs = TRUE, 
                           output = "dataframe",
                           betaToM = TRUE, 
                           file = NULL, 
                           ncores = 1)
}else if(Method %in% 'Aclust2'){
  
  #### Running Aclust2 on Simulated Data ####
  
  fit <- trigger.Aclust2(probe.vec = rownames(cpgLocs),
                         betas = bval_dat[rownames(bval_dat) %in% rownames(cpgLocs), ],
                         manifest = manifest,
                         minimum.cluster.size = minCpgs,
                         rthreshold = 0.25,
                         maxGap = maxGap,
                         type = "average", 
                         dist.type = "spearman",
                         missingness_max_prop = 0.2)
}

#### Save Results ####

sim.scene <- paste(Method, data, datatype, minCpgs, maxGap, sep = '_')
savenam <- paste(sim.scene, ".RData", sep = '')
save(fit, file = paste(dir, "Results/Real_Data", savenam, sep = '/'))