#### Load packages ####

packages <- c('pkgmaker', 'stringr', 'foreach', 'doParallel', 'data.table',
              'coMethDMR', 'tidyverse', 'readr', 'ClustGeo', 'dendextend',
              'doSNOW', 'gtools', 'dplyr', 'ChIPpeakAnno', 'pROC')
for (i in packages){
  print(i)
  print(packageVersion(i))
  suppressWarnings(
    suppressPackageStartupMessages(library(i, character.only = TRUE)))
}

#### Set Path ####

dir <- 'path_to_directory/Simulations/Null Simulations'

#### Load required functions ####

pkgmaker::source_files(paste(dir, "Codes", "fit_functions", sep = '/'),'*.R')

#### Load data ####

manifest_name <- c('HM450.h19.manifest.txt', 'EPIC.hg38.manifest.txt')

if(datatype %in% '450k'){
  bval_dat <- data.frame(fread(paste(dir, "Data/450K", paste(data, 'csv', sep = '.'), sep = '/')), row.names = 1)
  manifest <- data.frame(fread(paste(dir, "Data", manifest_name[1], sep = '/')))
}else{
  bval_dat <- data.frame(fread(paste(dir, "Data/EPIC", paste(data, 'csv', sep = '.'), sep = '/')), row.names = 1)
  manifest <- data.frame(fread(paste(dir, "Data", manifest_name[2], sep = '/')))
  manifest <- manifest[, c(1:3, ncol(manifest), 4:(ncol(manifest)-1))]
}

bval_dat[] <- lapply(bval_dat, function(x) as.numeric(x))


bval_dat <- na.omit(bval_dat)

if(nrow(manifest) != nrow(bval_dat)){
  manifest <- manifest[manifest$Probe_ID %in% rownames(bval_dat), ]
}
#### Get Regions with required number of CpGs ####

minCpgs.regions <- get.Genoregs(manifest, 
                                maxGap = maxGap, 
                                minCpgs = minCpgs, 
                                ceiling = 'equal', 
                                intergenic = TRUE)

#### Get Regions with required number of CpGs and required correlation range ####

numcores <- 10
qual.regions <- get.GenoregsCorr(minCpgs.regions, 
                                 bval_dat,
                                 method = 'spearman',
                                 minCorr = minCorr,
                                 maxCorr = 0.4,
                                 numcores = numcores)

#### Creating evaluation data with signals and noise ####

nreps <- 100
ncores <- 10
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
exports <- c('Sacoma', 'trigger.coMethDMR', 'trigger.Aclust2', 'perfMetrics', 
             'qual.regions', 'trigger.CoMeBack')
restmp <- foreach(r = 1:nreps, .combine = 'rbind', .packages = packages, .export = exports) %dopar% {
  set.seed(1233 + r)
  qual.regions <- qual.regions[sample(names(qual.regions), sample.out)]
  cpgLocs <- do.call('rbind', qual.regions)
  cpgLocs$Truth <- 'noise'
  cpgLocs <- split(cpgLocs, cpgLocs$seqnames)
  cpgLocs <- cpgLocs[mixedorder(names(cpgLocs))]
  cpgLocs <- lapply(cpgLocs, function(x) x[order(x$probeTarget), ])
  cpgLocs <- do.call('rbind', cpgLocs)
  rownames(cpgLocs) <- cpgLocs$Probe_ID
  cpgLocs.tmp <- split(cpgLocs, cpgLocs$bin)
  cpgs_ls <- lapply(cpgLocs.tmp, function(x) x$Probe_ID)
  
  #### Running Models on Evaluation Data ####
  
  if(Method %in% 'Sacoma'){
    
    #### Running Sacoma on Simulated Data ####
    
    fit <- Sacoma(bval.dnam = bval_dat,
                  annotated.dnam = cpgLocs,
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
  
  #### Obtaining Evaluation Metrics ####
  
  restmp <- permMetrics(Method = Method,
                        all.cpgs = cpgLocs,
                        predictions = fit,
                        nrep = r)
  
  return(restmp)
}
stopCluster(cl)

#### Save Results ####

sim.scene <- paste(Method, data, datatype, minCorr, minCpgs, maxGap, sample.out, sep = '_')
savenam <- paste(sim.scene, ".RData", sep = '')
save(restmp, file = paste(dir, "Results/Null", savenam, sep = '/'))
