######################################################################################
#### Function to make spatially contagious and correlated regions using coMethDMR ####
######################################################################################

trigger.coMethDMR <- function(dnam, 
                              annotated.dnam,
                              betaToM = FALSE,
                              method = "pearson",
                              rthresold = 0.4, 
                              minCpGs = 3, 
                              genome = "hg19",
                              arrayType = "450k",
                              CpGs_ls,
                              file = NULL,
                              returnAllCpGs = FALSE,
                              output = "dataframe",
                              ncores = 8){
  
  #### Applying coMethDMR to get regions ####
  
  start.time <- Sys.time()
  
  mem <- peakRAM::peakRAM({
  regions <- coMethDMR::CoMethAllRegions(dnam = dnam,
                                         betaToM = betaToM,
                                         method = method,
                                         rDropThresh_num = rthresold,
                                         minCpGs = minCpGs,
                                         genome = genome,
                                         arrayType = arrayType,
                                         CpGs_ls = CpGs_ls,
                                         file = file,
                                         returnAllCpGs = returnAllCpGs,
                                         output = output,
                                         nCores_int = ncores)
  })
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  
  #### Organizing output ####
  
  regions.all <- NULL
  for(lll in 1:length(regions)){
    tmp1 <- regions[[lll]]
    if(nrow(tmp1[tmp1$keep_contiguous == 1, ]) == 0){
      tmp1$prediction <- 'noise'
    }else{
      tmp1$prediction <- 'signal'
    }
    assignedbin <- annotated.dnam[annotated.dnam$Probe_ID == tmp1$CpG[1], 'bin']
    tmp1 <- data.frame(tmp1[, c(2:5, 8)], bin = assignedbin, coMethRegion = lll, Time.min = Time.min, 
                       peak.memory.gb = peak.memory.gb)
    regions.all <- data.frame(rbind(regions.all, tmp1))
  }
  names(regions.all)[1:3] <- c('Probe_ID', 'CHR', 'Position')
  return(regions.all)
}