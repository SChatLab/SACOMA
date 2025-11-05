######################################################################################
#### Function to make spatially contagious and correlated regions using coMethDMR ####
######################################################################################

trigger.Aclust2 <- function(probe.vec,
                            betas,
                            manifest,
                            minimum.cluster.size,
                            rthreshold,
                            maxGap,
                            type = "average", 
                            dist.type = "spearman",
                            missingness_max_prop = 0.2){
  
  #### Get Manifest ####
  
  start.time <- Sys.time()
  mem <- peakRAM::peakRAM({
    manifest.mod <- manifest %>%
      mutate(id = Probe_ID) %>%
      {rownames(.) <- NULL; .} %>%
      column_to_rownames(var = "id")
    
    #### Run Aclust2.0 region-maker function ####
    
    aclust.regions <- find_cluster_list(probe.vec = probe.vec,
                                        betas = betas,
                                        manifest = manifest.mod,
                                        minimum.cluster.size = minCpgs,
                                        thresh.dist = rthreshold,
                                        bp.thresh.dist = maxGap - 1,
                                        max.dist = maxGap,
                                        type = type,
                                        dist.type = dist.type,
                                        missingness_max_prop = missingness_max_prop)
  })
  
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  
  #### Obtaining output ####
  
  if(length(aclust.regions$clusters.list) < 1){
    outs.nosig <- rownames(betas)
    df.nosig <- data.frame(Probe_ID = outs.nosig, prediction = 'noise', Region = 0)
    df.all <- df.nosig
    df.all$Time.min <- Time.min
    df.all$peak.memory.gb <- peak.memory.gb
    return(df.all)
  }else{
    
    outs.all <- aclust.regions$clusters.list
    df.sig <- NULL
    for(l in 1:length(outs.all)){
      df.tmp <- data.frame(Probe_ID = outs.all[[l]], Region = l)
      df.sig <- data.frame(rbind(df.sig, df.tmp))
    }
    df.sig$prediction <- 'signal'
    outs.nosig <- rownames(betas)
    outs.nosig <- outs.nosig[!outs.nosig %in% df.sig$Probe_ID]
    df.nosig <- data.frame(Probe_ID = outs.nosig, prediction = 'noise', Region = 0)
    df.all <- data.frame(rbind(df.sig, df.nosig))
    df.all$Time.min <- Time.min
    df.all$peak.memory.gb <- peak.memory.gb
    return(df.all)
  }
}