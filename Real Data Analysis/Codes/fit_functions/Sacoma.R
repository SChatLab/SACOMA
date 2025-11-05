###################################################################################
#### Function to make spatially contagious and correlated regions using Sacoma ####
###################################################################################

Sacoma <- function(bval.dnam,
                   annotated.dnam,
                   bins_valid,
                   minCpGs, 
                   rthresold, 
                   method = "pearson",
                   ncores = 8){
  
  #### Applying Sacoma pipeline to get regions ####
  
  start.time <- Sys.time()
  mem <- peakRAM::peakRAM({
  if(ncores > 1){
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    regions.all <- foreach(r = 1:length(bins_valid), .combine = 'rbind', .packages = packages) %dopar%{
      
      #### Get Adjacency Information ####
      
      bval.dnam.sub <- bval.dnam[rownames(bval.dnam) %in% bins_valid[[r]], ]
      tmp <- annotated.dnam[rownames(annotated.dnam) %in% bins_valid[[r]], ]
      names(tmp)[c(1, 4)] <- c('CHR', 'POS')
      cor_mat <- cor(t(bval.dnam.sub), method = method)
      adjacency_matrix <- abs(cor_mat) >= rthresold
      if(all(adjacency_matrix)){
        region.tmp <- data.frame(Probe_ID = tmp$Probe_ID, CHR = tmp$CHR, POS = tmp$POS, bin = tmp$bin)
        region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
        region.tmp$prediction <- 'signal'
      }else{
        
        #### Obtaining feature and spatial distance matrices ####
        
        diag(adjacency_matrix) <- 1
        D0 <- as.dist(1 - adjacency_matrix)
        D1 <- dist(tmp[, c('start', 'end')])
        
        #### Cross-validation to obtain optimal alpha for spatial aware clustering ####
        
        range.alpha <- seq(0, 1, by = 0.1)
        cr <- ClustGeo::choicealpha(D0, D1, range.alpha, round(nrow(tmp)/3), graph = FALSE)
        diff.df <- as.data.frame(cr$Q)
        diff.df$dif <- abs(diff.df[, 1] - diff.df[, 2])
        alpha <- rownames(diff.df[diff.df$dif == min(diff.df$dif), ])
        if(length(alpha) > 1){
          alpha <- 0.5
        }else{
          alpha <- readr::parse_number(alpha)
        }
        
        #### Spatial Aware Clustering ####
        
        tree <- ClustGeo::hclustgeo(D0, D1, alpha = alpha)
        
        #### Collecting attributes for branches and leaves ####
        
        dend <- as.dendrogram(tree)
        labels <- dend %>% labels
        heights <- dend %>% dendextend::hang.dendrogram() %>% dendextend::get_leaves_attr("height")
        clust.members <- data.frame(Probe_ID = labels, heights = heights)
        branches <- dendextend::cutree_1h.dendrogram(dend, h = 0)
        branches <- data.frame(Probe_ID = names(branches), branch = as.numeric(branches))
        clust.members <- merge(clust.members, branches, id = 'Probe_ID')
        
        #### Selecting lowest height branches ####
        
        if(nrow(clust.members[clust.members$heights == min(clust.members$heights), ]) < nrow(clust.members[clust.members$heights <= 0, ])){
          clust.members <- clust.members[clust.members$heights <= 0, ]
        }else{
          clust.members <- clust.members[clust.members$heights == min(clust.members$heights), ]
        }
        
        #### Filtering lowest height branches which are not with minimum CpGs ####
        
        if(nrow(clust.members) < minCpGs){
          region.tmp <- data.frame(Probe_ID = tmp$Probe_ID, CHR = tmp$CHR, POS = tmp$POS, bin = tmp$bin)
          region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
          region.tmp$prediction <- 'noise'
          clust.members.new <- NULL
        }else if(nrow(clust.members[clust.members$heights == min(clust.members$heights), ]) < nrow(clust.members[clust.members$heights <= 0, ])){
          clust.members.new <- clust.members
        }else{
          clust.ls <- split(clust.members, clust.members$branch)
          clust.ls <- Filter(function(x) nrow(x) >= minCpGs, clust.ls)
          if(length(clust.ls) == 0){
            region.tmp <- data.frame(Probe_ID = tmp$Probe_ID, CHR = tmp$CHR, POS = tmp$POS, bin = tmp$bin)
            region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
            region.tmp$prediction <- 'noise'
            clust.members.new <- NULL
          }else{
            clust.members.new <- do.call('rbind', clust.ls)
          }
        }
        
        #### Organizing output ####
        
        if(!is.null(clust.members.new)){
          rownames(clust.members.new) <- 1:nrow(clust.members.new)
          rownames(tmp) <- 1:nrow(tmp)
          tmp.new <- tmp[tmp$Probe_ID %in% clust.members.new$Probe_ID, c('Probe_ID', 'CHR', 'POS', 'bin')]
          tmp.new <- tmp.new %>% mutate(grp = cumsum(c(TRUE, diff(as.numeric(rownames(tmp.new)))!= 1)))
          if(length(unique(tmp.new[, 'grp'])) == 1){
            region.tmp <- merge(tmp.new, clust.members.new, id = 'Probe_ID')
            region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
            region.tmp <- region.tmp[, c('Probe_ID', 'CHR', 'POS', 'bin')]
            region.tmp$prediction <- 'signal'
          }else{
            tmp.ls <- split(tmp.new, tmp.new$grp)
            tmp.ls <- Filter(function(x) nrow(x) >= minCpGs, tmp.ls)
            if(length(tmp.ls) == 0){
              region.tmp <- data.frame(Probe_ID = tmp$Probe_ID, CHR = tmp$CHR, POS = tmp$POS, bin = tmp$bin)
              region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
              region.tmp$prediction <- 'noise'
            }else{
              tmp.new <- do.call('rbind', tmp.ls)
              region.tmp <- merge(tmp.new, clust.members.new, id = 'Probe_ID')
              region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
              region.tmp <- region.tmp[, c('Probe_ID', 'CHR', 'POS', 'bin')]
              region.tmp$prediction <- 'signal'
            }
          }
        }
      }
      return(region.tmp)
    }
    stopCluster(cl)
  }else{
    regions.all <- NULL
    for(r in 1:length(bins_valid)){
      
      #### Get Adjacency Information ####
      
      bval.dnam.sub <- bval.dnam[rownames(bval.dnam) %in% bins_valid[[r]], ]
      tmp <- annotated.dnam[rownames(annotated.dnam) %in% bins_valid[[r]], ]
      names(tmp)[c(1, 4)] <- c('CHR', 'POS')
      cor_mat <- cor(t(bval.dnam.sub), method = method)
      adjacency_matrix <- abs(cor_mat) >= rthresold
      if(all(adjacency_matrix)){
        region.tmp <- data.frame(Probe_ID = tmp$Probe_ID, CHR = tmp$CHR, POS = tmp$POS, bin = tmp$bin)
        region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
        region.tmp$prediction <- 'signal'
      }else{
        
        #### Obtaining feature and spatial distance matrices ####
        
        diag(adjacency_matrix) <- 1
        D0 <- as.dist(1 - adjacency_matrix)
        D1 <- dist(tmp[, c('start', 'end')])
        
        #### Cross-validation to obtain optimal alpha for spatial aware clustering ####
        
        range.alpha <- seq(0, 1, by = 0.1)
        cr <- ClustGeo::choicealpha(D0, D1, range.alpha, round(nrow(tmp)/3), graph = FALSE)
        diff.df <- as.data.frame(cr$Q)
        diff.df$dif <- abs(diff.df[, 1] - diff.df[, 2])
        alpha <- rownames(diff.df[diff.df$dif == min(diff.df$dif), ])
        if(length(alpha) > 1){
          alpha <- 0.5
        }else{
          alpha <- readr::parse_number(alpha)
        }
        
        #### Spatial Aware Clustering ####
        
        tree <- ClustGeo::hclustgeo(D0, D1, alpha = alpha)
        
        #### Collecting attributes for branches and leaves ####
        
        dend <- as.dendrogram(tree)
        labels <- dend %>% labels
        heights <- dend %>% dendextend::hang.dendrogram() %>% dendextend::get_leaves_attr("height")
        clust.members <- data.frame(Probe_ID = labels, heights = heights)
        branches <- dendextend::cutree_1h.dendrogram(dend, h = 0)
        branches <- data.frame(Probe_ID = names(branches), branch = as.numeric(branches))
        clust.members <- merge(clust.members, branches, id = 'Probe_ID')
        
        #### Selecting lowest height branches ####
        
        if(nrow(clust.members[clust.members$heights == min(clust.members$heights), ]) < nrow(clust.members[clust.members$heights <= 0, ])){
          clust.members <- clust.members[clust.members$heights <= 0, ]
        }else{
          clust.members <- clust.members[clust.members$heights == min(clust.members$heights), ]
        }
        
        #### Filtering lowest height branches which are not with minimum CpGs ####
        
        if(nrow(clust.members) < minCpGs){
          region.tmp <- data.frame(Probe_ID = tmp$Probe_ID, CHR = tmp$CHR, POS = tmp$POS, bin = tmp$bin)
          region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
          region.tmp$prediction <- 'noise'
          clust.members.new <- NULL
        }else if(nrow(clust.members[clust.members$heights == min(clust.members$heights), ]) < nrow(clust.members[clust.members$heights <= 0, ])){
          clust.members.new <- clust.members
        }else{
          clust.ls <- split(clust.members, clust.members$branch)
          clust.ls <- Filter(function(x) nrow(x) >= minCpGs, clust.ls)
          if(length(clust.ls) == 0){
            region.tmp <- data.frame(Probe_ID = tmp$Probe_ID, CHR = tmp$CHR, POS = tmp$POS, bin = tmp$bin)
            region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
            region.tmp$prediction <- 'noise'
            clust.members.new <- NULL
          }else{
            clust.members.new <- do.call('rbind', clust.ls)
          }
        }
        
        #### Organizing output ####
        
        if(!is.null(clust.members.new)){
          rownames(clust.members.new) <- 1:nrow(clust.members.new)
          rownames(tmp) <- 1:nrow(tmp)
          tmp.new <- tmp[tmp$Probe_ID %in% clust.members.new$Probe_ID, c('Probe_ID', 'CHR', 'POS', 'bin')]
          tmp.new <- tmp.new %>% mutate(grp = cumsum(c(TRUE, diff(as.numeric(rownames(tmp.new)))!= 1)))
          if(length(unique(tmp.new[, 'grp'])) == 1){
            region.tmp <- merge(tmp.new, clust.members.new, id = 'Probe_ID')
            region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
            region.tmp <- region.tmp[, c('Probe_ID', 'CHR', 'POS', 'bin')]
            region.tmp$prediction <- 'signal'
          }else{
            tmp.ls <- split(tmp.new, tmp.new$grp)
            tmp.ls <- Filter(function(x) nrow(x) >= minCpGs, tmp.ls)
            if(length(tmp.ls) == 0){
              region.tmp <- data.frame(Probe_ID = tmp$Probe_ID, CHR = tmp$CHR, POS = tmp$POS, bin = tmp$bin)
              region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
              region.tmp$prediction <- 'noise'
            }else{
              tmp.new <- do.call('rbind', tmp.ls)
              region.tmp <- merge(tmp.new, clust.members.new, id = 'Probe_ID')
              region.tmp <- region.tmp[order(region.tmp$CHR, region.tmp$POS), ]
              region.tmp <- region.tmp[, c('Probe_ID', 'CHR', 'POS', 'bin')]
              region.tmp$prediction <- 'signal'
            }
          }
        }
      }
      regions.all <- data.frame(rbind(regions.all, region.tmp))
    }
  }
  
  #### Finalizing output ####
  })
  peak.memory.gb <- mem$Peak_RAM_Used_MiB * (2^20) / (10^9)
  stop.time <- Sys.time()
  Time.min <- round(difftime(stop.time, start.time, units = 'mins')[[1]], 3)
  reps <- regions.all %>% dplyr::select(bin, Probe_ID) %>% unique() %>% dplyr::group_by(bin) %>% dplyr::summarize(n = n())
  regions.all$Region <- rep(1:length(unique(regions.all[, 'bin'])), times = reps$n)
  regions.all$Time.min <- Time.min
  regions.all$peak.memory.gb <- peak.memory.gb
  return(regions.all)
}