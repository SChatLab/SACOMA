#### Obtain CpG-regions based on maximum base pair threshold ####

get.Genoregs <- function(manifest, 
                         maxGap = 1000, 
                         minCpgs = 3, 
                         ceiling = c('equal', 'gr.equal'), 
                         intergenic = TRUE){
  if(intergenic == FALSE){
    manifest <- dplyr::filter(manifest, manifest$UCSC_RefGene_Group != "")
  }
  manifest <- split(manifest, manifest$seqnames)
  manifest <- lapply(manifest, function(x) x[order(x$probeTarget), ])
  manifest <- do.call('rbind', manifest)
  manifest$bin <- bumphunter::clusterMaker(manifest$seqnames, manifest$probeTarget, assumeSorted = TRUE, maxGap = maxGap)
  manifest <- split(manifest, manifest$bin)
  if(ceiling %in% 'equal'){
    minCpgs.regions <- Filter(function(x) nrow(x) == minCpgs, manifest)
  }else{
    minCpgs.regions <- Filter(function(x) nrow(x) >= minCpgs, manifest)
  }
  minCpgs.names <- lapply(minCpgs.regions, function(x) paste(unique(x$seqnames), min(x$probeTarget), max(x$probeTarget), sep = ':'))
  names(minCpgs.regions) <- minCpgs.names
  return(minCpgs.regions)
}

#### Obtain correlation of CpG-regions based on maximum base pair threshold ####

get.GenoregsCorr <- function(minCpgs.regions, 
                             bval_dat,
                             method = 'spearman',
                             minCorr = 0.5,
                             maxCorr = 0.8,
                             numcores = 8){
  cl <- parallel::makeCluster(numcores)
  doParallel::registerDoParallel(cl)
  exports <- c('method', "bval_dat", 'minCorr', 'maxCorr')
  regsCorr <- foreach(l = 1:length(minCpgs.regions), .export = exports) %dopar%{
    x <- minCpgs.regions[[l]]
    bval.tmp <- bval_dat[rownames(bval_dat) %in% x$Probe_ID, ]
    cor.tmp <- cor(t(bval.tmp), method = method)
    cor.tmp <- abs(cor.tmp[upper.tri(cor.tmp, diag = FALSE)])
    cor.tmp1 <- median(cor.tmp)
    Truth <- ifelse(cor.tmp1 >= minCorr & cor.tmp1 <= maxCorr, 'signal', 'noise')
    x$cor.stats <- cor.tmp1
    x$Truth <- Truth
    return(x)
  }
  stopCluster(cl)
  names(regsCorr) <- names(minCpgs.regions)
  regsCorr.sub <- Filter(function(x) unique(x$Truth) %in% 'signal' , regsCorr)
  return(regsCorr.sub)
}

#### Obtain CpG list file for real data analysis ####

get.cpgList <- function(minCpgs.regions){
  cpgList <- lapply(minCpgs.regions, function(x) x$Probe_ID)
  return(cpgList)
}

#### Obtain Evaluations metrics for models ####

perfMetrics <- function(Method, all.cpgs, predictions, nrep){
  all.TPs <- all.cpgs[all.cpgs$Truth == 'signal', ]$Probe_ID
  resi <- predictions
  resi$adjPval <- ifelse(predictions$prediction == 'signal', 0.01, 1)
  adjPval <- resi$adjPval
  names(adjPval) <- resi$Probe_ID
  response <- rep(FALSE, length(adjPval))
  names(response) <- resi$Probe_ID
  response[names(response) %in% all.TPs] <- TRUE
  response <- factor(response, levels = c("TRUE", "FALSE"))
  predictor <- factor(adjPval <= 0.05, levels = c("TRUE","FALSE"))
  xtab <- table(predictor, response)
  cm <- caret::confusionMatrix(xtab)
  
  Sensitivity <- format(cm$byClass['Sensitivity'], scientific = FALSE)
  Sensitivity <- as.numeric(substr(as.character(Sensitivity), 1, 4))
  Specificity <- format(cm$byClass['Specificity'], scientific = FALSE)
  Specificity <- as.numeric(substr(as.character(Specificity), 1, 4))
  FDR <- format(1 - cm$byClass['Precision'], scientific = FALSE)
  FDR <- as.numeric(substr(as.character(FDR), 1, 4))
  F1 <- format(cm$byClass['F1'], scientific = FALSE)
  F1 <- as.numeric(substr(as.character(F1), 1, 4))
  MCC <- format(round(mltools::mcc(preds = predictor, actuals = response), 2), scientific = FALSE)
  MCC <- as.numeric(substr(as.character(MCC), 1, 4))
  
  #### ROC ####
  
  response <- rep(FALSE, length(adjPval))
  names(response) <- resi$Probe_ID
  response[names(response) %in% all.TPs] <- TRUE
  response <- factor(response)
  predictor <- adjPval
  proc_out1 <- pROC::roc(response = response, predictor = predictor, direction = ">", quiet = TRUE)
  ROC <- as.numeric(proc_out1$auc)
  
  #### pAUROC ####
  
  proc_out2 <- roc(response = response, predictor = predictor,
                   direction = ">",quiet = TRUE, partial.auc = c(1, 0.85), 
                   partial.auc.correct = TRUE)
  pAUROC <- as.numeric(proc_out2$auc)
  
  
  Time.min <- round(unique(predictions$Time.min), 2)
  PeakMemory_GB = unique(fit$peak.memory.gb)
  measures <- data.frame(Method = Method, Power = Sensitivity, FPR = 1 - Specificity, 
                         FDR = FDR, MCC = MCC, F1 = F1, ROC = ROC,
                         pAUROC = pAUROC, TP = xtab[1], TN = xtab[4], 
                         FP = xtab[3], FN = xtab[2], Time.min = Time.min, 
                         PeakMemory_GB = PeakMemory_GB, nrep = nrep)
  return(measures)
}

#### Function to Benchmark Models for Null Analysis ####

permMetrics <- function(Method, all.cpgs, predictions, nrep){
  Method <- Method
  n.all.cpgs <- as.numeric(nrow(all.cpgs))
  resi <- predictions
  all.FPs <- resi[resi$prediction == 'signal', ]$Probe_ID
  nDMcpgs.p <- as.numeric(length(all.FPs))
  
  Time.min <- round(unique(predictions$Time.min), 2)
  PeakMemory_GB = unique(predictions$peak.memory.gb)
  permfit <- data.frame(Method = Method,
                        nfeatures = nrow(resi),
                        nDMcpgs.p = nDMcpgs.p,
                        pDMcpgs.p = nDMcpgs.p/n.all.cpgs,
                        Time.min = Time.min,
                        PeakMemory_GB = PeakMemory_GB,
                        nrep = nrep)
  return(permfit)
}