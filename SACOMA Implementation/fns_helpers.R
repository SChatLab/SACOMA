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

#### Obtain CpG list ####

get.cpgList <- function(minCpgs.regions){
  cpgList <- lapply(minCpgs.regions, function(x) x$Probe_ID)
  return(cpgList)
}