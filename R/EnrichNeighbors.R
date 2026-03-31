#' @export
EnrichNeighbors <- function(seuratObj, neighborCol, group.by, exclude.seeds=F, seedCol=NULL, split.by=NULL){
  meta <- seuratObj@meta.data
  
  if(exclude.seeds==T){
    if(is.null(seedCol)){
      warning("Cannot exclude seeds when seedCol is set to NULL; proceeding with seeds included.")
    } else {
      s <- nrow(meta)
      meta <- meta[!meta[[seedCol]]==T,]
      message(paste0("Removed ", (s-nrow(meta)), "/", s, " rows containing seeds."))
    }
  }
  
  # make a new query column:
  if(!is.null(split.by)){
    meta$querycol <- paste0(gsub("\\-", "\\_", meta[[group.by]]), "-",  gsub("\\-", "\\_", meta[[split.by]]))
  } else {
    meta$querycol <- meta[[group.by]]
  }
  
  meta.neighb <- meta[meta[[neighborCol]]==T,]
  meta.bkgd <- meta[meta[[neighborCol]]==F,]
  
  df_out <- data.frame(Group=character(),
                       enrichRatio=numeric(),
                       nCells=numeric(),
                       pVal=numeric())
  
  # check each group for enrichment:
  message(paste0("Making ", length(unique(meta$querycol)), " enrichment queries..."))
  for(i in unique(meta$querycol)){
    
    # fisher test:
    #                 neighbor    background
    # in group          a             c
    # not in group      b             d
    a <- nrow(meta.neighb[meta.neighb$querycol == i,])
    b <- nrow(meta.neighb[! meta.neighb$querycol == i,])
    c <- nrow(meta.bkgd[meta.bkgd$querycol == i,])
    d <- nrow(meta.bkgd[! meta.bkgd$querycol == i,])
    p <- stats::fisher.test(matrix(rbind(a,b,c,d), nrow=2, ncol=2))$p.value
    e <- (a/b)/(c/d)
    df_out <- rbind(df_out, data.frame(Group=i, enrichRatio=e, nCells=sum(c(a,c)), pVal=p))
    
  }
  
  # adjust p-value if there are more than 3 comparisons:
  if(nrow(df_out) > 3){
    df_out$pAdj <- stats::p.adjust(df_out$pVal)
  }
  
  if(!is.null(split.by)){
    df_out <- tidyr::separate(df_out, Group, into=c("Group", "Split"), sep="\\-")
  }
  
  return(df_out)
}
