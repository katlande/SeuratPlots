#' @export
findNeighbors <- function(seuratObj, seed.col=NULL, nNeigh=8, split.by=NULL, 
                          xname=NULL, yname=NULL, neighborCol="isNeighbor"){
  
  meta <- seuratObj@meta.data
  
  if(any(colnames(meta) == neighborCol)){
    message(paste0("Column named '",neighborCol,"' already exists! Overwritting this column..."))
    meta <- meta[-c(which(colnames(meta) == neighborCol))]
  }
  
  
  meta$isNeighbourColNew <- FALSE
  colnames(meta)[which(colnames(meta) == seed.col)] <- "SEEDCOLUMN"
  
  # change x and y columns in meta data to user-specified value
  if(! is.null(xname)){ colnames(meta)[which(colnames(meta) == xname)] <- "x" }
  if(! is.null(yname)){ colnames(meta)[which(colnames(meta) == yname)] <- "y" }
  
  # if multiple FOVs/Samples/etc are present in the object, check each separately:
  if(! is.null(split.by)){
    colnames(meta)[which(colnames(meta) == split.by)] <- "SeedSample"
  } else {
    meta$SeedSample <- "Running cells as one area..."
  }
  
  for(i in unique(meta$SeedSample)){
    message(i)
    oneSamp <- subset(meta, SeedSample == i)
    hits <- subset(oneSamp, SEEDCOLUMN == TRUE)
    #tmp <- subset(oneSamp, SEEDCOLUMN == FALSE)
    
  }
  
  if(nrow(hits) >0){
    for(i in 1:nrow(hits)){
      
      rn_hits <- row.names(oneSamp)[as.numeric(RANN::nn2(query=hits[c("x", "y")][i,], 
                                                         data=oneSamp[c("x", "y")], 
                                                         k=(nNeigh+1))[[1]])]
      rn <- rn_hits[!rn_hits %in% row.names(hits)[i]]
      
      meta$isNeighbourColNew[row.names(meta)%in%rn] <- TRUE
    }
  } else {
    warning(paste("No cells were found to be 'TRUE' in seed column:", seed.col))
  }
  
  
  
  
  
  # change column names back:
  colnames(meta)[which(colnames(meta) == "SEEDCOLUMN")] <- seed.col
  colnames(meta)[which(colnames(meta) == "isNeighbourColNew")] <- neighborCol
  
  if(is.null(split.by)){
    meta <- meta[-c(which(colnames(meta) ==  "SeedSample"))]
  } else {
    colnames(meta)[which(colnames(meta) ==  "SeedSample")] <- split.by
  }
  
  if(! is.null(xname)){ colnames(meta)[which(colnames(meta) == "x")] <- xname }
  if(! is.null(yname)){ colnames(meta)[which(colnames(meta) == "y")] <- yname }
  
  seuratObj@meta.data <- meta
  return(seuratObj)
}

