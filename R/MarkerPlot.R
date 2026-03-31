# Internal function 1 - pulls necessary data from seurat object
#' @noRd
PullDotplotData <- function (object, features, assay = NULL, 
                             cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5,
                             dot.min = 0, dot.scale = 6, 
                             idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                             scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) {
  assay <- assay %||% Seurat::DefaultAssay(object = object)
  Seurat::DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = RColorBrewer::brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = ggplot2::scale_size, 
                       radius = ggplot2::scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = Seurat::CellsByIdentities(object = object, cells = colnames(object[[assay]]), 
                                        idents = idents))
  data.features <- Seurat::FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Seurat::Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- Seurat::FetchData(object = object, vars = split.by)[cells, 
                                                          split.by]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop(paste0("Need to specify at least ", length(x = unique(x = splits)), 
                    " colors using the cols parameter"))
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, 
                              sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat::PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[stats::hclust(d = stats::dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = log1p(data.use))
                               data.use <- Seurat::MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  color.by <- ifelse(test = split.colors, yes = "colors", 
                     no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  
  return(data.plot)
}

# Internal function 2 - clusters columns in plot
#' @noRd
cluster_mat <- function(mat, distance="correlation", method="ward.D2"){
  

  
  if(!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
    stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
  }
  
  
  if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & class(distance) != "dist"){
    stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
  }
  
  
  if(distance[1] == "correlation"){
    d = stats::as.dist(1 - stats::cor(t(mat)))
  }
  
  
  else{
    if(class(distance) == "dist"){
      d = distance
    }
    
    
    else{
      d = stats::dist(mat, method = distance)
    }
  }
  
  return(stats::hclust(d, method = method))
}


# sc marker plot function
# Annotations/genes - will be ordered alphabetically (unavoidable)
# Idents - ordered based on factor levels of input object Ident
#' @importFrom magrittr %>%
#' @export
MarkerPlot <- function(obj, genes, margin_factor=0.5, maxsize=4, label.fontsize=3, 
                       assay="RNA", show.annotations=T, cluster=T, bump_annot=0.5){
  
  genes <- stats::setNames(genes[c(1,2)], c("Gene", "Details"))
  genes <- genes[genes$Gene %in% row.names(obj@assays[[assay]]),]
  tmp <- genes
  tmp$count <- 1
  
  
  tmp %>%
    dplyr::group_by(Details) %>%
    dplyr::summarise(index=sum(count)) -> tmp
  tmp <- tmp[order(tmp$Details),]
  tmp$cumsum <- cumsum(tmp$index)
  tmp$annot_y <- tmp$cumsum - (0.5*(tmp$index))
  tmp$annot_y <- tmp$annot_y + bump_annot
  tmp$xpos <- length(unique(Seurat::Idents(obj)))+1
  
  intersects <- cumsum(as.numeric(table(genes$Details)))+0.5
  intersects <- intersects[1:(length(intersects)-1)]
  
  Seurat::DefaultAssay(obj) <- assay
  if(cluster == T){
    message("Grouping clusters by correlation...")
    test <- PullDotplotData(obj, genes$Gene)
    testmat <- reshape2::dcast(test[c(3:5)], id ~ features.plot) %>%
      tibble::column_to_rownames("id")
    
    tree_row <- cluster_mat(testmat)
    row_order <- row.names(testmat)[tree_row$order]
    Seurat::Idents(obj) <- factor(Seurat::Idents(obj), levels=row_order)
    
  }
  
  Seurat::DotPlot(obj, features = factor(genes$Gene[order(genes$Details)], 
                                         levels=genes$Gene[order(genes$Details)]))+
    ggplot2::scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(9, "RdBu")))+
    ggplot2::scale_size_continuous(range = c(0,maxsize))+
    ggplot2::theme_linedraw()+
    ggplot2::geom_vline(xintercept = intersects)+
    ggplot2::xlab("")+ 
    ggplot2::ylab(" ")+
    ggplot2::coord_flip(clip = "off", ylim = c(1, length(unique(Seurat::Idents(obj)))))+
    ggplot2::theme(legend.position = "left", 
                   legend.title = ggplot2::element_text(size=6),
                   axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1),
                   plot.title=ggplot2::element_text(hjust=0.5, face="bold"),
                   plot.margin=ggplot2::margin(t=1*margin_factor,
                                               r=5*margin_factor,
                                               b=1*margin_factor,
                                               l=1*margin_factor, unit = "cm")) -> d
  
  if(show.annotations==T){
    d+ggplot2::geom_text(data=tmp, 
                         mapping=ggplot2::aes(y=xpos, x=annot_y, label=Details), 
                         size=label.fontsize, hjust=0, vjust=0)->d
  }
  
  return(d)
}

# obj - seurat object
# genes - two column df with genes (col 1) and annotation label (col 2)
# margin factor - increase or decrease to scale the margins if the annotations are clipping
# maxsize - increase or decrease maximum dot size for aesthetics 
# label.fontsize - fontsize of annotation labels
