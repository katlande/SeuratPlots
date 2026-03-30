#' @export
PlotPolygons <- function(coordDF, fillVar, border_color="black", fillColor=NULL, fill_lims=NULL, legend="right", background=NULL, legend.name=NULL){
  
  
  if(is.null(background)){
    i <- which(colnames(coordDF) == fillVar)
    if(identical(i, integer(0))){
      warning(paste0("No column called ", fillVar, " is present in the meta.data."))
      stop()
    } else {
      colnames(coordDF)[i] <- "WORKINGFILLVAR"
    }
    
    if(is.null(legend.name)){
      legend.name <- fillVar 
    }
    
    ggplot2::ggplot(coordDF, ggplot2::aes(x = polygonX, y = polygonY, group = cell, fill=WORKINGFILLVAR)) +
      ggplot2::theme_void()+
      ggplot2::theme(plot.title=ggplot2::element_blank(), legend.position = legend)+
      ggplot2::geom_polygon(linewidth = 0.5, color=border_color) -> p
    
    if(is.null(fillColor)){
      colVec <- c("#FCBBA1", "#FB6A4A", "#CB181D", "#67000D")
    } else {
      colVec <- fillColor
    }
    
    if(is.null(fill_lims)){
      p+ ggplot2::scale_fill_gradientn(legend.name, colors=colVec) -> p
    } else{
      p+ ggplot2::scale_fill_gradientn(legend.name, colors=colVec, limits=fill_lims) -> p
    } 
  } else {
    
    ggplot2::ggplot(coordDF, ggplot2::aes(x = polygonX, y = polygonY, group = cell, fill=WORKINGFILLVAR)) +
      ggplot2::theme_void()+
      ggplot2::theme(plot.title=ggplot2::element_blank(), legend.position = legend)+
      ggplot2::geom_polygon(linewidth = 0.5, color=border_color, fill=background) -> p 
  }
  
  return(p)
  
}
