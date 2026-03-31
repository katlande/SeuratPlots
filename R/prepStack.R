# Internal function 1 - extracts image limits from coordinates file
#' @noRd
getFOVlims <- function(df){
  xmin<-min(df$x)
  xmax<-max(df$x)
  ymin<-min(df$y)
  ymax<-max(df$y)
  return(list(x=c(xmin,xmax),y=c(ymin,ymax)))
}

#' @export
prepStack <- function(p, coordDF, first=F){
  FOVlims <- getFOVlims(coordDF)
  
  prep <- p+
    ggplot2::coord_cartesian(ylim = FOVlims$y, xlim = FOVlims$x, expand = T)+
    ggplot2::theme(legend.position = "none", 
                   plot.title=ggplot2::element_blank(), 
                   plot.subtitle=ggplot2::element_blank(), 
                   plot.caption=ggplot2::element_blank()) # remove features that change the relative position of plots
  
  if(first==F){
    return(patchwork::inset_element(prep, left = 0, right=1, top=1, bottom=0, align_to = "panel"))
  } else {
    return(prep)
  }
}
