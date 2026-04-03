# colour darkness checker
# modified from: https://www.w3.org/TR/AERT/#color-contrast
#' @noRd
isDark <- function(colr) { (sum( grDevices::col2rgb(colr) * c(299, 587,114))/1000 < 100) }

#' @export
CellFraction <- function(obj, group.by, split.by, label=T, label.size=3, min.perc=5, 
                         tilt=T, colors=NULL, return.data=F){
  meta <- obj@meta.data
  
  colnames(meta)[which(colnames(meta) == group.by)] <- "GROUPBYVAR"
  colnames(meta)[which(colnames(meta) == split.by)] <- "SPLITBYVAR"
  meta$NCELLCOL <- 1
  
  meta %>%
    dplyr::group_by(SPLITBYVAR) %>%
    dplyr::mutate(total_cells_in_split=sum(NCELLCOL)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(SPLITBYVAR, GROUPBYVAR, total_cells_in_split) %>%
    dplyr::summarise(total_cells_in_group=sum(NCELLCOL)) -> summary
  
  if(return.data){
    return(summary)
  } else {
    
    summary$fraction <- summary$total_cells_in_group/summary$total_cells_in_split
    
    if(is.null(colors)){
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(6, "Spectral"))(length(unique(meta$GROUPBYVAR)))
    }
    
    # check if labels should be dark or light:
    data.frame(GROUPBYVAR=levels(as.factor(summary$GROUPBYVAR)),
               fillval=colors,
               labval=ifelse(lapply(colors, isDark), "white", "black")) %>%
      merge(summary, ., by="GROUPBYVAR", all.x=T, all.y=F) -> summary
    
    ggplot2::ggplot(summary, ggplot2::aes(x=SPLITBYVAR, y=fraction, fill=GROUPBYVAR))+
      ggplot2::geom_col(colour="black")+
      ggplot2::scale_y_continuous("Fraction of Cells", expand=c(0,0))+
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
            legend.title = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour = "black", fill = NA),
            panel.background = ggplot2::element_rect(fill = "white"),
            panel.grid.major.x = ggplot2::element_line(colour = "grey", linewidth = 0.25),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(colour = "grey", linewidth = 0.25),
            panel.grid.minor.y = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(colour = "black"),
            axis.title = ggplot2::element_text(colour = "black", face = "italic"),
            axis.ticks = ggplot2::element_line(colour = "black"),
            strip.background = ggplot2::element_rect(fill="black"),
            strip.text = ggplot2::element_text(colour="white", face="bold"),
            plot.title = ggplot2::element_text(hjust = 0.5, face="bold"),
            plot.subtitle = ggplot2::element_text(hjust = 0.5),
            legend.key = ggplot2::element_rect(fill = "white"))+
      ggplot2::scale_fill_manual(values=colors) -> g
    
    if(label==T){
      
      
      g+ggplot2::geom_text(ggplot2::aes(label = ifelse(fraction*100 >= min.perc, 
                                                       paste0(formatC(fraction*100, digits = 2), "%"),"")), 
                           position = ggplot2::position_stack(vjust = .5), 
                           size=label.size,
                           color=summary$labval) -> g
    }
    
    if(tilt==T){
      g <- g+ggplot2::theme(axis.text.x = ggplot2::element_text(colour = "black", angle=45, vjust=1, hjust=1))
    }
    
    return(g)
  }
}

