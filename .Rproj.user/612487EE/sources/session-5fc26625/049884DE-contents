#' @export
getPolygons <- function(obj, fov) {
  p <- methods::slot(methods::slot(obj@images[[fov]], "boundaries")$segmentation, "polygons")
  
  poly_list <- lapply(names(obj@images[[fov]]$segmentation), function(x) {
    p2 <- p[[x]]
    coords <- methods::slot(methods::slot(p2, "Polygons")[[1]], "coords")
    data.frame(cell = x, 
               polygonX = coords[, 1], 
               polygonY = coords[, 2])
  })
  
  do.call(rbind, poly_list) -> output
  row.names(output) <- c()
  
  dplyr::left_join(output, tibble::rownames_to_column(obj@meta.data, "cell"), by = "cell") -> output
  return(output)
}
