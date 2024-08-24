#' nmf_pull_components
#'
#' @param nmf.object NMF object outputed from either the NMF package or the CoGAPs package
#'
#' @return Returns a list of NMF outputs, including the A and P matrix, Scaled values, Maximum P values per column, colors assigned to each pattern.
#' @export
#'
#' @examples
nmf_pull_components <- function(nmf.object){

  output <- list()
  nmf <- nmf.object
  # Generate Objects:

  if (inherits(nmf, "CogapsResult")) {
    output[["A"]] <- nmf@featureLoadings
    output[["P"]] <- t(nmf@sampleFactors)
  } else {
    P <-  NMF::coef(nmf)
    A <- NMF::basis(nmf)
    rownames(P) = paste0("Pattern_", 1:nrow(P))
    colnames(A) = paste0("Pattern_", 1:ncol(A))
    output[["A"]] <- A
    output[["P"]] <- P
  }

  #Amplitude Weights
  output[["As"]] <- base::colSums(output[["A"]])

  #Normalized Pattern Matrix
  output[["Pweighted"]] <- output[["P"]]*output[["As"]]
  rownames(output[["Pweighted"]]) <- paste0("Pattern_", seq_along(1:nrow(output[["Pweighted"]])))

  #Scaled Pattern Matrix (to 100%)
  scaled =  output[["Pweighted"]] %*% diag(1/colSums( output[["Pweighted"]]))
  rownames(scaled) = paste0("Pattern_", seq_along(1:nrow(scaled)))
  colnames(scaled) = colnames( output[["Pweighted"]])
  output[["Pscaled"]] <- scaled

  #Max Pattern
  max_indices <- apply(output[["Pweighted"]], 2, which.max)
  max_indices.labels <- paste0("Pattern_",max_indices)
  names(max_indices.labels) <- names(max_indices)
  output[["PMax"]] <- max_indices

  Ps <- rowSums(output[["P"]] )
  Aweighted <- t(Ps * t(output[["A"]]))
  colnames(Aweighted) <- rownames(output[["Pweighted"]])
  output[["Aweighted"]] <- Aweighted
  output[["Ascaled"]] = Aweighted * (1/rowSums(Aweighted))
  output[["AP"]] =  output[["A"]] %*% output[["P"]]

  custom.cols <- c( "#E41A1C", "#377EB8", "#4DAF4A" , "#984EA3", "#FF7F00",  "#FFFF33",
                    "#A65628", "#F781BF","#999999", "black", "#a6cee3",
                    "#481568FF", "#74D055FF", "#C29037","#0B2D5D", "#FCC1DB" )
  #Assign Colors
  nPatterns = ncol(output[["A"]])
  if(nPatterns <= 9){
    pattern.cols <- RColorBrewer::brewer.pal(nPatterns, "Set1")
  }
  if(nPatterns >9 & nPatterns <= 16){
    pattern.cols <- custom.cols[1:nPatterns]
  }
  if(nPatterns > 16){
    cols <- colorRampPalette(custom.cols)
    pattern.cols <- cols(nPatterns)
  }
  names(pattern.cols) <- paste0("Pattern_",1:nPatterns)
  output[["pattern.cols"]] <- pattern.cols


  return(output)
}
