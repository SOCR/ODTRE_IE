library(data.tree)

#' Plot optimal DTR from a DTRtree function outputted object.
#' 
#' @param dtrtree an object outputted from DTRtree function
#' @param H patient history 
plot.DTR.tree <- function(dtrtree, H) {
  node.names <- paste("node", 1:dim(dtrtree)[1], sep = "")
  #Set up an vec of NA nodes and another vec of non-NA nodes
  nonnanodes <- data.frame(
    id = 1:length(which(!is.na(dtrtree[, 2]))),
    nodenumb = which(!is.na(dtrtree[, 2])),
    missing.node = 2
  )
  #plot root
  if (is.na(dtrtree[1, 2])) {
    node1 <-  Node$new(dtrtree[1, 5])
  } else{
    node1 <-
      Node$new(paste(names(H)[dtrtree[1, 2]], " \U2264 ", round(dtrtree[1, 3], 4), sep = ""))
  }
  
  #plot nodes and edges
  k = 1
  for (i in 2:dim(dtrtree)[1]) {
    if (!is.na(dtrtree[i, 2])) {
      if (nonnanodes[k, ]$missing.node > 0) {
        nonnanodes[k, ]$missing.node = nonnanodes[k, ]$missing.node - 1
        node <-
          get(node.names[nonnanodes[k, ]$nodenumb])$AddChild(paste(names(H)[dtrtree[i, 2]], " \U2264 ", round(dtrtree[i, 3], 4), sep = ""))
        if (nonnanodes[k, ]$missing.node == 1) {
          node$`type` <- "Yes"
        } else{
          node$`type` <- "No"
        }
        assign(node.names[i], node)
      } else{
        tempk = k
        while (nonnanodes[tempk, ]$missing.node == 0) {
          tempk = tempk - 1
        }
        nonnanodes[tempk, ]$missing.node = nonnanodes[tempk, ]$missing.node - 1
        node <-
          get(node.names[nonnanodes[tempk, ]$nodenumb])$AddChild(paste(names(H)[dtrtree[i, 2]], " \U2264 ", round(dtrtree[i, 3], 2), sep = ""))
        if (nonnanodes[tempk, ]$missing.node == 1) {
          node$`type` <- "Yes"
        } else{
          node$`type` <- "No"
        }
        assign(node.names[i], node)
      }
      k = k + 1
    } else {
      if (i > nonnanodes[dim(nonnanodes)[1], 2]) {
        #If our NA node is larger than the last non NA node -->we can do the minus.
        while (nonnanodes[k, ]$missing.node == 0) {
          k = k - 1
        }
        nonnanodes[k, ]$missing.node = nonnanodes[k, ]$missing.node - 1
        node <-
          get(node.names[nonnanodes[k, ]$nodenumb])$AddChild(as.numeric(dtrtree[i, 5]))
        if (nonnanodes[k, ]$missing.node == 1) {
          node$`type` <- "Yes"
        } else{
          node$`type` <- "No"
        }
        assign(node.names[i], node)
      } else{
        tempk = k
        while (nonnanodes[tempk, ]$missing.node == 0) {
          tempk = tempk - 1
        }
        nonnanodes[tempk, ]$missing.node = nonnanodes[tempk, ]$missing.node -
          1
        node <-
          get(node.names[nonnanodes[tempk, ]$nodenumb])$AddChild(as.numeric(dtrtree[i, 5]))
        if (nonnanodes[tempk, ]$missing.node == 1) {
          node$`type` <- "Yes"
        } else{
          node$`type` <- "No"
        }
        assign(node.names[i], node)
      }
    }
  }
  
  # Utility function to plot yes and no labels 
  # @param node input node for plotting labels.
  GetEdgeLabel <- function(node) {
    label = node$type
    return (label)
  }
  SetEdgeStyle(node1, fontname = 'helvetica', label = GetEdgeLabel)
  return(node1)
}
