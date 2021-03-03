#' Gene class labels
#'
#' The simulated class labels / class indicators for each gene node. The labels are signed based on simulated gene network
#'
#' @details Based on the simulated network configuration. We firstly apply the fast community detection algorithm to find all communities. And then we gradually merge two communities into one based on their smallest pair-wised betweenness defined by edge counts, untill finally we achive three large community representing three different gene classes. Finally, the largest community is assigned as null gene class, and then the upper/ down regulated are randomly chosen
#' @docType data
#' @keywords datasets
#' @name classLabels
#' @aliases classLabels
#' @usage data(classLabels)
#' @format A vector storing the class labels for each of the gene node. The vector is of length equal the total number of gene nodes. Each element can take value of 1, 2 or 3 as down-regulated gene class, null gene class (or not differentially expressed gene class) or up-regulated gene class
NULL
