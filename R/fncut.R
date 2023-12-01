#'  fncut auxiliary converts continuous data into two categories
#'
#'
#'  This is an internal function of the R package practicalSigni
#'  Psychologists use effect size to evaluate the practical
#'  importance of a treatment on a dependent variable using
#'  a binary [0,1] variable.  Assuming numerical data, we
#'  can always compute the median and regard values < or = the
#'  median as zero and other values as unity.
#'
#' @param x {numerical vector of data values}
#' @return x vector of zeros and ones split at the median.
#' @importFrom stats median
#' @author Prof. H. D. Vinod, Fordham University, NY
#'
#'
#' @export

fncut=function(x){
  mix=median(x)
  x[x<=mix]=0
  x[x>mix]=1
  return(x)}
