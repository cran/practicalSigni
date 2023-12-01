#' Compute the quantile for exact t test density using Taraldsen's methods
#'
#'
#' @param n {number of observastions, n-1 is degrees of freedom}
#' @param rho {True unknown population correlation coefficient, default=0}
#' @param cum {cumulative probability for which quantile is needed}
#' @return r quantile of Taraldsen's density for
#' correlation coefficient.
#' @note needs function hypergeo::hypergeo(). The quantiles are
#' rounded to 3 places and computed by numerical methods.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @importFrom hypergeo hypergeo
#' @seealso See Also as \code{\link{pvTarald}},
#' @references Taraldsen, G. "The Confidence Density for Correlation"
#' Sankhya: The Indian Journal of Statistics
#' 2023, Volume 85-A, Part 1, pp. 600-616.
#'
#' @export

qTarald=function(n,rho=0,cum){
  v=n-1
  r=seq(-1,1,by=0.001)
  Trm1=(v*(v-1)*gamma(v-1))/((sqrt(2*pi)*gamma(v+0.5)))
  Trm2=(1-r^2)^((v-1)/2)
  Trm2b=((1-rho^2)^((v-2)/2))*((1-rho*r)^((1-2*v)/2))
  Trm3b=hypergeo(3/2,-1/2,(v+0.5),(1+r*rho)/2)
  y0=Re(Trm1*Trm2*Trm2b*Trm3b)
  p=y0/sum(y0)
  cup=cumsum(p)
  loc=max(which(cup<cum))+1
  return(r[loc])}


