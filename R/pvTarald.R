#' Compute the p-value for exact correlation significance test
#' using Taraldsen's exact methods.
#'
#'
#' @param n {number of observations, n-1 is degrees of freedom}
#' @param rho {True unknown population correlation coefficient
#' in the r-interval [-1, 1], default=0}
#' @param obsr {observed r or  correlation coefficient}
#' @return ans is the p-value or probability from sampling distribution of observing a
#' correlation as extreme or more extreme than the input obsr or observed r.
#' @note needs function hypergeo from the package of that name.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @importFrom hypergeo hypergeo
#' @seealso See Also as \code{\link{qTarald}},
#' @references Taraldsen, G. "The Confidence Density for Correlation"
#' Sankhya: The Indian Journal of Statistics
#' 2023, Volume 85-A, Part 1, pp. 600-616.
#'
#' @export

pvTarald= function(n,rho=0,obsr){
  v=n-1
  r=seq(-1,1,by=0.001)
  if(v<=164)  Trm1=(v*(v-1)*gamma(v-1))/((sqrt(2*pi)*gamma(v+0.5)))
  if(v>164)  Trm1=(164*(163)*gamma(163))/((sqrt(2*pi)*gamma(163.5)))
  Trm2=(1-r^2)^((v-1)/2)
  if(rho!=0)  Trm2b=((1-rho^2)^((v-2)/2))*((1-rho*r)^((1-2*v)/2))
  if(rho==0)  Trm2b=1
  Trm3b=Re(hypergeo(3/2,-1/2,(v+0.5),(1+r*rho)/2))
  y0=Re(Trm1*Trm2*Trm2b*Trm3b)
  p=y0/sum(y0)
  cup=cumsum(p)
  loc=max(which(r<obsr))+1
  if(obsr<0) ans=cup[loc]
  if(obsr>=0) ans=1-cup[loc]
  return(ans)}
