#' Compute Effect Sizes for continuous or categorical data
#'
#' Psychologists' so-called "effect size" reveals
#' the practical significance of only one
#' regressor. This function generalizes their algorithm
#' to two or more regressors (p>2). Generalization first
#' converts the xi regressor into a categorical treatment variable
#' with only two categories. One imagines that observations
#' larger than the
#' median (xit> median(xi)) are "treated," and those
#' below the median are "untreated."
#' The aim is the measure the size of the
#' (treatment) effect of (xi) on y. Denote other variables
#' with postscript "o" as (xo). Since we have p regressors in
#' our multiple regression, we need to remove the nonlinear
#' kernel regression effect of
#' other variables (xo) on y while focusing on the effect of xi.
#' There are two options in treating (xo) (i) letting xo be
#' as they are in the data (ii) converting xo to binary
#' at the median. One chooses the first option (i) by setting the
#' logical argument ane=TRUE in calling the function.
#' ane=TRUE is the default. Set ane=FALSE for the second option.
#'
#' @param y { (T x 1) vector of dependent variable data values}
#' @param bigx { (T x p) data matrix of xi regressor variables associated
#'  with the regression}
#' @param ane {logical variable controls the treatment of other regressors.
#'  If ane=TRUE (default), other regressors are used in kernel regression
#'  without forcing them to be binary variables. When ane=FALSE,
#'  the kernel regression removes the effect of other regressors
#'  when other regressors are also binary type categorical variables,}
#' @return out vector with p values of t-statistics for p regressors
#' @note The aim is to answer the following question.
#' Which regressor has the largest
#' effect on the dependent variable? We assume that the signs
#' of regressors are already adjusted such that a numerically
#' larger effect size suggests that the corresponding regressor
#' is most important, having the largest effect size in explaining
#' y the dependent variable.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso \code{\link{pracSig13}}
#' @importFrom generalCorr kern
#' @examples
#' set.seed(9)
#'  y=sample(1:15,replace = TRUE)
#'  x1=sample(2:16, replace = TRUE)
#'  x2=sample(3:17, replace = TRUE)
#' effSizCut(y,bigx=cbind(x1,x2),ane=TRUE)
#'
#' @export


effSizCut=function(y,bigx, ane=TRUE){  #get t-stats
  p=NCOL(bigx)
  out=rep(NA,p)
bigx2=apply(bigx,2,fncut)
logi=apply(bigx2,2,as.logical)
  for ( i in 1:p){
  kxi=generalCorr::kern(dep.y=y, reg.x=bigx2[,i],residuals=TRUE)#r=resid
  rxi=residuals(kxi)#resid= (y - yhat) so yhat=(y-resid)
  xihat=y-rxi
if(!ane){kother=generalCorr::kern(dep.y=y, reg.x=bigx[,-i],residuals=TRUE)
  rother=residuals(kother)
  xotherhat=y-rother }
if(ane) xotherhat=rep(mean(y,na.rm=TRUE),length(y))
  mylogi=logi[,i]
  myxi=xihat[mylogi]
  xibar=mean(myxi,na.rm=TRUE)
  xivar=var(myxi,na.rm=TRUE)
  myxo=xotherhat[mylogi]
  xobar=mean(myxo,na.rm=TRUE) #o=other
  xovar=var(myxo,na.rm=TRUE)
  tim1=length(myxi)-1  #Ti-1, where ti=Ti, m1=minus 1
#  print(c("i,Ti-1= ",i,tim1))
  if(xivar<10e-9) denom=1 #var=variance
  if(xovar<10e-9) denom=1
#  if(denom==1) print(c("t-stat=NA, zero var. col.No.",i))
  if(denom !=1) denom=sqrt(xivar/tim1+xovar/tim1)
  #print(xivar,xovar)

  out[i]=(xibar-xobar)/denom  }#end for loop
  return(out)
}

