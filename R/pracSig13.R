#' Compute thirteen measures of practical significance
#'
#' Thirteen methods are denoted m1 to m13. Each
#' yields p numbers when there are p regressors denoted xi.
#' m1=OLS coefficient slopes. m2= t-stat of each slope.
#' m3= beta coefficients OLS after all variables have
#' mean zero and sd=1.
#' m4= Pearson correlation coefficient between y and xi (only two
#' variables at a time, assuming linearity).
#' Let r*(y|xi) denote the generalized correlation coefficient allowing
#' for nonlinearity from Vinod (2021, 2022). It does not equal
#' analogous r*(xi|y). The larger of the two,
#' max(r*(y|xi), r*(xi|y)), is given by
#' the function depMeas() from the 'generalCorr' package.
#' m5= depMeas, which allows nonlinearity. m5 is not comprehensive
#' because it measures only two  variables, y and xi, at a time.
#' m6= generalized partial correlation coefficient or
#' GPCC. This is the first comprehensive measure
#' of practical significance.
#' m7=a generalization of psychologists' "effect size" after incorporating
#' the nonlinear effect of other variables.
#' m8= local linear partial (dy/dxi) using the 'np' package
#' for kernel regressions and local linear derivatives.
#' m9= partial derivative (dy/dxi) using the 'NNS' package.
#' m10=importance measure using NNS.boost() function of 'NNS.'
#' m11=Shapley Value measure of importance (cooperative game theory).
#' m12 and m13= two versions of the random forest algorithm
#' measuring the importance of regressors.
#'
#' If m6, m10 slow down computations, we recommend setting
#' yes13[6]=0=yes13[10] to turn off slowcomputation of m6 and m10
#' at least initially to get quick answers for other m's.
#'
#' @param y {input dependent variable data as a vector}
#' @param bigx {input matrix of p regressor variables}
#' @param yes13 {vector of ones to compute respective
#' 13 measures m1 to m13. Default is all ones to compute all
#' e.g., yes13[10]=0 means do not compute the m10 method.}
#' @param verbo {logical to print results along the way default=FALSE}
#' @return output matrix (p x 13) containing m1 to m13
#' criteria (numerical measures of practical significance)
#' along columns
#' and a row for each regressor (excluding the intercept).
#' @note needs the function kern(), which requires package 'np'.
#' also needs 'NNS', 'randomForest', packages.
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also as \code{\link{effSizCut}},
#' @seealso See Also as \code{\link{reportRank}}
#' @importFrom stats coef cor lm median
#' @importFrom stats  residuals var
#' @importFrom utils combn
#' @importFrom generalCorr depMeas parcorVec kern kern2
#' @importFrom NNS dy.d_ NNS.boost
#' @importFrom randomForest randomForest importance
#' @importFrom np npreg npregbw
#' @importFrom ShapleyValue shapleyvalue
#' @references Vinod, H. D."Generalized Correlation and Kernel Causality with
#'  Applications in Development Economics" in Communications in
#'  Statistics -Simulation and Computation, 2015,
#'  \doi{10.1080/03610918.2015.1122048}
#' @keywords GPCC
#' @keywords Shapley Value
#' @keywords Random Forest
#' @note The machine learning methods are subject to random seeds.
#' For some seed values, m10 values from NNS.boost() become
#' degenerate and are reported as NA or missing. In that case
#' the average ranking output r613 from reportRank() needs manual adjustments.
#'
#' @references Vinod, H. D.", "Generalized Correlations and Instantaneous Causality
#'  for Data Pairs Benchmark," (March 8, 2015).
#'   \url{https://www.ssrn.com/abstract=2574891}
#' @references Vinod, H. D. “Generalized, Partial and Canonical
#' Correlation Coefficients,” Computational Economics
#' (2021) SpringerLink vol. 59, pp.1-28.
#' URL https://link.springer.com/article/10.1007/s10614-021-10190-x
#' @references Vinod, H. D. “Kernel regression coefficients for practical
#' significance," Journal of Risk and Financial Management
#' 15(1), 2022 pp.1-13. https://doi.org/10.3390/jrfm15010032
#' @references Vinod, H. D.", "Hands-On Intermediate Econometrics
#' Using R"  (2022) World Scientific Publishers: Hackensack, NJ.
#' \url{https://www.worldscientific.com/worldscibooks/10.1142/12831}
#'
#' @export



pracSig13=function(y,bigx, yes13=rep(1,13),verbo=FALSE){
  if (verbo){
  print("m1=ols coeff values")
  print("Linear regression focusing of statistical significance")}
  importance=NULL
  p=NCOL(bigx)
  options(np.messages=FALSE)
  mall=rep(NA,13)
  for(i in 1:13) mall[i]=paste("m",i,sep="")
  for(i in 1:13) assign(mall[i],rep(NA,p))
  if(verbo) print(c("initialize typical mi, m4=", m4))
  myp=1:p
  nam=colnames(bigx)
  inam=as.character(1:p);inam
  if(length(nam)==0) nam=paste("x",eval(myp),sep="")
  if(verbo)  print(nam)
  reg1=lm(y ~ bigx)
  if (verbo){
    print("compute m1=OLS coeff. values using original data")}
  if (yes13[1]==1){
  m1=coef(reg1)[2:(p+1)]} #omit intercept coeff magnitudes
  if (verbo){
    print(m1)}

  if (verbo){
    print("compute m2=t-stat values")}
  stu1x=summary(reg1)$coef[,3] #third col=t-values
  if (yes13[2]==1){
    m2=stu1x[2:(p+1)]
  if (verbo) print(m2)}

  if (verbo){
    print("m3=OLS regression coefficients on standardized data")}
  if (yes13[3]==1){
  sy=scale(y)
  sbigx=apply(bigx,2,scale) #data are standardized
  if (verbo){
  print("regression forced through origin, so NO intercept")}
  reg2=lm(sy~sbigx-1); m3=coef(reg2)
  if (verbo){
  print("standardized regr coefficients=practical importance beta")
  print(m3)}
  } #endif yes13[3]

  if (yes13[4]==1){
    m4=rep(NA,p)
  if (verbo){
    print("Pearson corr is m4")}
  for (i in 1:p) m4[i]=cor(y,bigx[,i])
  if (verbo){
    print(m4)}
  } # endif yes13[5]

  if (yes13[5]==1){
    m5=rep(NA,p)
  if (verbo){
    print("generalized corr max(r*ij, r*ji) is m5")}
  for (i in 1:p) m5[i]=generalCorr::depMeas(y,bigx[,i])
  if (verbo){
    print(m5)}
  }

  if (yes13[6]==1){
    if (verbo){
    print("m6=generalized partial correl. coeff")}
    if (p<5){
  m6=as.vector(generalCorr::parcorVec(cbind(y,bigx)))}
    if (p >= 5){
      m6=as.vector(generalCorr::parcorVecH(cbind(y,bigx)))}
    if (verbo){
    print(m6)}}

  if (yes13[7]==1){
  if (verbo){
    print("m7=effect size Psychology, force xi=(0,1) at median")}
  m7=effSizCut(y,bigx,ane=TRUE)
  if (verbo){
    print(m7)} }

  if (yes13[8]==1){
    if (verbo){
    print("now m8 for local linear partials dy/dxi using np package")}

  scaly=as.vector(scale(y))
  sbigx=apply(bigx,2,scale)
  m8=rep(NA,p) #place to store
  for ( i in 1:p){
    kxi=generalCorr::kern2(dep.y=y, reg.x=sbigx[,i],gradients=TRUE)
    m8[i]=apply(kxi$grad,2,mean) }
  if (verbo){
    print("np estimate of mean of derivatives standardized units")
  print(round(m8,3))}
} #endif m8 calculations

  if (yes13[9]==1){
    if (verbo){
    print("now m9 dy/dxi using NNS package")}
  scaly=scale(y)
  sbigx=apply(bigx,2,scale)
  nnsout=NNS::dy.d_(x=sbigx, y=scaly[,1], wrt=1:p,
               messages=FALSE, eval.points = "obs")
  lapply(nnsout["First",], mean)
  m9x= lapply(nnsout["First",], mean)
  m9=as.numeric(m9x)
  if (verbo){
    print(m9)}
} #endif m9 calculations

  if (yes13[10]==1 & verbo) {
    print("now m10 boost importance using NNS package")}

    if (yes13[10]==1){
  NNSb=NNS::NNS.boost(IVs.train = bigx, DV.train = y,status=FALSE,
                 feature.importance = FALSE)
  m10=as.numeric(NNSb$feature.weights)
  if(verbo)print(m10)
  }# endif m10

  if (yes13[11]==1){
    if (verbo){
    print("now m11 Shapley measure of importance")}
    bigxx=data.frame(bigx)
  m11=as.numeric(shapleyvalue(y,bigxx)[2,])
  if (verbo)   print(m11)
  } # endif m11

  if (yes13[12]==1 | yes13[13]==1){
  if (verbo){
  print("now m12, m13 random forest measures of importance")}
  y.rf = randomForest(x=bigx,y=y, importance=TRUE)
  imp12=round(importance(y.rf), 3)
  if(yes13[12]==1) m12=as.numeric(imp12[,1])
  if (verbo){
    print(m12)}
  if(yes13[13]==1) m13=as.numeric(imp12[,2])
  if (verbo){
    print(m13)}
} #endif m12 m13

  out=matrix(NA,nrow=p, ncol=13)

  for (j in 1:13){
  if (j==1 & yes13[j]==1)  out[,j]=m1
  if (j==2 & yes13[j]==1)  out[,j]=m2
  if (j==3 & yes13[j]==1)  out[,j]=m3
  if (j==4 & yes13[j]==1)  out[,j]=m4
  if (j==5 & yes13[j]==1)  out[,j]=m5
  if (j==6 & yes13[j]==1)  out[,j]=m6
  if (j==7 & yes13[j]==1)  out[,j]=m7
  if (j==8 & yes13[j]==1)  out[,j]=m8
  if (j==9 & yes13[j]==1)  out[,j]=m9
  if(length(m10)==p){
  if (j==10 & yes13[j]==1)  out[,j]=m10}
  if (j==11 & yes13[j]==1)  out[,j]=m11
  if (j==12 & yes13[j]==1)  out[,j]=m12
  if (j==13 & yes13[j]==1)  out[,j]=m13
  }
  rownames(out)=nam
  return(out)
}
