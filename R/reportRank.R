#' Function to report ranks of 13 criteria for practical significance
#'
#' This function generates a report based on the
#' regression of y on bigx.
#' It acknowledges that some methods for evaluating the importance
#' of  regressor in explaining y may give the importance value
#' with a wrong (unrealistic) sign. For example, m2 reports t-values. Imagine
#' that due to collinearity m2 value is negative when the correct
#' sign from prior knowledge of the subject matter is that the
#' coefficient should be positive, and hence the t-stat should be positive.
#' The wrong sign means the importance of regressor in explaining y
#' should be regarded as relatively less important. The larger the
#' absolute size of the t-stat, the less its true importance in
#' measuring y. The ranking of coefficients computed here
#' suitably deprecates the importance of the regressor
#' when its coefficient has the wrong sign (perverse
#' direction).
#'
#'
#'
#' @param y { (T x 1) vector of dependent variable data values}
#' @param bigx { (T x p) data marix of xi regressor variables associated
#'  with the regression}
#' @param yesLatex {default 1 means print Latex-ready Tables}
#' @param yes13 {default vector of ones to compute all 13 measures.}
#' @param bsign {A (p x 1) vector of right signs of regression
#' coefficients. Default is bsign=0 means the right sign is the same
#' as the sign of the covariance, cov(y, xi) }
#' @param dig {digits to be printed in latex tables, default, dig=d33}
#' @param verbo {logical to print results by pracSig13, default=FALSE}
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @importFrom stats coef cor lm median residuals var
#' @seealso \code{\link{pracSig13}}
#' @note The machine learning methods are subject to random seeds.
#' For some seed values, m10 values from NNS.boost() rarely become
#' degenerate and are reported as NA or missing. In that case
#' the average ranking output r613 here needs adjustment.
#' @return
#'  \item{v15}{practical significance index values (sign adjusted)
#'  for m1 to m5 using older linear and /or bivariate methods}
#'  \item{v613}{practical significance index values
#'  for m6 to m13 newer
#'  comprehensive and nonlinear methods}
#'  \item{r15}{ranks and average rank for m1 to m5
#'  using older linear and /or bivariate methods}
#'  \item{r613}{ranks and average rank for m6 to m13 newer
#'  comprehensive and nonlinear methods}
#'
#'
#' @examples
#' \donttest{
#' set.seed(9)
#' y=sample(1:15,replace = TRUE)
#' x1=sample(2:16, replace = TRUE)
#' x2=sample(3:17, replace = TRUE)
#' x3=sample(4:18,replace = TRUE)
#' reportRank(y,bigx=cbind(x1,x2,x3))
#' }
#'
#' @export


reportRank=function(y, bigx, yesLatex=1, yes13=rep(1,13),
                    bsign=0,dig=3,verbo=FALSE) {
  d3=dig
  p=NCOL(bigx)
  m1=NULL;m2=NULL;m3=NULL;m4=NULL;m5=NULL;m6=NULL;m7=NULL
  m8=NULL;m9=NULL;m10=NULL;m11=NULL;m12=NULL;m13=NULL
  nam=colnames(bigx)
  if (bsign==0) {
    bsign=rep(1,p)
    for (i in 1:p) bsign[i]=sign(cor(y,bigx[,i])) }
  p13=pracSig13(y,bigx,verbo=verbo,yes13=yes13)
  if(verbo){
  print("values before sign adjustment for m1 to m9")
  print("Sign always positive for importance indexes m10 t0 m13")
  print(p13[,1:6],digits=d3)
  print(p13[,7:13],digits=d3)}
  if(verbo){
  if(yesLatex==1){
  print("values before sign adjustment for m1 to m9")
    print(xtable::xtable(p13[,1:6],digits=d3))
    print(xtable::xtable(p13[,7:13],digits=d3)) } }
  sg19=diag(bsign)%*%p13[,1:9]
  PSM=cbind(sg19,p13[,10:13])
  colnames(PSM)=c("m1","m2","m3","m4","m5","m6",
              "m7","m8","m9","m10","m11","m12","m13")
  if(verbo){
    print("values AFTER sign adjustment for m1 to m9")
  print(PSM[,1:6],digits=d3)
  print(PSM[,7:13],digits=d3)
  if(yesLatex==1){
    print("values before sign adjustment for m1 to m9")
  print(xtable::xtable(PSM[,1:6],digits=d3))
  print(xtable::xtable(PSM[,7:13],digits=d3)) }}
n13=NCOL(PSM)
for(j in 1:13) assign(paste("m",j,sep=""),PSM[,j])
#attach(data.frame(PSM))
if(verbo) print(nam)
v15=cbind(m1,m2,m3,m4,m5)  #v=values 15=m1 to m5
rownames(v15)=nam
if(verbo){
print("values of m1 to m5 after sign adjustments if any")
print(v15)
if(verbo){
if(yesLatex==1) print(xtable::xtable(v15,digits=4))
  print("above v15 values of m1 to m5")}}
#now compute rank numbers for m1 to m5
nr=NROW(v15)
rank15=((nr+1)-apply(v15,2,rank,na.last=TRUE, ties.method="average"))
if(verbo) print("ranks values as they are")
avrank15=round(apply(rank15,1,mean),1)
r15=cbind(rank15,avrank15)
if(verbo) {print(r15)
if(yesLatex==1) print(xtable::xtable(r15,digits=2))}
v613=cbind(m6,m7,m8,m9,m10,m11,m12,m13)
#detach(data.frame(PSM))
rownames(v613)=nam
if(verbo) {print(v613)
if(yesLatex==1) print(xtable::xtable(v613,digits=4))}
nr=NROW(v613)
rank613=((nr+1)-apply(v613,2,rank,na.last=TRUE,
                      ties.method="average"))
if(verbo){
print("sign adjusted (if necessary) ranks of p regressors")}
avrank613=round(apply(rank613,1,mean),1)
r613=cbind(rank613,avrank613)
if(verbo){
print(r613)
if(yesLatex==1) print(xtable::xtable(r613,digits=2))
}
list(v15=v15, v613=v613, r15=r15, r613=r613)
}

