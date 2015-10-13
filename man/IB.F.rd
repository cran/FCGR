\name{IB.F}
\alias{IB.F}
\title{Bootstrap confidence bands for fatigue lifetime}
\usage{
IB.F(z, nB, alpha = 0.05, method = c("SEP-lme_bkde", "SEP-lme_kde",
     "PB-nlme"))
}
\arguments{
\item{z}{cracks.growth object.}
\item{nB}{Number of bootstrap resampling.}
\item{alpha}{Confidence level.}
\item{method}{Character string showing the distribution estimates 
method: "SEP-lme_bkde", "SEP-lme_kde" or "PB-nlme.}
}
\value{
Return a list with the following values:
\item{Mat.F.B }{Matrix that contents the fatigue lifetimes corresponding
to each bootstrap resampling.}
\item{I.Bootstrap}{Data frame that contents the bootstrap confidence bands
for lifetime distribution, at a confidence level of 95 percent (by default). It is 
composed by two columns corresponding to the bands limits: low,  up.}
}
\description{
It performs bootstrap confidence bands for fatigue lifetime. The lifetime 
matrix is calculated by bootstrap resampling by means the above mentioned 
methodologies (see craks.growth). The confidence bands are estimated by 
the quantile based method.
}
\details{
IB.F is performed from the output of cracks.growth function.
}
\examples{
\dontrun{
## Using the Alea.A dataset
data(Alea.A)
x <- Alea.A
## Critical crack length
aF <- 1.6    
## Censoring time         
T_c <- 0.12           
## cracks.growth function applied to Alea.A data
cg <- cracks.growth(x, aF, T_c, method = c("SEP-lme_bkde", "SEP-lme_kde",
                     "PB-nlme"), nBKDE = 5000, nKDE = 5000, nMC = 5000) 
## z is a cracks.growth object
z <- cg
## Number of bootstrap resamplings
nB <- 100 
## Application of IB.F function to cg object
ic.b <- IB.F(z, nB, alpha = 0.05, method = c("SEP-lme_bkde", "SEP-lme_kde",
           "PB-nlme")) 
## ic.b values obtainde by the "SEP-lme_bkde" model
names(ic.b)
# [1] "Mat.F.B"     "I.Bootstrap" 
## Chart with the empirical and estimated distribution functions,
## with bootstrap confidence bands at 95%.
# Observations from which the distribution function is estimated 
F1.F <- z$F.est[,2]  
plot( ic.b$I.Bootstrap$low,F1.F, col=2, type="l", lty=2, lwd=2, 
          xlim=c(0.05,0.18),
          main="Plot: distributions of failure times\n  confidence intervals",
          xlab="million cycles",  ylab="probability",  cex.lab=1.7,
          cex.main=2, las=1)
lines(ic.b$I.Bootstrap$up, F1.F, col=2, lty=2, lwd=2)
points(z$F.est, pch=20)
points(z$F.emp, col=4, pch=20, cex=1.5)
legend("topleft", c("Empirical", "Estimated","Bootstrap (95 percent)"),
              col=c("blue","black","red"),  lty=c(1,1,1), pch=c(20,20,20),
              cex=1.5, bty="n")
## Graph with confidence bands
matplot(ic.b$Mat.F.B, F1.F,  main="Bootstrap resampling lines",  
         type="l", lwd=2, xlim=c(0.05,0.18), xlab="million cycles", 
         ylab="probability", cex.lab=1.7,  cex.main=2, las=1)
}
}
\author{
Antonio Meneses \email{antoniomenesesfreire@hotmail.com}, Salvador
Naya \email{salva@udc.es}, Javier Tarrio-Saavedra
\email{jtarrio@udc.es}, Ignacio Lopez-Ullibarri \email{ilu@udc.es}  
}
\references{
Meeker, W., Escobar, L. (1998) Statistical Methods for Reliability Data. 
John Wiley & Sons, Inc. New York.\cr

Pinheiro JC., Bates DM. (2000) Mixed-effects models in S ans S-plus. 
Statistics and Computing.  Springer-Verlang. New York.\cr

Paris, P.C. and Erdogan, F. (1963) A critical analysis of crack 
propagation laws. J. Basic Eng., 85, 528.\cr
}

\keyword{IB.F}
