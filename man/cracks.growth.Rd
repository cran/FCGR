\name{cracks.growth}
\alias{cracks.growth}
\title{Fatigue crack growth in reliability analysis}
\usage{
cracks.growth(x, aF, T_c, method = c("SEP-lme_bkde", "SEP-lme_kde",
              "PB-nlme"), nMC = 5000, nBKDE = 5000, nKDE = 5000)
}
\arguments{
\item{x}{Matrix or data frame composed of three columns: times or                
number of cycles, crack lengths and specimen number.}

\item{aF}{Critical crack length for which the material failure is produced.}
  
\item{T_c}{Censoring time or frequency.}
  
\item{method}{A string of characters: 
"SEP-lme_bkde" (default methodology) indicates that a mixed 
effects linear regression model is applied to crack growth data 
and the lifetime density is estimated by bkde method.
"SEP-lme_kde" indicates that a mixed effects linear regression 
model is applied to crack growth data and the lifetime 
distribution is estimated by kde method.
"PB-nlme" indicates that a mixed effects nonlinear regression 
model is applied to crack growth data and the lifetime 
the parameters are estimated maximum likelihood 
methodologies, and the lifetime distribution by Monte Carlo.}

\item{nMC}{Number of Monte Carlo estimates, by default 5000.}

\item{nBKDE}{Number of bkde estimates, by default 5000.}

\item{nKDE}{Number of kde estimates, by default 5000.}
}
\value{
Return a list with the following values:

\item{data}{Data frame with the data corresponding to number of cycles, 
crack length, and sample.} \item{a.F}{Critical crack length.} 
\item{Tc}{Censoring time.} \item{param}{Data frame with the estimates 
of Paris law parameters: C and m} \item{crack.est}{Data frame with time, 
crack growth estimates, and corresponding sample or specimen.} 
\item{sigma}{Residual standard deviation.} \item{residuals}{Residuals 
resulting from the crack length fitting.} \item{crack.pred}{Data frame 
with time, crack growth predictions out of the experimental time interval, 
and corresponding sample or specimen.} \item{F.emp}{Data frame with the 
empirical lifetime distribution and the corresponding time: time, Fe. } 
\item{bw}{Bandwidth used in bkde and kde methods.} \item{F.est}{Data frame 
with the estimated lifetime distribution and the corresponding time: time, F.}
\item{nBKDE}{Number of bkde estimates.}
\item{nKDE}{Number of kde estimates.}
\item{nMC}{Number of Monte Carlo estimates.}
}
\description{
It provides the lifetime distribution of metallic materials that fail due
to the crack growth produced by mechanical fatigue efforts. The crack 
growth trends are fitted by linear or nonlinear mixed effects regression 
models in order to make predictions about the material lifetime. The 
lifetime is defined as the time passed before the material does not meet 
the specification requirements and it is conditioned by a critical crack 
length that induces the material failure. Three different methods can be 
applied to estimate the fatigue lifetime distribution: "SEP-lme_bkde" 
and "SEP-lme_kde" are nonparametric while "PB-nlme" corresponds to the 
parametric approach proposed by Pinheiro and Bates (2000).
}
\details{
This function provides a simultaneous fitting of crack growth data 
corresponding to different specimens when these are subjected to 
mechanical fatigue efforts. For this purpose, mixed effects linear 
models (lme) with B-spline smoothing are applied. Since the failure 
is defined at a specific critical crack length, predictions of 
material lifetime are obtained assuming the linearized Paris-Erdogan 
law and the material lifetime distribution is estimated. There are 
available three different techniques to estimate the lifetime 
distribution: the binned kernel density estimate (bkde), the 
kernel estimator for the distribution function (kde) computed 
by Quintela del Rio and Estevez-Perez (2012), and in addition the 
parametric method proposed by Pinheiro and Bates (2000) based on 
mixed effects nonlinear regression (nlme), maximum likelihood and 
Monte Carlo simulation. 
}
\examples{
\dontrun{
## Using the Alea.A dataset
data(Alea.A)
x <-  Alea.A
## Critical crack length
aF <- 1.6    
## Censoring time         
T_c <- 0.12           
## cracks.growth function applied to Alea.A data
cg <- cracks.growth (x, aF, T_c, method = c("SEP-lme_bkde", "SEP-lme_kde",
                   "PB-nlme"), nBKDE = 5000, nKDE = 5000, nMC = 5000) 
## cracks.growth values using the "SEP-lme_bkde" by default method.
names(cg)
# [1]	 "data"       "a.F"        "Tc"         "param"      "crack.est"
# [6] 	"sigma"      "residuals"  "crack.pred" "F.emp"      "bw"
#[11]	 "F.est"      "nBKDE"
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

Quintela-del-Rio, A. and Estevez-Perez, G. (2012) Nonparametric Kernel 
Distribution Function Estimation with kerdiest: An R Package for Bandwidth 
Choice and Applications, Journal of Statistical Software 50(8), 1. 
URL http://www.jstatsoft.org/v50/i08/.\cr
}

\keyword{cracks.growth}
