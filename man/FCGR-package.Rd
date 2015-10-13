\docType{package}
\name{FCGR-package}
\alias{FCGR-package}
\title{Fatigue Crack Growth in Reliability}
\description{
Fatigue Crack Growth in Reliability estimates the distribution of 
material lifetime due to mechanical fatigue efforts. The FCGR package
provides simultaneous crack growth curves fitting to different specimens
in materials under mechanical stress efforts. Linear mixed-effects models
(LME) with smoothing B-Splines and the linearized Paris-Erdogan law are
applied. Once defined the fail for a determined crack length, the distribution 
function of failure times to fatigue is obtained. The density function 
is estimated by applying nonparametric binned kernel density estimate 
(bkde) and the kernel estimator of the distribution function (kde). The 
results of Pinheiro and Bates method based on nonlinear mixed-effects 
regression (nlme) can be also retrieved. The package contains the 
crack.growth, PLOT.cg, IB.F, and Alea.A  (database) functions.}
\details{
\tabular{ll}{ Package: \tab FCGR\cr Type: \tab Package\cr Version: \tab 1.0\cr
Date: \tab 2015-09-29\cr License: \tab GPL >= 2\cr } 
The main functions and data frame are \code{cracks.growth}, \code{IB.F},
\code{PLOT.cg}, \code{Alea.A} 
}
\author{
Antonio Meneses \email{antoniomenesesfreire@hotmail.com}, Salvador
Naya \email{salva@udc.es}, Javier Tarrio-Saavedra
\email{jtarrio@udc.es}, Ignacio Lopez-Ullibarri \email{ilu@udc.es}

Maintainer: Antonio Meneses \email{antoniomenesesfreire@hotmail.com}
}
\references{
Meeker, W., Escobar, L. (1998) Statistical Methods for Reliability 
Data. John Wiley & Sons, Inc. New York.\cr

Pinheiro JC., Bates DM. (2000) Mixed-effects models in S ans S-plus. 
Statistics and Computing.  Springer-Verlang. New York.\cr

Paris, P.C. and Erdogan, F. (1963) A critical analysis of crack 
propagation laws. J. Basic Eng., 85, 528.\cr
}
\keyword{package}
