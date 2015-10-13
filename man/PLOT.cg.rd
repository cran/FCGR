\name{PLOT.cg}
\alias{PLOT.cg}
\title{Fatigue Crack Growth in Reliability plots}
\usage{
PLOT.cg(x)
}
\arguments{
\item{x}{cracks.growth object}
}
\value{
Return the following values: \item{plot.data}{Exploratory chart.} 
\item{plot.pred}{Plot for fatigue lifetimes estimates and predictions.} 
\item{plot.F}{Plot for the empirical distribution and lifetimes 
distribution estimates of fatigue lifetimes.} \item{plot.resid}{Residuals chart.}
}
\description{
It provides graphical outputs composed of the trends corresponding to 
the crack length growth due to mechanical fatigue, the crack length 
estimates by the models, crack length predictions, and lifetime 
distribution estimates. 
}
\details{
Specifically, the following graphs are provided: exploratory dataset graph, 
plot with the crack length estimates and predictions, residuals graph, 
empirical and estimated lifetime distribution plot obtained by SEP-lme_bkde, 
SEP-lme_kde or PB-nlme methods.}
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
cg <- cracks.growth (x, aF, T_c, method = c("SEP-lme_bkde", "SEP-lme_kde",
                     "PB-nlme"), nBKDE = 5000, nKDE = 5000, nMC = 5000) 
## PLOT.cg applied to cg object.
PLOT <- PLOT.cg(cg) 
names(PLOT)
## [1]  "plot.data"  "plot.pred"  "plot.F"     "plot.resid"
## Exploratory chart for the Alea.A dataset
PLOT$plot.data(main = "Plot:  crack growth", xlab = "million cycles",
                             ylab = "cracks(inches)",  cex.lab=1.8,
                             cex.main = 2)
text(0.02, x$a.F + 0.05, "Failure", cex = 1.8)
text(0.095, 0.95, "Censoring time->", cex = 1.5)
## Plot for fatigue lifetimes estimates and predictions.
PLOT$plot.pred(xlab = "million cycles", ylab = "cracks(inches)",
 main = "Plot: crack growth, estimation and prediction\n failure times (red)",
 cex.lab = 1.8, cex.main = 1.5)
text(0.02,x$a.F+0.05, "Failure", cex = 1.8)
text(0.085,0.95, "Censoring time->", cex = 1.5)  
## Plot for the empirical distribution and lifetimes distribution estimates 
## of  fatigue lifetimes
PLOT$plot.F(main = "Plot: distributions of failure times",
             xlab = "million cycles", ylab = "probability",
             cex.lab = 1.7, cex.main=2)
text(0.14, 0.1, "<-Censoring time", cex = 1.5)
legend("topleft", c("Empirical", "Estimated"), col = c("blue","black"),
        pch=c(20,20), cex=1.5, bty="n")
## Residuals chart.
PLOT$plot.resid(main = "Plot: residual", xlab = "fitted", ylab = "residuals",
                 cex = 1.5, col = "blue", cex.lab = 1.7, cex.main = 2)
}
}
\author{
Antonio Meneses \email{antoniomenesesfreire@hotmail.com}, Salvador
Naya \email{salva@udc.es}, Javier Tarrio-Saavedra
\email{jtarrio@udc.es}, Ignacio Lopez-Ullibarri \email{ilu@udc.es}
}

\references{
Meeker, W., Escobar, L. (1998) Statistical Methods for Reliability 
Data. John Wiley & Sons, Inc. New York.\cr

Pinheiro JC., Bates DM. (2000) Mixed-effects models in S ans S-plus. 
Statistics and Computing.  Springer-Verlang. New York.\cr

Paris, P.C. and Erdogan, F. (1963) A critical analysis of crack 
propagation laws. J. Basic Eng., 85, 528.\cr
}
\keyword{PLOT.cg()}
