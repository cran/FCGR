\docType{data}
\name{Alea.A}
\alias{Alea.A}
\title{Crack growth of aluminum A-alloy due mechanical fatigue efforts}

\format{This data frame is composed of the following columns:
\describe{
\item{\code{cycles}}{21 vectors with times or cycles corresponding to each
                     specimen (overall 262 cycles).}
\item{\code{cracks}}{21 vectors with the crack lengths corresponding 
                     to each specimen.}
\item{\code{sample}}{21 different specimens, each one repeated the number 
                     of elements of "cycles" or "cracks".}
}}
\source{
Meeker, W., Escobar, L. (1998) Statistical Methods for Reliability Data.
John Wiley & Sons, Inc. New York.
}
\description{
Crack growth of aluminum A-alloy due mechanical fatigue efforts. Alea.A is 
a data frame composed of 262 rows y 3 columns.
}
\details{
Alea.A is composed of the crack length growth and number of cycles of 21 
different specimens. It is detailed in Hudak et al. (1978) and referred 
by Meeker and Escobar (1998) in table C.14. 
}
\examples{
data(Alea.A)
}
\keyword{datasets}
