\name{summary.genotypeD}
\alias{summary.genotypeD}
\title{genotypeD Object}
\description{creates and object containing allele and genotype frequency for genotypeD objects}
\usage{
summary.genotypeD(object, ...)
}
\arguments{
  \item{object}{genotypeD object}
  \item{...}{further arguments to be passed}
}
\value{
  \item{locus}{locus information field (if present)} 
  \item{allele.names}{vector of allele names: 0 and 1}
  \item{allele.freq}{estimated allele frequencies with finite sample size correction (Lynch \& Milligan 1994)}
  \item{genotype.freq}{frequencies of observed genotypes (phenotypes)}
}
\references{Lynch M. \& Milligan B.G. (1994) Molecular Ecology 3 91-99}
\author{Jarrod Hadfield \email{j.hadfield@sheffield.ac.uk}}
\seealso{\code{genotype}, \code{summary.genotypeD}}
\examples{
l1<-rbinom(100,1,0.5)
l1<-genotypeD(l1)
summary(l1)
}
\keyword{misc}

