\name{extractA}
\alias{extractA}
\title{Allele Frequencies}
\description{extracts allele frequencies from genotype data}
\usage{extractA(G, marker.type="MS", ...)}
\arguments{
  \item{G}{data frame or list of \code{genotype} objects}
  \item{marker.type}{optional "MS" or "AFLP" for codominant or dominant markers respectively}
  \item{...}{further arguments to be passed}
}
\value{list of allele frequnecies at each loci}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{genotype.list}}, \code{genotype}}
\examples{
  data(WarblerG)

  A<-extractA(WarblerG)
  A[[1]]
}
\keyword{manip}
