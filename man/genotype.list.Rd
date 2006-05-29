 \name{genotype.list}
\alias{genotype.list}
\title{Genotype Objects for all Loci}
\description{
  Creates a \code{list} of \code{genotype} objects from a \code{matrix} or \code{data.frame} of multilocus genotypes.
}
\usage{
genotype.list(G, ...)
}
\arguments{
  \item{G}{matrix or data.frame of multilocus genotypes with individuals down the rows and loci across columns. Adjacent columns are taken to be the same locus}
  \item{...}{further arguments to be passed}
}
\value{
  list of \code{genotype} objects for all loci
}
\author{Jarrod Hadfield \email{j.hadfield@sheffield.ac.uk}}
\seealso{\code{genotype}}
\examples{

  data(WarblerG)

  G<-genotype.list(WarblerG)
  summary(G[[1]])

}
\keyword{manip}
