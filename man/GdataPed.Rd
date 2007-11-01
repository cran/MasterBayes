\name{GdataPed}
\alias{GdataPed}
\alias{is.GdataPed}
\title{GdataPed Object}
\description{An object containing genotype data and the categories over which error rates may vary.}
\usage{
GdataPed(G, id = NULL, categories = NULL, perlocus=FALSE, marker.type="MS", ...)
}
\arguments{
  \item{G}{a list of \code{genotype} objects for each locus, or a \code{data.frame} to be coerced using \code{\link{genotype.list}}}
  \item{id}{a vector of individual identifiers associated with each genotype, indivdiuals can have more than one observed genotype.  If \code{G} is a \code{data.frame} to be coerced and has a column name \code{id}, this will be used.}
 \item{categories}{an optional vetcor indicating subsets of genotypes that have different error rates.  If \code{G} is a \code{data.frame} to be coerced and has a column name \code{categories}, this will be used.}
 \item{perlocus}{if \code{TRUE} different error rates are estimated for each locus} 
 \item{marker.type}{"MS" or "AFLP" for codominant or dominant markers respectively}
 \item{...}{further arguments to be passed}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
data(WarblerG)
GdP<-GdataPed(WarblerG)
}
\keyword{classes}

