\name{simgenotypes}
\alias{simgenotypes}
\title{Genotype and Genotyping Error Simulation}
\description{
Simulates genotypes given a pedigree and allele frequencies.  Option exists to simulate observed genotypes given Wangs's (2004) model of genotyping error.
}
\usage{
simgenotypes(A, E1 = 0, E2 = 0, pedigree, no_dup = 1, ...)
}
\arguments{
  \item{A}{list of allele frequencies at each locus}
  \item{E1}{Allelic Dropout Rate}
  \item{E2}{Stochastic Error Rate}
  \item{pedigree}{pedigree in 3 columns: id, dam, sire. Base individuals have NA as parents. All parents must be in id, and each indivdiual must either have both parents in id, or both parents as base.}
  \item{no_dup}{integer: number of times genotypes are to be observed}
  \item{...}{Further arguments to be passed}
}
\value{
  \item{G}{list of genotype objects; true genotypes for each locus}
  \item{Gid}{vector of id names indexing \code{G}}
  \item{Gobs}{list of genotype objects; observed genotypes for each locus}
  \item{id}{vector of \code{id} names indexing \code{Gobs}}
}
\references{Wang J.L. 2004 Genetics 166 4 1963-1979}
\author{Jarrod Hadfield \email{j.hadfield@sheffield.ac.uk}}
\seealso{\code{genotype}}
\examples{
pedigree<-cbind(1:10, rep(NA,10), rep(NA, 10))

gen_data<-simgenotypes(A=list(loc_1=c(0.5, 0.2, 0.1, 0.075, 0.025)), 
 E1=0.1, E2=0.1, pedigree=pedigree, no_dup=1)

summary(gen_data$G[[1]])
summary(gen_data$Gobs[[1]])
}
\keyword{datagen}

