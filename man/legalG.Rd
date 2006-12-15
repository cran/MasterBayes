\name{legalG}
\alias{legalG}
\title{Legal Genotype Configurations}
\description{
A function for checking whether a set of genotypes have a positive probability given the pedigree. If not, a legal configuration is found using heuristic methods.  Missing genotypes are also replaced with compatible genotypes.
}
\usage{
legalG(G, A, ped, time_born=NULL, marker.type="MS", ...)
}
\arguments{
  \item{G}{list of \code{genotype} objects}
  \item{A}{list of allele frequencies}
  \item{ped}{pedigree with id in the first column, dam in the second, and sire in the third. The genotypes must be in the same order as the id column}
  \item{time_born}{an optional vector for ordering a pedigree more efficiently (see \code{\link{order.ped}})} 
  \item{marker.type}{"MS" or "AFLP" for codominant or dominant markers respectively}
  \item{...}{further arguments to be passed}
}
\value{
  \item{G}{a list of \code{genotype} objects with positive likelihood given the pedigree}
  \item{legal}{logical; TRUE if the the genotype configuration passed to \code{legalG} had a positive likelihood}
}
\author{Jarrod Hadfield \email{j.hadfield@sheffield.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
data(WarblerG)
A<-extractA(WarblerG[,16:17])

pedigree<-matrix(NA, 8,3)
pedigree[,1]<-1:8
pedigree[,2][5:8]<-c(1,1,2,2)
pedigree[,3][5:8]<-c(3,3,4,4)

G<-simgenotypes(A, E1=0, E2=0.3, pedigree=pedigree, no_dup=1)

newG<-legalG(G=G$Gobs,A=A,ped=pedigree)
newG$valid

# The input genotypes had a zero probability given the pedigree 
# (because of genotype error) but the outout genotypes have 
# positive probability

legalG(newG$G,A,pedigree)$valid
}
\keyword{manip}
\keyword{models}
