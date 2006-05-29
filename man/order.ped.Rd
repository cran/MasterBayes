\name{order.ped}
\alias{order.ped}
\title{Orders a Pedigree}
\description{Orders a pedigree so parents come before offspring}
\usage{
order.ped(ped, time_born=NULL, ...)
}
\arguments{
  \item{ped}{pedigree with id, dam and sire in ech column}
  \item{time_born}{an optional vector of birth dates by which the pedigree can be ordered)} 
  \item{...}{further arguments to be passed}
}
\value{
  an ordered pedigree pedigree with id, dam and sire in each column}
}
\author{Jarrod Hadfield \email{j.hadfield@sheffield.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
pedigree<-matrix(NA, 8,3)
pedigree[,1]<-1:8
pedigree[,2][5:8]<-c(1,1,2,2)
pedigree[,3][5:8]<-c(3,3,4,4)

pedigree<-pedigree[sample(1:8),]

pedigree2<-order.ped(pedigree)
}
\keyword{manip}
