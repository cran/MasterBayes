\name{fillX.G}
\alias{fillX.G}
\title{Mendelian Transition Probalities}
\description{
   This function is primarily intended for use within \code{\link{getXlist}}, and fills in the design matrices of the model with the genetic likelihoods.  The model of genotyping error is taken from CERVUS (Kalinowski, 2006; Marshall, 1998).
}
\usage{
fillX.G(X.list, A, G, E2=0.005)
}
\arguments{
  \item{X.list}{list of design matrices for each offspring derived using \code{\link{getXlist}}}
  \item{A}{list of allele frequencies}
  \item{G}{list of genotype objects; rows must correspond to individuals in the vector \code{X.list$id}}
  \item{E2}{per-allele genotyping error rate. \code{E2}(2-\code{E2}) is the per-genotype rate defined in Kalinowski (2006)}
  \item{...}{further arguments to be passed}
}
\value{
  list of design matrices of the form \code{X.list} containing genetic likelihoods for each offspring.
}
\references{Kalinowski S.T. \emph{et al} (2006) Molecular Ecology \emph{in press}}
\author{Jarrod Hadfield \email{j.hadfield@sheffield.ac.uk}}
\note{If a  \code{\link{GdataPed}} object is passed to \code{\link{getXlist}} then the genetic likelihoods will be calculated by default.
}
\seealso{\code{\link{getXlist}}}
\examples{
data(WarblerG)
A<-extractA(WarblerG)

ped<-matrix(NA, 5,3)
ped[,1]<-1:5
ped[,2]<-c(rep(NA, 4), 1)
ped[,3]<-c(rep(NA, 4), 2)

genotypes<-simgenotypes(A, ped=ped)

sex<-c("Female", "Male", "Female", "Male","Female")
offspring<-c(0,0,0,0,1)

data<-data.frame(id=ped[,1], sex, offspring)

res1<-expression(varPed(x="offspring", restrict=0))

PdP<-PdataPed(formula=list(res1), data=data)
GdP<-GdataPed(G=genotypes$Gobs, id=genotypes$id)

X.list<-getXlist(PdP)
# creates design matrices for offspring (in this case indivdiual "5")

X.list.G<-fillX.G(X.list, A=A, G=genotypes$Gobs, E2=0.005)
# genetic likelihoods are arranged sires within dams 

X.list.G$X$"5"$dam.id
X.list.G$X$"5"$sire.id

# so for this example we have parental combinations 
# ("1","2"), ("1","4"), ("3","2"), ("2","4"):

X.list.G$X$"5"$G

# The true parents have the highest likelihood in this case

}
\keyword{manip}
