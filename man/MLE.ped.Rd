\name{MLE.ped}
\alias{MLE.ped}
\title{Maximum Likelihood Estimation of the Pedigree}
\description{  
Finds the MLE pedigree using the genetic data only.  An approximation is used for genotyping error.
}
\usage{
MLE.ped(X.list, ped=NULL, USdam=FALSE, nUSdam=NULL, USsire=FALSE, 
   nUSsire=NULL, threshold=0, ...)
}
\arguments{
  \item{X.list}{list of design matrices for each offspring derived using \code{\link{getXlist}}}
  \item{ped}{optional pedigree with id, dam and sire in ech column}
  \item{USdam}{logical or character; if \code{TRUE} a single undiferentiated population of unsampled females exists. If \code{USdam} is a character vector it must have the same length as \code{id} with factor levels representing sub-populations (in time or space) over which the number of unsampled females vary.}  
  \item{nUSdam}{numeric vector for number of unsampled females}
  \item{USsire}{logical or character; if \code{TRUE} a single undiferentiated population of unsampled males exists. If \code{USsire} is a character vector it must have the same length as \code{id} with factor levels representing sub-populations (in time or space) over which the number of unsampled males vary.}  
  \item{nUSsire}{numeric vector for number of unsampled males}
  \item{threshold}{threshold probability under which ML parents are replaced by NA}
  \item{...}{further arguments to be passed}
}
\details{
   ML estimation of the pedigree is based on the Mendelian transition probabilities in the presence of genotyping error as outlined in Kalinwoski (2006). The probability that the ML parents are the true parents is simply the Mendelian transition probability for those parents divided by the sum of the transition probabilities for the remaining potential parents, both sampled and unsampled.  If \code{ped} exists and the dam column contains known dam assignemnts and the sire column contains only NA's, then the ML sires will be returned conditional on the dam assignements being true. ML dam estimation with known sires can be performed in the same way.  Individuals whose parents cannot be assigned with the required level of certainty (\code{threshold}), or whose parents belong to the base or unsampled population, have NA in the dan and sire columns. 
}
\value{
 pedigree with id in the first column, and ML estimates for dam and sire in the second and third columns
}
\references{
Hadfield J.D. \emph{et al} (2006) Molecular Ecology 15 3715-31
Marshall J.D. \emph{et al} (1998) Molecular Ecology 7 639-655
Kalinowski S.T. \emph{et al}, Molecular Ecology \emph{in press}
}
\author{Jarrod Hadfield \email{j.hadfield@ed.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\examples{
data(WarblerP)
data(WarblerG)

GdP<-GdataPed(WarblerG)

res1<-expression(varPed("offspring", restrict=0))
res2<-expression(varPed("terr", gender="Female", relational="OFFSPRING", 
  restrict=TRUE))

PdP<-PdataPed(formula=list(res1,res2), data=WarblerP, USsire=TRUE)

X.list<-getXlist(PdP=PdP, GdP=GdP, E2=0.005)

ped<-MLE.ped(X.list, USsire=TRUE, nUSsire=10, threshold=0.75)
}
\keyword{models}
\keyword{optimize}
