\name{varPed}
\alias{varPed}
\title{Transforms Variables for a Multinomial Log-Linear Model}
\description{Creates offspring specific design matrices the columns of which refer to the explanatory variables of the liner model.
}
\usage{
varPed(x, gender=NULL, lag=c(0,0), relational=FALSE, 
  lag_relational=c(0,0), restrict=NULL, keep=FALSE, 
  USvar=NULL, merge=FALSE, ...)
}
\arguments{
  \item{x}{predictor variable; numeric or factor}
  \item{gender}{the gender of the parent to which \code{x} applies}
  \item{lag}{numeric vector of length 2. The time interval over which \code{x} is evaluated relative to a record of the offspring.}
  \item{relational}{a character string. If "Offspring", calculates the distance between \code{x} in the parents and \code{x} in the offspring. If "Mate", calculates the distance between \code{x} in the two parental sexes. If \code{NULL}, \code{x} is untransformed.} 
  \item{lag_relational}{numeric vector of length 2. If \code{relational} is not empty then the time interval over which \code{x} is evaluated in the relational category relative to the offspring record can be changed}
  \item{restrict}{Designates parents with a zero prior probability of parentage. logical or character string. If \code{relational="Offspring"} then \code{restrict} may be \code{TRUE} or \code{FALSE}, and parents that and \code{x} is a \code{factor} then only parents with zero distance are retained as plausible candidates. If a character string only parents for which \code{x} matches are retained.}
  \item{keep}{logical; if \code{TRUE} then the design matrices for parents excluded using the argument \code{restrict} are retained in the estimation of beta}
  \item{USvar}{if \code{NULL}, unsampled parents take the population mean value of the parameter.  If \code{x} is a factor then \code{USvar} can be a level of that factor to which unsampled parents belong.  If \code{x} is numeric then \code{USvar} can be the value for unsampled parents. Is not implemented for \code{relational=="MATE"} or interactions between male and female variables.}  
  \item{merge}{logical; if \code{TRUE} then beta is the log odds ratio of an offspring's parent belonging to category \eqn{A} compared to category \eqn{B}, where \eqn{A} and \eqn{B} are levels of \code{x}.  If \code{FALSE} then beta is the log odds ratio of an individual belonging to category \eqn{A} being the parent of an offspring compared to an individual of category \eqn{B}.  Is not implemented for \code{relational=="MATE"} or interactions between male and female variables, and is only valid for  categorical variables.}
  \item{...}{further arguments to be passed}
}
\details{
  The design matrix for each offspring represents the state of each parental (dam/sire) combination for each explanatory variable. The number of rows in the design matrix (the number of parental combinations) is free to vary across offspring, but the number of explanatory variables remain the same.  As with standard generalised linear modelling the columns of the design matrices take on numerical values or inidicator values for continuous and categorical variables, respectively.  When \code{relational=FALSE}, elements of the design matrices refering to specific parental combinations will not vary across offspring (unless longitudinal data are being used) and the associated vector of parameters will relate the explanatory variables to overall fecundity. For these variables the model is essentially the multinomial analogue of the more familiar Poisson model often used to analyse such data.  However, the counts of the multinomial are not known with certainty because uncertainty exists around the maternity and/or paternity of each offspring.

  Additional variables can fitted that relate specific parental combinations to specific offspring, or specific dams to specific sires.  Elements of the design matrices refering to specific parental combinations are then free to vary across offspring.  The most obvious variable of this type is the mendelian transition probability obtained from the genetic data themsleves. However, by specifying \code{relational="OFFSPRING"}, or  \code{relational="MATE"} non-genetic variables are free to vary across offspring. When \code{x} is \code{numeric} the Euclidean distances between parents and offspring, or between mates enter into the design matrix.   When \code{x} is a \code{factor} then an indicator variable is set up indicating whether parent and offspring, or mate, factor levels match.  Often, each offspring will have a variable number of candidate parents as some parents may be excluded \emph{a priori}. When \code{x} is a \code{factor} and both \code{relational="OFFSPRING"} and \code{restrict=TRUE}, only those potential parents that have factor levels matching the offspring factor level are retained.  When \code{relational=FALSE}, \code{restrict} can take on factor levels which exclude parents that have non-matching factor levels. 

If a time variable (\code{timevar}) is not passed to \code{PdataPed} the data are assumed to be cross-sectional and each indivdiual only respresented once.  If a time variable (\code{timevar}) is passed to \code{PdataPed} then \code{lag} and \code{lag_mate} can be set so that time specific covariates are used.  \code{lag} designates time units realtive to the offspring record when \code{relational=FALSE}; for example, if \code{lag=c(0,0)} the value of \code{x} is taken for that parent during the same time period as the offspring record.   If \code{relational="OFFSPRING"} or \code{relational="MATE"} then \code{lag} determines the time units relative to the record of the offspring or mate to which the focal inidvidual is being compared.  This record can be specified by using \code{lag_relational}, which is always relative to the offspring record. Negative lags refer to previous time intervals (e.g. \code{lag=c(-1,-1)} takes \code{x} from the previous time step), and if the elements of \code{lag} or \code{lag_relational} differ then the average value of \code{x} during this period is taken (e.g \code{lag=c(-1,0)} averages \code{x} in the record matching and preceding the offspring record).  This is not applicable when \code{x} is a \code{factor} unless \code{restrict=TRUE} in which case parents are retained when factor levels match for any times in the specified interval.\ 

Below are models that can be fitted using \code{varPed}, where \code{x} is a univariate continuous variable:

\code{varPed(x, gender="Female")}
\deqn{p^{(o)}_{i,j} \propto \textrm{exp}(\beta_{1}x_{i}...)}{p(i,j) = exp(b*x(i)...)}

\code{varPed(x, gender="Male")}
\deqn{p^{(o)}_{i,j} \propto \textrm{exp}(\beta_{1}x_{j}...)}{p(i,j) = exp(b*x(j)...)}

\code{varPed(x)}
\deqn{p^{(o)}_{i,j} \propto \textrm{exp}(\beta_{1}(x_{i}+x_{j})...)}{p(i,j) = exp(b*(x(i)+x(j))...)}

\code{varPed(x, gender="Female", relational="OFFSPRING")}
\deqn{p^{(o)}_{i,j} \propto \textrm{exp}(\beta_{1}(|x_{i}-x_{o}|)...)}{p(i,j) = exp(b*abs(x(i)-x(o))...)}

\code{varPed(x, gender="Female", relational="MATE")}
\deqn{p^{(o)}_{i,j} \propto \textrm{exp}(\beta_{1}(|x_{i}-x_{j}|)...)}{p(i,j) = exp(b*abs(x(i)-x(j))...)}

\code{varPed(x, gender="Female", lag=c(-1,-1))}
\deqn{p^{(o)}_{i,j} \propto \textrm{exp}(\beta_{1}x_{i,t-1}...)}{p(i,j) = exp(b*x(i,t-1)...)}

\code{varPed(x, gender="Female", lag=c(-1,-1), relational="OFFSPRING")}
\deqn{p^{(o)}_{i,j} \propto \textrm{exp}(\beta_{1}(|x_{i,t-1}-x_{o,t}|)...)}{p(i,j) = exp(b*abs(x(i,t-1)-x(o,t))...)}

\code{varPed(x, gender="Female", lag=c(0,0), relational="MATE",}

  \code{   lag_relational=c(-1,-1))}
\deqn{p^{(o)}_{i,j} \propto \textrm{exp}(\beta_{1}(x_{i,t}+x_{j,t-1})...)}{p(i,j) = exp(b*(x(i,t)+x(j,t-1))...)}

For a categorical variable with two levels (\code{A} and \code{B}) the model specified by \code{varPed(x, gender="Female")} takes on the form

\deqn{p^{(o)}_{i,j} \propto \textrm{exp}(\beta_{1}\delta_{i}...)}{p(i,j) = exp(b*I(i)...)} 

where \eqn{\delta_{i}}{I(i)} is an indicator variable taking the value 1 if \eqn{x_{i}}{x(i)} is equal to the first level of \code{x} and zero otherwise. \eqn{\beta_{1}}{beta} is then the log odds ratio of the two levels of \code{x} with respect to maternity.  If \code{merge=TRUE} is specified then \eqn{\beta_{1}}{beta} may vary across offspring, and \eqn{\beta_{m}}{b_m} is estimated. \eqn{\beta_{m}}{b_m} is related to \eqn{\beta_{1}}{b}:

\deqn{\beta_{m}  = \textrm{logit}\left[\frac{\theta N_{A}}{\theta N_{A} + (1-\theta)N_{B}}\right]}{b_m = logit(((theta*N_A)/(N_A*theta+N_B*(1-theta)))}  

where \eqn{\theta}{theta} is the inverse logit transformation of \eqn{\beta_{1}}{b}, and \eqn{N_{A}}{N_A} and \eqn{N_{B}}{N_B} are the number of potential mothers that have level \code{A} and \code{B} for \code{x}. If \eqn{N_{A}}{N_A} and \eqn{N_B}{N_B} are invariant over offspring the models are functionally equivalent.


The denominator of the multinomial likelihood is the summed linear predictors of all possible parents (after setting up a contrast with the baseline parents).  Designating the first set of parents as baseline, the contrast for each set of parents is simply:

\deqn{\eta^{(o)}_{i,j}  = \textrm{log}\left[\frac{p^{(o)}_{i,j}}{p^{(o)}_{1,1}}\right]}{eta(i,j)=log(p(i,j)/p(1,1))}  

and the likelihood of \eqn{\beta}{b} is 

\deqn{Pr(x|\bf{\beta})  = \prod^{n_{o}}_{o}\left[\frac{\textrm{exp}(\eta^{(o)}_{d,s})}{\sum^{n^{(o)}_{i}}_{i=1}\sum^{n^{(o)}_{j}}_{j=1}\textrm{exp}(\eta^{(o)}_{i,j})}\right]}{Pr(x| b) = prod(no)(exp(eta(d,s))/sum(ni*nj)(exp(eta(i,j))))}  

where \eqn{n_{o}}{no}, \eqn{n^{(o)}_{i}}{ni} and \eqn{n^{(o)}_{j}}{nj} are the number of offspring, the number of potential mothers for offspring \eqn{o}{o}, and the number of potential fathers for offspring \eqn{o}{o}, respectively.  \eqn{d}{d} and \eqn{s}{s} are the actual parents of offspring \eqn{o}{o}. The set of possible parents  in the denominator of the multinomial likelihood are those that are not excluded using the argument \code{restrict}. However, if the argument \code{keep=TRUE} is used then the denominator of the likelihood will include excluded parents depsite the fact that \eqn{d\neq i}{d!=i} and \eqn{s\neq j}{s!=j}.
}
\value{
  list containing the design matrix for variable \code{x}, the identity of retained parents and the gender of the parents}
}
\references{Hadfield J.D. \emph{et al} (2006) Molecular Ecology 15 3715-31}
\author{Jarrod Hadfield \email{j.hadfield@sheffield.ac.uk}}
\seealso{\code{\link{MCMCped}}}
\keyword{models}
\keyword{manip}
