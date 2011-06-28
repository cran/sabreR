\name{print.sabre}
\alias{print.sabre}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Defining and Fitting SABRE Models}
\description{
  'print.sabre' is used to print information about a fitted sabre model.
}
\usage{



\method{print}{sabre}(x, estimates = TRUE, iterations = FALSE, variables = FALSE, settings = TRUE, shm = TRUE, rem = TRUE,\dots) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{


  \item{x}{The sabre model to be printed (mandatory)}
  
  \item{estimates}{Output will include model estimates for standard
    homogenous model (if shm = TRUE) and random effects model (if rem =
    TRUE).}
  
  \item{iterations}{Output will include the iterations from the model fitting for standard
    homogenous model (if shm = TRUE) and random effects model (if rem
    = TRUE).}
  
  \item{variables}{Output will include the model variable names and types for standard
    homogenous model (if shm = TRUE) and random effects model (if rem =
    TRUE).}

  \item{settings}{Output will include the model settings and likelihood
    figures for standard  homogenous model (if shm = TRUE) and random effects model (if rem = TRUE).}  
  
  \item{shm}{Enables information regarding the standard homogenous model
    when set to TRUE}

  \item{rem}{Enables information regarding the random effects model
    when set to TRUE} 
  
  \item{...}{other arguments to be passed to or from other methods}

}

\value{
Displays the details of the fitted model.
}

\author{Prof. R. Crouchley
  Centre for e-Science
  Lancaster University
  Lancaster
  United Kingdom
  e-mail : asarc@exchange.lancs.ac.uk}

\examples{

# load data ...
data(drvisits)
# ... and attach it
attach(drvisits)

# the first model
sabre.model.1<-sabre(numvisit~reform+age+educ+married+badh+loginc+summer,
                     case=obs,
                     first.family="poisson")

# the second model
sabre.model.2<-sabre(numvisit~reform+age+educ+married+badh+loginc+summer,
                     case=id,
                     first.family="poisson")

# compare them
print(sabre.model.1)
print(sabre.model.2)

}
\keyword{multivariate} % at least one, from doc/KEYWORDS


