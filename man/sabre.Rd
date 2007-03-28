\name{sabre}
\alias{sabre}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Defining and Fitting SABRE Models}
\description{
  'sabre' is used to define and fit a SABRE model.
}
\usage{
sabre(model.formula.uni, model.formula.bi = NULL, model.formula.tri = NULL,
    case, alpha = 0.01, approximate = 5, max.its = 100, arithmetic.type = "fast", 
    offset = "", convergence = 5e-05, correlated = "yes", left.end.point = NULL, 
    right.end.point = NULL, first.family = "binomial", second.family = "binomial", 
    third.family = "binomial", first.link.function = "logit", 
    second.link.function = "logit", third.link.function = "logit", 
    first.mass = 12, second.mass = 12, third.mass = 12, ordered = FALSE, 
    first.scale = -10000, second.scale = -10000, third.scale = -10000, 
    first.rho = 0, second.rho = 0, third.rho = 0, first.sigma = 1, 
    second.sigma = 1, third.sigma = 1, tolerance = 1e-06, equal.scale = FALSE, 
    depend = FALSE, only.first.derivatives = FALSE, adaptive.quad = FALSE, 
    fixed.effects = FALSE) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{


  \item{model.formula.uni}{Univariate model formula (mandatory)}
  
  \item{model.formula.bi}{Bivariate model formula}
  
  \item{model.formula.tri}{Trivariate model formula}
  
  \item{case}{List of variables which contain the case structure (one for univariate, two for
    bivariate or three for trivariate (mandatory)}
  
  \item{alpha}{Value of the orthogonality constant, which is used to trigger special action during the model fitting process.
    If the orthogonality criterion (which measures the propensity of the
    step direction to be orthogonal to the gradient) at any iteration is
    less than alpha, then the diagonal of the estimated Hessian matrix is
    doubled (after first ensuring that all elements on the diagonal are
    positive) and the iteration repeated.
    A low value of alpha will cause action to be taken only in conditions of
    bad orthogonality. A high value will cause special action to be taken
    often, and is not usually recommended.}

  \item{approximate}{approximate is the number of iterations (positive integer) which are to be performed using the Meilijson approximation to the Hessian matrix.
    This first derivative based approximation is more robust than the matrix
    of true second derivatives when the parameter estimates are a long way from their
    maximum likelihood solutions. For this reason, the first number
    iterations of any mixture model fit are carried out using the
    approximation to the Hessian.
    If the true Hessian is positive semi-definite on the (number+1)st
    iteration, then the algorithm switches to this matrix. Otherwise, the
    approximate method is
    retained until the true Hessian attains positive semi-definiteness. All further iterations continue to make use of the matrix of true second derivatives.
    
    Note that, even if convergence is achieved within the first number
    iterations, the parameter estimate standard errors are still calculated
    from the true Hessian since
    it alone provides consistent estimates of their values.}
  
  
  \item{max.its}{Sets to number the maximum number of iterations which are to be
    performed by the model fitting algorithm. For mixture models, the
    first approximate iterations use the Meilijson approximation to the
    Hessian, while the remaining
    (max.its - approximate) iterations use the matrix of true second derivatives.}

  \item{arithmetic.type}{Determines which method will be used in the calculation of the likelihood and its derivatives for mixture models.
    If arithmetic.type is "fast", then real numbers are stored in standard
    floating point format, which may potentially give rise to underflow
    problems for long sequences of observations for the same case.
    If arithmetic.type = "accurarte", then real numbers are stored in
    mantissa-exponent format, which prevents underflow but may significantly
    decrease the speed of the fitting algorithm.
    If underflow does occur while using the fast method, an error message will be returned and the fitting abandoned.}

  \item{offset}{Specifies the (quoted) name as an a priori known component to be included
    in the linear predictor during model fitting. The offset itself is an
    explanatory variable - which may be both subject and time period specific - whose coefficient is fixed at 1.0, and is thus constant throughout the fitting process.
    In order to include an offset in a homogeneous model, it is necessary
    to have a single mass point and no end points (see left.end.point,
    right.end.point, first.mass, second.mass and third.mass)}

  <
  \item{convergence}{Sets the upper limit in difference in log-likelihood between subsequent iterations that
    has to be acheived for the model fitting to be deemed as converged. }


  \item{correlated}{Specifies whether a correlated or uncorrelated model is to be fitted. Applies only to bivariate and trivariate models. }

  \item{left.end.point}{Specificies the initial (numeric) value for the left end point. Setting
    left.end.point to NULL (the default) indicates that a left end point is
    not to be used. If a numeric value is specified then a left end point is
    used by the fitting algorithm starting with this initial value. Can only
    be used with univariate binary.}

  \item{right.end.point}{Specificies the initial (numeric) value for the right end point. Setting
    right.end.point to NULL (the default) indicates that a right end point is
    not to be used. If a numeric value is specified then a right end point is
    used by the fitting algorithm starting with this initial value. Can only
    be used with univariate binary or univariate Poisson models.}


  \item{first.family}{Specifies the distribution of the first dependent
  variable for univariate, bivariate or trivariate models. Values can be "binomial", "gaussian" or "poisson".}

  \item{second.family}{Specifies the distribution of the second dependent
  variable for univariate, bivariate or trivariate models. Values can be
  "binomial", "gaussian" or "poisson".}

  \item{third.family}{Specifies the distribution of the third dependent
  variable for trivariate models. Values can be
  "binomial", "gaussian" or "poisson".}  

  \item{first.link.function}{Specifies the first link function for
    univariate, bivariate or trivariate models. Values can be "l"
    (logit), "p" (probit), or "c" (complementary log-log).}

  \item{second.link.function}{Specifies the second link function for
    bivariate or trivariate models. Values can be "l"
    (logit), "p" (probit), or "c" (complementary log-log).}
  
   \item{third.link.function}{Specifies the third link function for
     trivariate models. Values can be "l"
     (logit), "p" (probit), or "c" (complementary).}

   \item{first.scale}{Sets the initial value (positive real number) for
     the first scale parameter in univariate, bivariate or trivariate
     models.}

   \item{second.scale}{Sets the initial value (positive real number) for
     the second scale parameter in bivariate and trivariate
     models.}
      
   \item{third.scale}{Sets the initial value (positive real number) for
     the third scale parameter in trivariate models.}

   \item{first.mass}{Specifies the number of quadrature points for the
     first response in univariate, bivariate or trivariate
     models. Values can be
     1,2,4,6,8,10,12,14,16,20,24,32,48,64,72,80,88,96,104,112,128,144,160,176,193,208,224,240 or 256.}

   \item{second.mass}{Specifies the number of quadrature points for the
     second response in bivariate or trivariate
     models. Values can be 1,2,4,6,8,10,12,14,16,20,24,32,48,64,72,80,88,96,104,112,128,144,160,176,193,208,224,240 or 256.}

   \item{third.mass}{Specifies the number of quadrature points for the
     third response in trivariate
     models. Values can be 1,2,4,6,8,10,12,14,16,20,24,32,48,64,72,80,88,96,104,112,128,144,160,176,193,208,224,240 or 256.}
   
   \item{first.sigma}{Sets the initial value for the residual standard
     deviation (positive real number) in univariate, bivariate or trivariate models where the first
     response is linear.}

   \item{second.sigma}{Sets the initial value for the residual standard
     deviation (positive real number) in bivariate or trivariate models where the second
     response is linear.}

   \item{third.sigma}{Sets the initial value for the residual standard
     deviation (positive real number) in trivariate models where the third
     response is linear.}

   \item{first.rho}{Sets the initial value for the correlation parameter
     in bivariate models, or the first correlation parameter (corr(1,2)) in trivariate models.}

   \item{second.rho}{Sets the initial value for the second correlation parameter (corr(1,3)) in trivariate models.}

   \item{third.rho}{Sets the initial value for the third correlation parameter (corr(2,3)) in trivariate models. }


   \item{ordered}{Specifies an ordered response model if
     ordered="yes". This is an ordered probit
     if the link function is "probit" or an ordered logit if the link function is "logit".}

   
   \item{depend}{If depend=TRUE then specifies a univariate random effects model with two scale parameters (first.scale & second.scale).}
   
   \item{equal.scale}{If TRUE, specifies a bivariate random effects model in
     which both scale parameters are equal, so that the estimated random
     effects parameters in a subsequent model fit are first.scale = second.scale and correlated="yes".}


   \item{only.first.derivatives}{If TRUE, uses only the Meilijson
     approximation to the Hessian matrix (see approximate).}

   \item{adaptive.quad}{If TRUE, the employs adaptive quadrature.}
   
   \item{tolerance}{tolerance (numeric) sets the tolerance used in the
     matrix inversion routines to detect for extrinsic aliasing and lack
     of positive semi-definiteness.}
   
   
   \item{fixed.effects}{If fixed.effects=TRUE then a fixed effects model is fitted. The fixed
   effects model uses implicit dummy variables, with a dummy variable
   included for each unique value of the corresponding variable. The
   parameter estimates for these dummy variables are omitted from the
   output.
   Note that time-constant covariates (including the constant itself)
   should not be included in the model as they cannot be estimated
   and will be aliased out of the model.}
   
}

\value{
Displays the progress of the model fitting. To see the details of the
fitted model use print.
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
\keyword{models} % at least one, from doc/KEYWORDS



