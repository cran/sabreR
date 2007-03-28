# load library
library(sabreR)
# load data ...
data(growth)
# ... and attach it
attach(growth)

# create a model
sabre.model.1<-sabre(proximity~1,
	             case=teacher,
                     first.family="gaussian",
                     first.mass=64,
                     first.scale=0.5)

# create a different model
# the explanatory variables
sabre.model.2<-sabre(proximity~d1+d2+d3+d4-1,
                     case=teacher,
                     first.family="gaussian",
                     first.mass=64,
                     first.scale=0.5,
                     first.sigma=0.25)




