
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


