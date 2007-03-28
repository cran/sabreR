# load data
data(tu)
# attach it
attach(tu)
# create the model
sabre.model.1<-sabre(TU~YEAR+AGE+EVNO+SUPR+HRS+factor(NOEM)+SEX1+PROM+factor(SC80)-1,case=CASE)
# examine the results
sabre.model.1
