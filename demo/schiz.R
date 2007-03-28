
# load data ...
data(schiz)
# ... and attach it
attach(schiz)

# the first model
sabre.model.1<-sabre(impso~weeksqrt+treatment+interact-1,
                     case=id,ordered=TRUE)

# examine the result
print(sabre.model.1)
