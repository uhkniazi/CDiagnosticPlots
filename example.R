# Name: Example.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 25/05/2017
# Desc: example usage of cdiagnostic plots


source('CDiagnosticPlots.R')

# utility function to load object
f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}

## load the test data, some microarray test data
lData = f_LoadObject('lData.test_data.rds')
# it is a list with various components
names(lData)

mCounts = lData$data

## create the object
oDiag = CDiagnosticPlots(mCounts, 'first Test')
fBatch.1 = lData$batch
# check for batch effects with some covariates
boxplot.median.summary(oDiag, fBatch.1, legend.pos = 'topright', axis.label.cex = 1)

plot.mean.summary(oDiag, fBatch.1, axis.label.cex = 1)

plot.sigma.summary(oDiag, fBatch.1, axis.label.cex = 1)

plot.missing.summary(oDiag, fBatch.1, axis.label.cex = 1)

plot.PCA(oDiag, fBatch.1)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag)
l
# set all parameters to false
l2 = lapply(l, function(x) x = F)

# recalculate with new parameters
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag, l2)
plot.PCA(oDiag.2, fBatch.1)

plot.dendogram(oDiag, fBatch.1, labels_cex = 1)

## check for extreme values
oDiag.2 = calculateExtremeValues(oDiag.2)
m = mGetExtremeValues(oDiag.2)
## which sample has the most extreme values
apply(m, 2, function(x) sum(x > 0))
## which variable was extreme most times
v = apply(m, 1, function(x) sum(x > 0))
v[which(v > 1)]
mCounts[v>1,]

