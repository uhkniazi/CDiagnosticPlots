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
plot.mean.summary(oDiag, fBatch.1)

plot.sigma.summary(oDiag, fBatch.1)

plot.missing.summary(oDiag, fBatch.1)

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