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

## load the test data 1 and 2, same data normalized in different ways
lData.1 = f_LoadObject('lData.test_data.rds')
lData.2 = f_LoadObject('lData.test_data_second.rds')

# it is a list with various components
names(lData.1)
names(lData.2)

## lets check the 2 normalization methods and see how they compare
## create the object
oDiag.1 = CDiagnosticPlots(lData.1$data, 'method 1')
oDiag.2 = CDiagnosticPlots(lData.2$data, 'method 2')

# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
# e.g. in this case it is different lanes/machines
fBatch = factor(lData.1$batch)

par(mfrow=c(1,2))

## compare the 2 methods using various plots
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.PCA(oDiag.2, fBatch, cex.main=1)

plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## how about the extreme values
oDiag.1 = calculateExtremeValues(oDiag.1)
oDiag.2 = calculateExtremeValues(oDiag.2)
m1 = mGetExtremeValues(oDiag.1)
m2 = mGetExtremeValues(oDiag.2)

## samples with most extreme values
apply(m1, 2, function(x) sum(x > 0))
apply(m2, 2, function(x) sum(x > 0))

## variables that are contributing to this
v1 = apply(m1, 1, function(x) sum(x > 0))
v2 = apply(m2, 1, function(x) sum(x > 0))

which(v1 > 0)
which(v2 > 0)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

## this should give an error as scaling can't be done
## if all the vector 0 for PCA
oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)
## reset flag for jittering
l$PCA.jitter = T
## should work this time
oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)
plot.PCA(oDiag.1, fBatch)
plot.PCA(oDiag.1.2, fBatch)
plot.dendogram(oDiag.1, fBatch)
plot.dendogram(oDiag.1.2, fBatch)





