# Name: CDiagnosticPlots.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 25/05/2017
# Desc: class to create a CDiagnosticPlots Object


library(methods)
library(LearnBayes)
library(car)
logit.inv = function(p) {exp(p)/(exp(p)+1) }


###### Class CDiagnosticPlots
# Name: CDiagnosticPlots
# Desc: Perform some diagnostics for High dimensional data matrix 
# Slots:
# mData = data matrix
# csTitle = character string with title for plots
# lData = list to hold various types of data and results
# lParam = list of parameters to change for PCA and other Plots
setClass('CDiagnosticPlots', slots=list(mData='matrix', csTitle='character', lData='list', lParam='list'))


## object construction function i.e. steps performed to build object, called via constructor
setGeneric('CDiagnosticPlotsBuild', function(obj, ...)standardGeneric('CDiagnosticPlotsBuild'))
setMethod('CDiagnosticPlotsBuild', signature = 'CDiagnosticPlots', definition = function(obj, ...){
  #### define internal functions
  ## calculate summaries
  
  ## Mean and SD at population level
  logPostNorm = function(theta, data){
    # we define the sigma on a log scale as optimizers work better
    # if scale parameters are well behaved
    s = exp(theta[2])
    m = theta[1]
    d = data$vector # data 
    log.lik = sum(dnorm(d, m, s, log=T))
    log.prior = 1
    log.post = log.lik + log.prior
    return(log.post)
  }
  
  ## Log posterior for data missing and present
  logPostBin = function(theta, data){
    ## theta contains parameters we wish to track
    th = logit.inv(theta['theta'])
    suc = data$success
    fail = data$fail
    # define likelihood function
    lf = function(s, f, t) return(dbinom(s, s+f, t, log=T))
    # calculate log posterior
    log.prior = dbeta(th, 10, 10, log = T) # 
    val = lf(suc, fail, th) + log.prior
    return(val)
  }
  
  ## summary for the means and sd
  mean.sd = function(mData){
    lRet = apply(mData, 2, function(inData) {
      # choose a starting value
      start = c('mu'=mean(inData), 'sigma'=log(sd(inData)))
      lData = list('vector'=inData)
      ## try the laplace function from LearnBayes
      fit = laplace(logPostNorm, start, lData)
      return(fit)
    })
    ## get sd and average for the parameters
    getsds = function(f){
      m = f$mode['sigma']
      se = sqrt(diag(f$var))['sigma']
      m.up = exp(m+1.96*se)
      m.down = exp(m-1.96*se)
      m = exp(m)
      ret= c(m, m.up, m.down)
      names(ret) = c('m', 'm.up', 'm.down')
      return(ret)
    }
    getms = function(f){
      m = f$mode['mu']
      se = sqrt(diag(f$var))['mu']
      m.up = m+1.96*se
      m.down = m-1.96*se
      ret= c(m, m.up, m.down)
      names(ret) = c('m', 'm.up', 'm.down')
      return(ret)
    }
    dfSD = t(sapply(lRet, getsds))
    dfMS = t(sapply(lRet, getms))
    return(list('mean'=dfMS, 'sigma'=dfSD))
  }
  
  ## summary of missing data
  missing.bin = function(mData){
    lRet = apply(mData, 2, function(inData) {
      inData[is.na(inData) | !is.finite(inData)] = 0
      inData = as.logical(inData)
      # choose a starting value
      start = c('theta'=logit(0.5))
      lData = list('success'=sum(inData), fail=sum(!inData))
      ## try the laplace function from LearnBayes
      fit = laplace(logPostBin, start, lData)
      return(fit)
    })
    ## get sd and average for the parameters
    getms = function(f){
      m = f$mode['theta']
      se = sqrt(diag(f$var))['theta']
      m.up = m+1.96*se
      m.down = m-1.96*se
      ret = c(logit.inv(m), logit.inv(m.up), logit.inv(m.down))
      names(ret) = c('psuc', 'psuc.up', 'psuc.down')
      return(ret)
    }
    dfMS = t(sapply(lRet, getms))
    return(list('present'=dfMS))
  }
  
  getSummary = function(mData, type){
    lSummary = type(mData)
  }
  
  ## perform PCA
  getPCA = function(mData, scaleSubjects=T, scaleVariables=T){
    ## add a normal jitter to each cell to remove zeros
    n = dim(mData)[1] * dim(mData)[2]
    mData = mData + rnorm(n)
    
    #### standardize samples first
    if (scaleSubjects){
      s = apply(mData, 2, sd)
      mData = sweep(mData, 2, s, '/')
    }
    mData = t(mData)
    # set scaling to TRUE to scale variables i.e. genes/variables in columns
    pr.out = prcomp(mData, scale = scaleVariables)
    return(pr.out)
  }
  
  getHClust = function(mData, scaleSubjects=T, scaleVariables=T){
    ## add a normal jitter to each cell to remove zeros
    n = dim(mData)[1] * dim(mData)[2]
    mData = mData + rnorm(n)
    
    #### standardize samples first
    if (scaleSubjects){
      s = apply(mData, 2, sd)
      mData = sweep(mData, 2, s, '/')
    }
    if (scaleVariables){
      s = apply(mData, 1, sd)
      mData = sweep(mData, 1, s, '/')
    }
    hc = hclust(dist(t(mData)))
    return(hc)
  }
  
  ########## end internal functions
  
  ## steps to build the object
  # step 1, calculate various summaries for one dimensional statistics
  lData = getSummary(obj@mData, mean.sd)
  lData = append(lData, getSummary(obj@mData, missing.bin))
  l = getPCA(obj@mData, obj@lParam$PCA.scaleSubjects, obj@lParam$PCA.scaleVariables)
  lData$PCA = l
  l = getHClust(obj@mData, obj@lParam$HC.scaleSubjects, obj@lParam$HC.scaleVaribles)
  lData$HC = l
  obj@lData=lData
  return(obj)
})


## constructor
## ARGS: mData = data matrix with samples in columns (subject space) and variables in rows (variable space)
CDiagnosticPlots = function(mData, csTitle){
  # check some diagnostics here and create object
  if (class(mData) != 'matrix') stop('CDiagnosticPlots: mData is not object of class matrix')
  
  ## create the object
  ob = new('CDiagnosticPlots', mData=mData, csTitle=csTitle, lData=list(), lParam=list())
  ob@lParam = CDiagnosticPlotsGetParameters(ob)
  return(CDiagnosticPlotsBuild(ob))
}

## setting and getting parameters
setGeneric('CDiagnosticPlotsGetParameters', function(obj, ...)standardGeneric('CDiagnosticPlotsGetParameters'))
setMethod('CDiagnosticPlotsGetParameters', signature = 'CDiagnosticPlots', definition = function(obj, ...){
  # parameters empty set first time
  if(length(obj@lParam) == 0) {
    obj@lParam$PCA.scaleSubjects = T
    obj@lParam$PCA.scaleVariables = T
    obj@lParam$HC.scaleSubjects = T
    obj@lParam$HC.scaleVaribles = T
  }
  return(obj@lParam)
})

setGeneric('CDiagnosticPlotsSetParameters', function(obj, lParam)standardGeneric('CDiagnosticPlotsSetParameters'))
setMethod('CDiagnosticPlotsSetParameters', signature = 'CDiagnosticPlots', definition = function(obj, lParam){
  # check if the names match
  if (!all(names(CDiagnosticPlotsGetParameters(obj)) %in% names(lParam))) {
    warning('CDiagnosticPlotsSetParameters: the names of parameters do not match')
    return(obj)
  }
  # or else put new parameters
  obj@lParam = lParam
  return(obj)
})

###################################### analysis related functions
## plotting functions

setGeneric('plot.mean.summary', function(obj, fBatch, legend.pos='bottomright', ...)standardGeneric('plot.mean.summary'))
setMethod('plot.mean.summary', signature = 'CDiagnosticPlots', definition = function(obj, fBatch, legend.pos='bottomright', ...){
  df = obj@lData$mean
  # order the data according to batches
  i = order(fBatch)
  col.p = rainbow(nlevels(fBatch))
  col = col.p[as.numeric(fBatch)[i]]
  plot(df[i,'m'], pch=20, col=col, ylim=c(min(df[,'m.down']), max(df[,'m.up'])), main=paste('Mean Summary', obj@csTitle),
       xlab='Samples', xaxt='n', ylab='Average')
  # order the data matrix according to batches
  df = df[i,]
  for(l in 1:nrow(df)){
    lines(x=c(l, l), y=df[l,c(2,3)], col=col[l], lwd=0.5)
  }
  legend(legend.pos, legend = levels(fBatch), fill=col.p, ncol=nlevels(fBatch))
})

setGeneric('plot.sigma.summary', function(obj, fBatch, legend.pos='bottomright', ...)standardGeneric('plot.sigma.summary'))
setMethod('plot.sigma.summary', signature = 'CDiagnosticPlots', definition = function(obj, fBatch, legend.pos='bottomright', ...){
  df = obj@lData$sigma
  # order the data according to batches
  i = order(fBatch)
  col.p = rainbow(nlevels(fBatch))
  col = col.p[as.numeric(fBatch)[i]]
  plot(df[i,'m'], pch=20, col=col, ylim=c(min(df[,'m.down']), max(df[,'m.up'])), main=paste('Sigma Summary', obj@csTitle),
       xlab='Samples', xaxt='n', ylab='Average')
  # order the data matrix according to batches
  df = df[i,]
  for(l in 1:nrow(df)){
    lines(x=c(l, l), y=df[l,c(2,3)], col=col[l], lwd=0.5)
  }
  legend(legend.pos, legend = levels(fBatch), fill=col.p, ncol=nlevels(fBatch))
})

setGeneric('plot.missing.summary', function(obj, fBatch, legend.pos='bottomright', ...)standardGeneric('plot.missing.summary'))
setMethod('plot.missing.summary', signature = 'CDiagnosticPlots', definition = function(obj, fBatch, legend.pos='bottomright', ...){
  df = obj@lData$present
  # order the data according to batches
  i = order(fBatch)
  col.p = rainbow(nlevels(fBatch))
  col = col.p[as.numeric(fBatch)[i]]
  plot(df[i,'psuc'], pch=20, col=col, ylim=c(min(df[,'psuc.down']), max(df[,'psuc.up'])), main=paste('Proportion of Data Present Summary', 
                                                                                                     obj@csTitle),
       xlab='Samples', xaxt='n', ylab='Average Proportion')
  # order the data matrix according to batches
  df = df[i,]
  for(l in 1:nrow(df)){
    lines(x=c(l, l), y=df[l,c(2,3)], col=col[l], lwd=0.5)
  }
  legend(legend.pos, legend = levels(fBatch), fill=col.p, ncol=nlevels(fBatch))
})


setGeneric('plot.PCA', function(obj, fBatch, legend.pos='bottomright', csLabels=NULL, ...)standardGeneric('plot.PCA'))
setMethod('plot.PCA', signature = 'CDiagnosticPlots', definition = function(obj, fBatch, legend.pos='bottomright', csLabels=NULL, ...){
  pr.out = obj@lData$PCA
  col.p = rainbow(nlevels(fBatch))
  col = col.p[as.numeric(fBatch)]
  plot(pr.out$x[,1:2], col=col, pch=20, xlab='Z1', ylab='Z2',
       main=paste('PCA comp 1 and 2', obj@csTitle))
  if (is.null(csLabels)) csLabels = colnames(obj@mData)
  text(pr.out$x[,1:2], labels = csLabels, pos = 1, cex=0.6)
  legend(legend.pos, legend = levels(fBatch), fill=col.p, ncol=nlevels(fBatch), cex=0.8)
})

setGeneric('plot.dendogram', function(obj, fBatch, legend.pos='topright', csLabels=NULL, ...)standardGeneric('plot.dendogram'))
setMethod('plot.dendogram', signature = 'CDiagnosticPlots', definition = function(obj, fBatch, legend.pos='topright', csLabels=NULL, ...){
  hc = obj@lData$HC
  col.p = rainbow(nlevels(fBatch))
  ## plotting a hc object with colours is not straightforward so following example from 
  ## http://stackoverflow.com/questions/18802519/label-and-color-leaf-dendrogram-in-r
  if(!require(dendextend)) stop('CDiagnosticPlots: plot.dendogram, need library dendextend for coloured dendogram')
  dend = as.dendrogram(hc)
  # Assigning the labels of dendrogram object with new colors:
  labels_colors(dend) = col.p[as.numeric(fBatch)][order.dendrogram(dend)]
  labels_cex(dend) = 0.25
  # Plotting the new dendrogram
  plot(dend, main=paste('Hierarchical clustering of distance matrix for', obj@csTitle), xlab='', sub='Coloured on Batch')
  legend(legend.pos, legend = levels(fBatch), fill=col.p, ncol=nlevels(fBatch), cex=0.8)
})







