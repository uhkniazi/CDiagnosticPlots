# Name: CDiagnosticPlots.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 25/05/2017
# Desc: class to create a CDiagnosticPlots Object


library(methods)
library(LearnBayes)

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

## Mean and SD with group and population levels
## Mean and SD at population level
logPostNormGroup = function(theta, data){
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

getSummary = function(mData, type=mean.sd){
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
  
  lSummary = mean.sd(mData)
}



