# Name: CDiagnosticPlots.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 25/05/2017
# Desc: class to create a CDiagnosticPlots Object


library(methods)
library(LearnBayes)
library(car)
logit.inv = function(p) {exp(p)/(exp(p)+1) }

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



