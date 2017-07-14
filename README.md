# CDiagnosticPlots
Draw diagnostic plots for data matrix using various statistics and groupings or batches of interest.

## Class definition
**Name: CDiagnosticPlots**  
**Desc: Perform some diagnostics for High dimensional data matrix**  
**Slots:**  
**mData = data matrix**  
**csTitle = character string with title for plots**  
**lData = list to hold various types of data and results**  
**lParam = list of parameters to change for PCA and other Plots**  
```R
setClass('CDiagnosticPlots', slots=list(mData='matrix', csTitle='character', lData='list', lParam='list'))
```  

## Constructor  
**Name: CDiagnosticPlots**  
**ARGS:**  
mData = data matrix with samples in columns (subject space) and variables in rows (variable space)  
**NOTE** It is suggested that data is logged and zeros will be considered as missing values.  
csTitle = the title for the plots, a character string  
**Desc:**  The constructor internally checks for data types, creates the object, sets the parameters by calling CDiagnosticPlotsGetParameters and builds the object using CDiagnosticPlotsBuild function.  

## Object Building  
**Name: CDiagnosticPlotsBuild**  
**ARGS:**  
obj = Object of type CDiagnosticPlots  
**DESC:**  This function is not called directly by the user, and is called from the constructor. It has some internal functions and performs some steps using these functions.  
**Internal Functions**  
```R
## Mean and SD at population level
  logPostNorm = function(theta, data)
## Log posterior for data missing and present
  logPostBin = function(theta, data)
## summary for the means and sd of each sample in matrix
  mean.sd = function(mData)
  ## extracts data frames with following information for each sample
  list('mean'=dfMS, 'sigma'=dfSD, 'tpar'=lTpar)
## summary of missing data
  missing.bin = function(mData)
## perform PCA
  getPCA = function(mData, scaleSubjects=T, scaleVariables=T)
## perform HC
  getHClust = function(mData, scaleSubjects=T, scaleVariables=T)
```  
**Object Building Steps**  
1.  Get the summary for each statistic for each sample.  
2.  Fill the lData slot with the statistics and sub objects.  

## Plotting, Accessor and Utility Functions.  
**Name: CDiagnosticPlotsGetParameters & CDiagnosticPlotsSetParameters**  
**Type:** Slot accessor function.  
  
**Name: boxplot.median.summary**  
**Type:** Plotting.  

**Name: plot.mean.summary**  
**Type:** Plotting.  

**Name: plot.sigma.summary**  
**Type:** Plotting.  

**Name: plot.missing.summary**  
**Type:** Plotting.  
**NOTE:** NA, Infinite and 0 are considered as missing values. 

**Name: plot.PCA**  
**Type:** Plotting.  

**Name: plot.dendogram**  
**Type:** Plotting.  

**Name: plot.heatmap**  
**Type:** Plotting.  

**Name: calculateExtremeValues & mGetExtremeValues**  
**Type:** Utility.  
**DESC:** This function may be useful in some cases, and will mark high data values within the calculated parameter space and model. We assume each sample is normally distributed - we check if any observed data points are larger than than what would be expected under that model, under repeated simulations. (Use the results as a guideline only)
