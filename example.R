# Name: Example.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 25/05/2017
# Desc: example usage of cdiagnostic plots



# utility function to load object
f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}

## load the test data
lData = f_LoadObject('lData.publish.rds')
# it is a list with various components
names(lData)

mCounts = lData$expression

