confint<-function(x,alpha=0.05)
{
  library( Rmisc )
  a = data.frame(t(CI(x, ci=alpha)[c(2,1,3)]))
  return(a)
}

confidencexy = function(x,level=0.95)
{
  x.level = confint(x,level)  
  return(x.level)
}