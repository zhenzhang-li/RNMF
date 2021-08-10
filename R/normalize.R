f1 = function(x)
{
    if(sum(x)>0)
    {
      x = x/sum(x)
    }
    return(x)
}
normalize = function(X)
{
  if(is.data.frame(X))
  {
	X = as.matrix(X)
  }
  if(is.vector(X))
  {
	X = t(X)
  } 
  da = apply(X, 2, f1)
  if(is.data.frame(da))
  {
	da = as.matrix(da)
  }
  if(is.vector(da))
  {
	da = t(da)
  }    
  return(da)
}
