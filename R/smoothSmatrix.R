smoothSmatrix <- function(v, p, eps)
{
  n = nrow(v)
  m = ncol(v) 
  if(is.data.frame(v))
  {
	v = as.matrix(v)
  }
  if(is.vector(v))
  {
	v = t(v)
  } 
  if(is.data.frame(p))
  {
	p = as.matrix(p)
  }
  if(is.vector(p))
  {
	p = t(t(p))
  }   
  r = ncol(p)
  w = normalize(p)
  w1 = w
  h = ginv(w)%*%v
  h[h<eps] = eps   
  h = normalize(h)    
  niter = 10
  for( i in 1:niter )
  {
    x1 = matrix(rep(apply(w, 2, sum), m), r, m)
    atmp = w%*%h
    atmp[atmp<eps] = eps 
    h = h*(t(w)%*%(v/atmp))/x1           
  }   
  h[h<eps] = 0  
  h = normalize(h)  
  return(h)
}