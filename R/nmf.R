nmf <- function(v, r, eps, alpha, beta)
{
  if(is.data.frame(v))
  {
	v = as.matrix(v)
  }
  if(is.vector(v))
  {
	v = t(t(v))
  }
  n = nrow(v)
  m = ncol(v) 
  w = t(rdirichlet(r, rep(1, n))) 
  h = t(rdirichlet(m, rep(1, r))) 
  if(is.data.frame(w))
  {
	w = as.matrix(w)
  }
  if(is.vector(w))
  {
	w = t(t(w))
  } 
  if(is.data.frame(h))
  {
	h = as.matrix(h)
  }  
  if(is.vector(h))
  {
    h = t(h)
  }      
  xin = c(as.vector(w), as.vector(h))
  res = nlm(fobjective, xin, r, n, m, v, beta, alpha, gradtol = 1e-4, stepmax = 1000, steptol = 1e-4, iterlim = 100) 
  xin = abs(res$estimate)
  xin[xin<eps] = eps
  w = matrix(xin[1:(n*r)], n, r)  
  h = matrix(xin[-c(1:(n*r))], r, m)  
  w = normalize(w) 
  h = normalize(h)  
  niter = 20000  
  for( i in 1:niter )
  {
    x1 = matrix(rep(apply(w, 2, sum), m), r, m)
	atmp = w%*%h
	atmp[atmp<eps] = eps
    h = h*(t(w)%*%(v/atmp))/x1      
    x2 = matrix(rep(apply(h, 1, sum), n), n, r, byrow = TRUE)
	atmp = w%*%h
	atmp[atmp<eps] = eps	
    w = w*((v/atmp)%*%t(h))/x2        
  } 
  w[w<eps] = eps
  h[h<eps] = eps  
  w = normalize(w) 
  h = normalize(h)   
  return(list('P'=w,'S'=h))
}