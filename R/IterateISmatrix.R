IterateISmatrix = function(X, S, P, error=10^-10, myexp=10^-9, eps = 0)
{
  if(is.data.frame(X))
  {
	X = as.matrix(X)
  }
  if(is.vector(X))
  {
	X = t(t(X))
  } 
  if(is.data.frame(P))
  {
	P = as.matrix(P)
  }
  if(is.vector(P))
  {
	P = t(t(P))
  } 
  P = normalize(P)
  if(is.data.frame(S))
  {
	S = as.matrix(S)
  }  
  if(is.vector(S))
  {
    S = t(S)
  }  
  S = normalize(S)
  SNEWtmp = S
  SNEW = matrix(0, dim(S)[1], dim(S)[2])
  istep = 0 
  while( sum((SNEW-S)^2) >= error && istep <= 100000000 )
  { 
    S = SNEWtmp      
    V1 = t(P)%*%X
    V2 = t(P)%*%P%*%S + myexp  
    V2[V2 == 0] = eps
    SNEW = S*V1/V2
    SNEW = normalize(SNEW)	
	SNEWtmp = SNEW
    istep = istep + 1
  } 
  cat("Iterative error:",sum((S-SNEW)^2)," &  Iteration:",istep,"\n")
  return(S)  
}
