optimizePS = function(X, k, S, P, error=10^-50, myexp=10^-9)
{
  S0 = matrix(0,dim(S)[1],dim(S)[2])
  P0 = matrix(0,dim(P)[1],dim(P)[2])
  Sold = S0
  Pold = P0  
  
  # update P
  istep = 0 
  while((sum((P-Pold)^2)>=error || sum((S-Sold)^2)>=error) && istep <= 100000)
  {
    V = X%*%t(S)
    V0 = P%*%S%*%t(S)
    V1 = t(P)%*%X
    V2 = t(P)%*%P%*%S + myexp
    Sold = S0
    Pold = P0    
    for( i in 1:k )
    {
      for( j in 1:dim(P)[1])
      {
        if(V0[j,i] == 0)
        {
          P0[j,i] = P[j,i]     
        }else{
          P0[j,i] = P[j,i]*V[j,i]/V0[j,i]        
        }
      }
      for( z in 1:dim(S)[2])
      {
        if(V2[i,z] == 0)
        {
          S0[i,z] = S[i,z]
        }else{
          S0[i,z] = S[i,z]*V1[i,z]/V2[i,z]       
        }        
      }     
    }
    P = normalize(P0) 
    S = normalize(S0) 
    istep = istep + 1
  }   
  res = NULL 
  if( sum(apply(abs(diff(as.matrix(P))),2,sum) == 0) == 0 )
  {  
    colnames(S) = rownames(S) = colnames(P) = rownames(P) = NULL
    res = list('P'=P,'S'=S,'iterations'=istep)
  }else{
    res = list('P'=NULL,'S'=NULL,'iterations'=istep)
  }
  return(res)  
}

IterateInitialvalue = function ( X=NULL, k=NULL, S=NULL, P=NULL, error=10^-10, myexp=10^-9)
{
  if(is.null(X))
  {
    stop("Input is empty!")
  }
  if(is.null(k))
  {
    k = 1
  }  
  if(is.data.frame(X))
  {
	X = as.matrix(X)
  }
  if(is.vector(X))
  {
	X = t(t(X))
  }   
  X =  normalize(X)
  if(storage.mode(X)!="double")
  {
    storage.mode(X) <- "double"
  }  
  mn <- dim(X)
  m <- mn[1]
  n <- mn[2]
  if(is.null(S))
  {
    S = t(rdirichlet(n, rep(1,k) ))
  }
  S = normalize(S) 
  if(is.null(P))
  {
    P = X%*%ginv(S)
    ix = P < 0
    P[ix] = 0
  } 
  P = normalize(P)
  
  # optimize
  res = optimizePS(X, k, S, P, error=error, myexp=myexp)
  return(list('P'=res$P,'S'=res$S))
}
