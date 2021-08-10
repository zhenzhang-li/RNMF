S0choose = function( X, res, tail = 5 )
{
  P = res$P
  S = res$S
  k = res$k
  kset = sort(unique(k))
  KS = NULL 
  for(k0 in kset)
  {
    ix = which( k == k0 )   
    if(length(ix)>0)
    {
      ctmp = NULL
      for(j in 1:length(ix))
      {
        ctmp = c(ctmp, sum((X-P[[ix[j]]]%*%S[[ix[j]]])^2))		
      }
      tail = ifelse(tail<=length(ix),tail,length(ix))
      index = sort.list(ctmp)[1:tail]
      KS = c(KS, ix[index])
    }      
  }
  llist = list()
  for(i in 1:length(KS))
  {
    llist$P[[i]] = res$P[[KS[i]]]
    llist$S[[i]] = res$S[[KS[i]]]
    llist$E[i] = res$E[KS[i]]
    llist$k[i] = res$k[KS[i]]
  }	
  return(llist)
}