corf = function(A, B)
{
  m = dim(A)[2]
  mname = colnames(A)
  n = dim(B)[2]
  nname = colnames(B)
  corM = matrix(0, m, n)  
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      corM[i,j] = as.numeric(as.character(A[,i]%*%B[,j]/(sqrt(A[,i]%*%A[,i])*sqrt(B[,j]%*%B[,j])))) 
    }
  }
  rownames(corM) = mname
  colnames(corM) = nname
  return(corM)
}