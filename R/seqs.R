seqs = function(x)
{
  x = strsplit(x,'')[[1]]
  xnew = NULL
  for( i in 1:length(x) )
  {
    if(x[i] %in% 'C')
    {
      ach = 'G'
    }
    if(x[i] %in% 'G')
    {
      ach = 'C'
    }
    if(x[i] %in% 'A')
    {
      ach = 'T'
    }
    if(x[i] %in% 'T')
    {
      ach = 'A'
    }
    if(x[i] %in% 'N')
    {
      ach = 'N'
    }      
    xnew = c(xnew, ach)
  }
  fseq = NULL
  for( i in 1:length(xnew) )
  {
	fseq = paste(fseq,xnew[i],sep="")
  }
  return(reverse(fseq))
}
