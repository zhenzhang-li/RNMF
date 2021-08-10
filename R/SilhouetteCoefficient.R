SilhouetteCoefficient = function(res)
{
  S = res$S
  k = res$k
  kset = unique(sort(k))
  processStabAvg = rep(1, max(kset))

  for( i in 1:length(kset) )
  {
    k0 = kset[i]
    index = which(k == k0)
    kocount = length(index)
    if(kocount>0)
    {
      data = array(NA,dim=c(kocount,dim(S[[1]])[2],k0))
      for(j in 1:k0)
      {
        data[1,,j] = S[[index[1]]][j,]
        if(kocount>1)
		{
			for(m in 2:kocount)
			{         
			  corset = rep(0,k0)
			  for(ki in 1:k0)
			  {
				corset[ki] = sum((data[1,,j]-S[[index[m]]][ki,])^2)
			  }
			  indexcor = which.min(corset) 
			  data[m,,j] = S[[index[m]]][indexcor,]
			}  				
		}  		      
      }
      
      Sseq = NULL
      for(j in 1:k0)
      {
        ai = NULL
        bi = NULL
        for(n in 1:kocount)
        {
          sami = data[n,,j]
          aiset = NULL
          for(m in 1:kocount)
          {           
            if(n!=m)
            {
              aiset = c(aiset, as.numeric(dist(rbind(sami, data[m,,j]))))
            }else{
			  aiset = 0
			}
          }
          ai = c(ai, mean(aiset))
          
          biset = NULL          
          if(k0==1)
          {
            biset = c(biset, 1)
          }else{
            for(ij in 1:k0)
            {
              if(j!=ij)
              {
                biset2 = NULL
                for(m in 1:kocount)
                {           
                  biset2 = c(biset2, as.numeric(dist(rbind(sami, data[m,,ij]))))
                }
                biset = c(biset, mean(biset2))
              }
            } 
          }
          bi = c(bi, min(biset))
        }
        for(ij in 1:length(ai))
        {
          if(ai[ij]<bi[ij])
          {
            Sseq = c(Sseq, 1 - ai[ij]/bi[ij])
          }
          if(ai[ij]==bi[ij])
          {
            Sseq = c(Sseq, 0)
          }
          if(ai[ij]>bi[ij])
          {
            Sseq = c(Sseq, bi[ij]/ai[ij] - 1)
          }
        }
      }     
      processStabAvg[i] = mean(Sseq)
    }
  }
  return(processStabAvg)
}
