##   P
#Subtype	SBS1	SBS2
#A[C>A]A	0.0322862036350085	0.000563894105147364
#A[C>A]C	0.013632971643522	0.00268369719318389

## 	S
#SampleID	FP1705100059DN01	FP1705100061DN01
#SBS1	0.564108835202935	0.156024536682433


samMatrixSigs <- function(Pfile, Sfile, AnalCOSMICSigType = 'SBS')
{
  library( data.table )
  P = as.data.frame(fread( Pfile ))
  S = as.data.frame(fread( Sfile ))
  m = dim(P)[1]
  k = dim(P)[2]
  n = dim(S)[2]
  dir.create('samplesResults', showWarnings = FALSE, recursive = TRUE, mode = "0777")
  samMatrixsigsData = NULL
  try(load(paste("samplesResults/", AnalCOSMICSigType, ".samMatrixsigsData.Rdata", sep="")), silent=TRUE)
  if(is.null(samMatrixsigsData))
  {  
	fud = function(da)
	{
	  f1 = function(x)
	  {
	    x = as.numeric(x)
		if(sum(x)>0)
		{
		  x = x/sum(x)
		}
		return(x)
	  }
	  da[,2:ncol(da)] = apply(da[,-1], 2, f1)
	  return(da)
	}  
    SampleIDs = colnames(S)[-1]
    Stmp = S[, -1]
    nl = n - 1
    for( i in 1:nl )
    {
      dir.create(paste('samplesResults/', SampleIDs[i], "/", AnalCOSMICSigType, sep=""), showWarnings = FALSE, recursive = TRUE, mode = "0777")
      samMatrixsigsData$SampleID[[i]] = SampleIDs[i]
      Pptmp = NULL
      for( j in 1:m )
      {
        x = P[j, -1]*Stmp[, i]
		#asum = sum(as.numeric(x))
		#if(asum == 0)
		#{
		#  asum = 1
		#}
        #Pptmp = rbind(Pptmp, c(P[j, 1], as.numeric(x)/asum))
		Pptmp = rbind(Pptmp, c(P[j, 1], as.numeric(x)))
      }
      colnames(Pptmp) = colnames(P)
	  Pptmp = fud(Pptmp)	  
      samMatrixsigsData$P[[i]] = Pptmp
      write.table(Pptmp, file=paste('samplesResults/', SampleIDs[i], "/", AnalCOSMICSigType, "/", SampleIDs[i], ".P.txt",sep=""), quote=F, col.names=T, row.names=F, se="\t")
    }
    save(samMatrixsigsData, file=paste("samplesResults/", AnalCOSMICSigType, ".samMatrixsigsData.Rdata", sep=""))
  }
  
  return(samMatrixsigsData)
}
