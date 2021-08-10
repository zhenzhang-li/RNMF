fd = function (x1, x2, nsim = 10000) 
{
    n1 <- length(x1)
    n2 <- length(x2)
    n <- n1 + n2
    x <- c(x1, x2)
    dbar <- mean(x2) - mean(x1)
    z <- array(, nsim)
    for (i in 1:nsim) 
    {
      mn <- sample(n, n2, replace = FALSE)
      dbardash <- mean(x[mn]) - mean(x[-mn])
      z[i] <- dbardash
    }
    pval <- (sum(z >= abs(dbar)) + sum(z <= -abs(dbar)))/nsim
    return(pval)
}

fu = function(i, G, f, sample.id, Hugo, SigConName, Mutfreq)
{
  sam = NULL
  keep = 'N'
  ix = f[, c(Hugo)] == G[i]
  sam = unique(as.character(f[ix, c(sample.id)]))
  # remove those gene with lower than mutation frequecy we difine!
  if(length(sam) > floor(length(SigConName)*Mutfreq))
  {   
    keep = 'Y'
  }
  return(keep)
}

genePerMutSigs <- function(Sfile = NULL, choose.Sigs = 'SBS1', Gfile = NULL, file = NULL, sample.id = 'Tumor_Sample_Barcode', Hugo = 'Hugo_Symbol', Mutfreq = 0.04, threshold = 0.06, Qvalue = 0.1, plot = TRUE) 
{
  if(is.null(Sfile))
  {
    stop("Please enter a abundance fractions matrix S with a format like the output result of RNMF software!\n")
  } 
  
  if(is.null(Gfile))
  {
    stop("Please enter a abundance fractions matrix of gene for all samples with a format like the output result of RNMF software!\n")
  }

  if(is.null(file))
  {
    stop("Please enter a non-silent mutation dataset, such as the annotation result file in MAF format!\n")
  }

  library(data.table)
  S = as.data.frame(fread(Sfile))
  G = as.data.frame(fread(Gfile))
  f = as.data.frame(fread(file)) 
  dir.create('samplesResults/genePerMutSigs/', showWarnings = FALSE, recursive = TRUE, mode = "0777")
  f = f[, c(sample.id, Hugo)]
  fss = unique(f[, sample.id])
  indexLL = NULL
  indexLL1 = NULL  
  for( i in 1:length(fss) )
  {
	ix = which(colnames(S) == fss[i])
	iy = which(colnames(G) == fss[i]) 
	if(sum(ix)>0)
	{
		indexLL = c(indexLL, ix)	
	}
	if(sum(iy)>0)
	{
		indexLL1 = c(indexLL1, iy)	
	}	
  }
  indexLL = c(1, indexLL)
  S = S[, indexLL]  
  indexLL1 = c(1, indexLL1)
  G = G[, indexLL1]
  index = which(S[, 1] == choose.Sigs)
  S = S[index, -1]
  SigConName = colnames(S)  
  ng = nrow(G)  
  tmp = matrix(1:ng, ng, 1)
  tmpYN = apply(tmp, 1, fu, G$GeneName, f, sample.id, Hugo, SigConName, Mutfreq)
  ix = tmpYN == 'Y'
  G = G[ix,]
  GGNAME = colnames(G)[-1]	  
  ng = nrow(G)
  Gmatrix = matrix(0, ng, 3)
  Gmatrix[, 2] = 1
  Gmatrix[, 3] = 1  
  rownames(Gmatrix) = G[, 1]
  colnames(Gmatrix) = c('Median','Pvalue','Qvalue')
  for(i in 1:ng)
  {
    gname = G[i, 1]
    ix = f[, c(Hugo)] == gname
    if(sum(ix)>0)
    {
	  gsamN = f[ix, c(sample.id)]  
	  index_s1 = NULL  
	  index_med = NULL  	  
      for(j in 1:length(gsamN))
      {
        iy = which(SigConName == gsamN[j])
		iy1 = which(GGNAME == gsamN[j])
		index_s1 = c(index_s1, iy)
		index_med = c(index_med, iy1)		
      }
	  index_s2 = setdiff(1:length(SigConName), index_s1)
      s1 = as.numeric(as.character(S[1, index_s1]))	  
	  s2 = as.numeric(as.character(S[1, index_s2]))
      Gmatrix[i, 2] = fd(s1, s2, nsim = 10000)  
	  ix = G[, 1] == gname
	  Nnum = as.numeric(as.character(G[ix, -1]))
	  Gmatrix[i, 1] = median(Nnum[index_med])	  
    }
  }
  Gmatrix[, 3] = p.adjust(as.numeric(as.character(Gmatrix[, 2])),method="fdr",n=length(as.numeric(as.character(Gmatrix[, 2]))))
  Gmatrix = cbind(t(t(rownames(Gmatrix))), Gmatrix)
  colnames(Gmatrix)[1] = 'GeneName'
  write.table(Gmatrix, file=paste('samplesResults/genePerMutSigs/Gene.final.mutation_signature.',choose.Sigs,'.xls',sep=""),quote=F,col.names=T,row.names=F,sep="\t")
    
  if(plot)
  {
    file2 = paste('samplesResults/genePerMutSigs/Gene.final.mutation_signature.',choose.Sigs,'.xls',sep="")
    da = read.table(file2,header=T,sep="\t")
	eps = .Machine$double.eps
	ix = da$Pvalue < eps
	da$Pvalue[ix] = eps
    x = da$Median
    y = -log10(da$Pvalue)
    redcol<-colorRampPalette(c("grey80", "black"))(100)
    pdf(paste('samplesResults/genePerMutSigs/Gene.final.mutation_signature.',choose.Sigs,'.pdf',sep=""),6,6)
    plot(0,0,type='l',xlim=c(0, 1),ylim=c(0,ceiling(max(y)+1)),xlab='Median contribution abundance of gene with non-silent mutations',ylab="-log10( PerMutation P value )")
    n=length(y)
    col=rep(NA,n)
    for(i in 1:n){
      col[i]=redcol[ceiling(da$Pvalue[i]*100)]
    }
    ix = da$Qvalue < Qvalue & da$Median >= threshold
    col[ix]='red'
    points(x,y,pch=20,col=col,cex=2^(1-da$Pvalue))
    if(sum(ix)>0)
    {
      text(x[ix]+0.08,y[ix]+0.05,labels=as.character(da$GeneName[ix]),col='blue',cex=1,font=3)
    }
    dev.off()  
  }

}
