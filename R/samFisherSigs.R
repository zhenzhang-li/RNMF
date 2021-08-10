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

samFisherSigs <- function(Sfile = NULL, choose.Sigs = 'SBS1', file = NULL, sample.id = 'Tumor_Sample_Barcode', Hugo = 'Hugo_Symbol', Mutfreq = 0.04, threshold = 0.06, Qvalue = 0.1, plot = TRUE) 
{
  if(is.null(Sfile))
  {
    stop("Please enter a abundance fractions matrix S with a format like the output result of RNMF software!\n")
  }
  
  if(is.null(file))
  {
    stop("Please enter a non-silent mutation dataset, such as the annotation result file in MAF format!\n")
  }
  
  library(data.table)
  S = as.data.frame(fread(Sfile))
  f = as.data.frame(fread(file)) 
  dir.create('samplesResults/samFisherSigs/', showWarnings = FALSE, recursive = TRUE, mode = "0777")
  f = f[, c(sample.id, Hugo)]
  fss = unique(f[, sample.id])
  indexLL = NULL
  for( i in 1:length(fss) )
  {
	ix = which(colnames(S) == fss[i])
	if(sum(ix)>0)
	{
		indexLL = c(indexLL, ix)	
	}
  }
  indexLL = c(1, indexLL)
  S = S[, indexLL]  
  index = which(S[, 1] == choose.Sigs)
  S = S[index, -1]
  SigConName = colnames(S)
  G = as.character(unique(f[, c(Hugo)])) 
  ng = length(G)
  tmp = matrix(1:ng, ng, 1)
  tmpYN = apply(tmp, 1, fu, G, f, sample.id, Hugo, SigConName, Mutfreq)
  ix = tmpYN == 'Y'
  G = G[ix]
  G = as.character(G)
  ng = length(G)
  Gmatrix = matrix(0, ng, 3)
  Gmatrix[, 2] = 1
  Gmatrix[, 3] = 1  
  rownames(Gmatrix) = G
  colnames(Gmatrix) = c('Median','Pvalue','Qvalue')
  for(i in 1:ng)
  {
    gname = G[i]
    ix = f[, c(Hugo)] == gname
    if(sum(ix)>0)
    {
      gsamN = f[ix, c(sample.id)]
      Highcon_Mutated = 0
      Highcon_Wildtype = 0
      Lowcon_Mutated = 0
      Lowhcon_Wildtype = 0
      index_s1 = NULL  
      for(j in 1:length(gsamN))
      {
        iy = which(SigConName == gsamN[j])
        index_s1 = c(index_s1, iy)
      }
      index_s2 = setdiff(1:length(SigConName), index_s1)
      s1 = as.numeric(as.character(S[1, index_s1]))	  
      s2 = as.numeric(as.character(S[1, index_s2]))
      ixs1 = s1 >= threshold 
      Highcon_Mutated = sum(ixs1)
      Lowcon_Mutated = length(s1) - Highcon_Mutated
      ixs2 = s2 >= threshold 
      Highcon_Wildtype = sum(ixs2)
      Lowhcon_Wildtype = length(s2) - Highcon_Wildtype
      Gmatrix[i, 2] = fisher.test(matrix(c(Highcon_Mutated,Highcon_Wildtype,Lowcon_Mutated,Lowhcon_Wildtype),2,2),alternative = 'greater')$p.value           
      Gmatrix[i, 1] = median(s1)
    }
  }
  Gmatrix[, 3] = p.adjust(as.numeric(as.character(Gmatrix[, 2])),method="fdr",n=length(as.numeric(as.character(Gmatrix[, 2]))))
  Gmatrix = cbind(t(t(rownames(Gmatrix))), Gmatrix)
  colnames(Gmatrix)[1] = 'GeneName'
  write.table(Gmatrix, file=paste('samplesResults/samFisherSigs/Sam.final.mutation_signature.',choose.Sigs,'.xls',sep=""),quote=F,col.names=T,row.names=F,sep="\t")
  
  if(plot)
  {
    file2 = paste('samplesResults/samFisherSigs/Sam.final.mutation_signature.',choose.Sigs,'.xls',sep="")
    da = read.table(file2,header=T,sep="\t")
    eps = .Machine$double.eps
    ix = da$Pvalue < eps
    da$Pvalue[ix] = eps
    x = da$Median
    y = -log10(da$Pvalue)
    redcol<-colorRampPalette(c("grey80", "black"))(100)
    pdf(paste('samplesResults/samFisherSigs/Sam.final.mutation_signature.',choose.Sigs,'.pdf',sep=""),6,6)
    plot(0,0,type='l',xlim=c(0,1),ylim=c(0,ceiling(max(y)+1)),xlab='Median mutational exprosure of samples with non-silent mutations',ylab="-log10( Fisher's exact P value )")
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
