InverseNMF = function(originalGenomes, sampleNames, subtypes, AnalCOSMICSigType = c('SBS','DBS','ID'), genome.build = c("Ch37","Ch38"), cutoff = 0.06, SBS.version = c("V2","V3"), Pmatrix = NULL, steptol = 10^-10, plot = TRUE)
{
  library( data.table )
  library( MCMCpack )
  library( reshape )
  library( ggplot2 )
  originalGenomes = as.data.frame(fread(originalGenomes, header=F))
  sampleNames = as.data.frame(fread(sampleNames, header=F))
  subtypes = as.data.frame(fread(subtypes, header=F))
  dir.create('samplesResults', showWarnings = FALSE, recursive = TRUE, mode = "0777")
  
  if(length(genome.build)>1)
  {
    genome.build = "Ch37"
  }
  if(length(SBS.version)>1)
  {
    SBS.version = "V3"
  }  
  if(length(AnalCOSMICSigType)>1)
  {
    AnalCOSMICSigType = "SBS"
  } 
  n = nrow(subtypes)  
  if(is.null(Pmatrix))
  {
    if(AnalCOSMICSigType == 'ID')
    {
      data(COSMIC.ID.GRCh37)
      COSMIC = COSMIC[1:n, ]
    }  
    if(AnalCOSMICSigType == 'DBS')
    {
	  data(list=paste("COSMIC.", AnalCOSMICSigType, ".GR", genome.build, sep = ""))	  
    }	
    if(AnalCOSMICSigType == 'SBS')
    {
      if(SBS.version == 'V2')
      {
        data("COSMIC.SBS.V2")
      }else{
		data(list=paste("COSMIC.", AnalCOSMICSigType, ".GR", genome.build, sep = ""))
      }
    }
  }else{
    COSMIC = as.data.frame(fread(Pmatrix))
    colnames(COSMIC)[1] = 'Subtype'
  }
  Ptmp = NULL
  for( i in 1:n )
  {
    ix = COSMIC$Subtype == subtypes$V1[i]
    Ptmp = rbind(Ptmp, COSMIC[ix,])
  }
  P = Ptmp[, -1] 
  Y = as.matrix(originalGenomes)
  X = normalize(Y)  
  error = steptol
  myexp = 10^-9 
  eps = .Machine$double.eps
  res = NULL
  try(load(paste("samplesResults/", AnalCOSMICSigType, ".Smatrix.RData", sep="")), silent=TRUE)  
  if(is.null(res))
  {	
    S = smoothSmatrix(X, P, eps)
    S = IterateISmatrix(X, S, P, error = error, myexp = myexp, eps = eps)
    Scutoff = S
    ix = Scutoff < cutoff
    Scutoff[ix] = 0
    Scutoff = normalize(Scutoff)
    colnames(S) = sampleNames[, 1]
    S = cbind(t(t(colnames(P))), S)
    colnames(S)[1] = 'SampleNames'	
    colnames(Scutoff) = sampleNames[, 1]
    Scutoff = cbind(t(t(colnames(P))), Scutoff)
    colnames(Scutoff)[1] = 'SampleNames'		
    res = list('S'=S, 'Scutoff'=Scutoff)
    save(res, file=paste("samplesResults/", AnalCOSMICSigType, ".Smatrix.RData", sep=""))	
  }else{
    S = res$S
    Scutoff = res$Scutoff
  }
  write.table(Ptmp, file=paste('samplesResults/', AnalCOSMICSigType, ".Pmatrix.txt",sep=""), quote=F, col.names=T, row.names=F, se="\t")  
  write.table(t(S), file=paste('samplesResults/', AnalCOSMICSigType, ".Smatrix.txt",sep=""), quote=F, col.names=F, row.names=T, se="\t") 
  write.table(t(Scutoff), file=paste('samplesResults/', AnalCOSMICSigType, ".Scutoffmatrix.txt",sep=""), quote=F, col.names=F, row.names=T, se="\t") 
  
  if(plot)
  {
    S0 = as.data.frame(fread(paste('samplesResults/', AnalCOSMICSigType, ".Smatrix.txt",sep="")))
    Scutoff0 = as.data.frame(fread(paste('samplesResults/', AnalCOSMICSigType, ".Scutoffmatrix.txt",sep="")))
        
    dataS0 = melt(S0, id = 'SampleNames')
    #png(filename = paste('samplesResults/',AnalCOSMICSigType,'.MutationalSignaturesAnalysis.png',sep=""),
    #width = 5500, height = 800)	
    pim <- ggplot(dataS0, aes(SampleNames, value, fill=variable)) +
      geom_bar(stat="identity") +
      guides(fill=guide_legend(reverse=F)) +
      scale_y_continuous(expand=c(0,0)) + 
      labs(x = "SampleID", y = "Sample contribution", title = "Mutational Signatures Analysis") + 
      theme(axis.title = element_text(size = 20), axis.text.x = element_text(angle = 90, size = 8, color = 'grey'), axis.text.y = element_text(angle = 0, size = 13, color = 'black'), panel.background = element_blank(), axis.line = element_line(colour = "grey"))
    pdf(file = paste('samplesResults/',AnalCOSMICSigType,'.MutationalSignaturesAnalysis.pdf',sep=""), 16, 8)  	    
    print(pim)
    dev.off()
    
    dataS0 = melt(Scutoff0, id = 'SampleNames')
    #png(filename = paste('samplesResults/',AnalCOSMICSigType,'.MutationalSignaturesAnalysisCutOff.png',sep=""),
    #width = 5500, height = 800)	
    pim <- ggplot(dataS0, aes(SampleNames, value, fill=variable)) +
      geom_bar(stat="identity") +
      guides(fill=guide_legend(reverse=F)) +
      scale_y_continuous(expand=c(0,0)) + 
      labs(x = "SampleID", y = "Sample contribution", title = "Mutational Signatures Analysis") + 
      theme(axis.title = element_text(size = 20), axis.text.x = element_text(angle = 90, size = 8, color = 'grey'), axis.text.y = element_text(angle = 0, size = 13, color = 'black'), panel.background = element_blank(), axis.line = element_line(colour = "grey"))
    pdf(file = paste('samplesResults/',AnalCOSMICSigType,'.MutationalSignaturesAnalysisCutOff.pdf',sep=""), 16, 8)
    print(pim)
    dev.off()
    
    A = S0
    A = A[,-1]
    n = ncol(A) + 1
    #png(filename = paste('samplesResults/',AnalCOSMICSigType,'.COSMICSignaturesWeight.png',sep=""),
    #width = 2400, height = 800)  
    pdf(file = paste('samplesResults/',AnalCOSMICSigType,'.COSMICSignaturesWeight.pdf',sep=""), 16, 4)
    plot(0,0,ylim=c(0,1),xlim=c(0,n),axes=F,type='l',xlab='',ylab='',cex=8)
    for(i in 1:dim(A)[2])
    {
      for(j in 1:dim(A)[1])
      {
        points(i,A[j,i],col='blue',pch=20)
      }
      a=A[,i]
      segments(i-0.5,median(a),i+0.5,median(a),col='black',lwd=2)
    }
    axis(2,at=seq(0,1,by=0.2),labels=seq(0,1,by=0.2), col="black",las=1,line=-4,cex.axis=1.2,font=2)
    axis(1,at=seq(1,n-1,by=1),labels=seq(1,n-1,by=1), col="black",las=1,cex.axis=1.2,font=2)
    mtext('COSMIC Signatures', side = 1, line = 4, outer = FALSE, at = NA,
          adj = NA, padj = NA, cex = 2, col = NA, font = NA)
    mtext('Relative weight', side = 2, line = 1.5, outer = FALSE, at = NA,
          adj = NA, padj = NA, cex = 2, col = NA, font = NA)    
    dev.off()
    
    A = Scutoff0
    A = A[,-1]
    n = ncol(A) + 1
    #png(filename = paste('samplesResults/',AnalCOSMICSigType,'.COSMICSignaturesWeightcutoff.png',sep=""),
    #width = 2400, height = 800)  
    pdf(file = paste('samplesResults/',AnalCOSMICSigType,'.COSMICSignaturesWeightcutoff.pdf',sep=""), 16, 4)
    plot(0,0,ylim=c(0,1),xlim=c(0,n),axes=F,type='l',xlab='',ylab='',cex=8)
    for(i in 1:dim(A)[2])
    {
      for(j in 1:dim(A)[1])
      {
        points(i,A[j,i],col='blue',pch=20)
      }
      a=A[,i]
      segments(i-0.5,median(a),i+0.5,median(a),col='black',lwd=2)
    }
    axis(2,at=seq(0,1,by=0.2),labels=seq(0,1,by=0.2), col="black",las=1,line=-4,cex.axis=1.2,font=2)
    axis(1,at=seq(1,n-1,by=1),labels=seq(1,n-1,by=1), col="black",las=1,cex.axis=1.2,font=2)
    mtext('COSMIC Signatures', side = 1, line = 4, outer = FALSE, at = NA,
          adj = NA, padj = NA, cex = 2, col = NA, font = NA)
    mtext('Relative weight', side = 2, line = 1.5, outer = FALSE, at = NA,
          adj = NA, padj = NA, cex = 2, col = NA, font = NA)    
    dev.off()          
  }
}
