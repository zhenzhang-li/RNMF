similarityAB = function(Pfile1, Pfile2, AnalCOSMICSigType = c('SBS','DBS','ID'), plot = TRUE, fontsize = 6, filename = NULL)
{
  library( data.table )
  library( reshape )
  library( pheatmap )
  if(is.null(Pfile1))
  {
    stop("Please enter mutational signature matrix A with a format like the output result of RNMF software!\n")
  }
   if(is.null(Pfile2))
  {
    stop("Please enter mutational signature matrix B with a format like the output result of RNMF software!\n")
  }   
  P1 = as.data.frame(fread(Pfile1, header=T))
  colnames(P1)[1] = 'Subtype'
  P2 = as.data.frame(fread(Pfile2, header=T))
  colnames(P2)[1] = 'Subtype'  
  dir.create('samplesResults', showWarnings = FALSE, recursive = TRUE, mode = "0777")  

  if(length(AnalCOSMICSigType)>1)
  {
    AnalCOSMICSigType = "SBS"
  }  
  if(AnalCOSMICSigType == 'ID')
  {
    data(IDsRaw)
    Subtype = ID
  }  
  if(AnalCOSMICSigType == 'DBS')
  {    
    data(DBsRaw)
    Subtype = DBs    
  }    
  if(AnalCOSMICSigType == 'SBS')
  {    
    data(SBsRaw)
    Subtype = SBS    
  }   
  n = nrow(Subtype)
  Ptmp = NULL
  Cotmp = NULL
  for(i in 1:n)
  {
    ix = P1$Subtype == Subtype[i,]
    iy = P2$Subtype == Subtype[i,]
    Ptmp = rbind(Ptmp, P1[ix,])
    Cotmp = rbind(Cotmp, P2[iy,])
  }
  coefficient = corf(Ptmp[,-1], Cotmp[, -1])    
  coefficienttmp = cbind(t(t(rownames(coefficient))), coefficient)
  colnames(coefficienttmp)[1] = 'Sigs'
  write.table(coefficienttmp, file=paste('samplesResults/', AnalCOSMICSigType, '.SimilarityHeatmap.txt',sep=""),quote=F, sep="\t", row.names = F, col.names = T)
  
  if(is.null(filename))
  {
    filename = paste('samplesResults/', AnalCOSMICSigType, '.SimilarityHeatmap.pdf', sep="")
  }
  if(plot)
  {    
    pheatmap(t(coefficient), display_numbers=T, fontsize=fontsize, filename = filename)
  }
}
