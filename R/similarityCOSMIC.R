similarityCOSMIC = function(Pfile, AnalCOSMICSigType = c('SBS','DBS','ID'), genome.build = c('Ch37', 'Ch38'), SBS.version = c("V2","V3"), plot = TRUE, fontsize = 6, filename = NULL)
{
  library( data.table )
  library( reshape )
  library( pheatmap )
  if(is.null(Pfile))
  {
    stop("Please enter mutational signature matrix P with a format like the output result of RNMF software!\n")
  }   
  P = as.data.frame(fread(Pfile, header=T))
  colnames(P)[1] = 'Subtype'  
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
  if(AnalCOSMICSigType == 'ID')
  {
    data(COSMIC.ID.GRCh37)
    COSMIC = COSMIC[1:nrow(P), ]   
    data(IDsRaw)
    Subtype = ID
  }  
  if(AnalCOSMICSigType == 'DBS')
  {
	data(list=paste("COSMIC.", AnalCOSMICSigType, ".GR", genome.build, sep = ""))
    data(DBsRaw)
    Subtype = DBs    
  }  
  if(AnalCOSMICSigType == 'SBS')
  {
    if(SBS.version == 'V2')
    {
      data("COSMIC.SBS.V2")
    }else{      
	  data(list=paste("COSMIC.", AnalCOSMICSigType, ".GR", genome.build, sep = ""))
    }

    data(SBsRaw)
    Subtype = SBS  
  }  
  n = nrow(Subtype)
  Ptmp = NULL
  Cotmp = NULL
  for(i in 1:n)
  {
    ix = P$Subtype == Subtype[i,]
    iy = COSMIC$Subtype == Subtype[i,]
    Ptmp = rbind(Ptmp, P[ix,])
    Cotmp = rbind(Cotmp, COSMIC[iy,])
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
