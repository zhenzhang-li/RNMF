MafSigsInputID = function(mut.ref, genome.build = c("Ch37","Ch38"), sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', ID.mc.cores = 2)
{
  library( BSgenome )
  library( GenomeInfoDb )
  if(length(genome.build)>1)
  {
    genome.build = 'Ch37'
  }  
  if(genome.build %in% c('Ch37','hg19'))
  {
    library( BSgenome.Hsapiens.UCSC.hg19 )
    data(Ch37)
  }else{
    library( BSgenome.Hsapiens.UCSC.hg38 )
    data(Ch38)	
  }
  
  if(exists("mut.ref", mode = "list"))
  {
    mut.full <- mut.ref
  } else {
    if(file.exists(mut.ref)){
      library( data.table )  
      mut.full <- as.data.frame(fread(mut.ref))
    } else {
      stop("mut.ref is neither a file nor a loaded data frame")
    }
  }  
  mut <- mut.full[,c(sample.id, chr, pos, ref, alt)]
  ix = grepl('_', mut[,chr]) | grepl('M', mut[,chr])
  mut = mut[!ix,]
  mut[,ref] = toupper(mut[,ref])
  mut[,alt] = toupper(mut[,alt])  
  mut[,chr] = sub('chr', '', mut[,chr])
  mut[,chr] = sub('Chr', '', mut[,chr])
  mut[,chr] = sub('CHR', '', mut[,chr])  
  sampleAll <- unique(mut[, sample.id])
  indel_types = c('1:Del:C:0', '1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5',
                  '1:Del:T:0', '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5',
                  '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
                  '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', 
                  # >1bp INDELS
                  '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5',
                  '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5',
                  '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5',
                  '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5',
                  '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5', 
                  '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5', 
                  '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5',
                  '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5',
                  #MicroHomology INDELS
                  '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3',
                  '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5', '2:Ins:M:1', 
                  '3:Ins:M:1', '3:Ins:M:2', '4:Ins:M:1', '4:Ins:M:2', '4:Ins:M:3', '5:Ins:M:1', 
                  '5:Ins:M:2', '5:Ins:M:3', '5:Ins:M:4', '5:Ins:M:5', 'complex', 'non_matching')
  #indel_types = indel_types[1:83]	
  
  final.matrix = matrix(0, ncol = length(indel_types), nrow = length(sampleAll))
  colnames(final.matrix) = indel_types
  rownames(final.matrix) = sampleAll  

  ix = mut[, ref] %in% c('A', 'T', 'C', 'G') & mut[, alt] %in% c('A', 'T', 'C', 'G')   
  mut = mut[!ix,] 
  if(nrow(mut) > 0)
  {  
	  mut = mut[order(mut[,pos]),]
	  mut = mut[order(mut[,chr]),] 
	  mut = mut[order(mut[,sample.id]),] 
	  mut[,chr] = paste('chr', mut[,chr], sep='') 
		
	  
	  # remove the same deletions from one sample
	  # ACTGACGACGCT -->  1) ACTG[ACG]ACGCT;  2) ACTGACG[ACG]CT 
	  # ACTGACGTTACGCT -->  1) ACTGACG[T]TACGCT;  2) ACTGACGT[T]ACGCT 
	  mut$input = 'NO'
	  jdd = function(x)
	  {
		lengxx = length(strsplit(as.character(x[4]),'')[[1]])
		lengxx1 = length(strsplit(as.character(x[5]),'')[[1]])
		if(lengxx1==1)
		{
			if(lengxx<100 && as.character(x[4]) != '-' )
			{
				x[6] = 'YES'
			}
		} 
		return(x)
	  }
	  tmpmmut = as.data.frame(t(apply(mut,1,jdd)))
	  ix = tmpmmut[,c('input')] == 'YES'
	  tmpmmutkepp = tmpmmut[!ix,1:5]
	  tmpmmut = tmpmmut[ix,1:5]
	  newmut = NULL  
	#  for(samc in unique(tmpmmut[,sample.id]))
	#  {
	#	  ix = tmpmmut[,sample.id] == samc
	#	  mutsam = tmpmmut[ix,]	
	#	  for( ch in unique(mutsam[,chr]))
	#	  {
	#	    cat("Remove the same deletions from sample",samc,"on",ch,"\n")
	#		ix = mutsam[,chr] == ch
	#		tmpdata = mutsam[ix,]
	#		nlen = nrow(tmpdata)
	#		if(nlen>1)
	#		{
	#			  # make a index matrix
	#			  matrixindex = matrix(0,2,nlen-1)
	#			  matrixindex[1,] = 1:(nlen-1)
	#			  matrixindex[2,] = 2:nlen
	#			  indexix = NULL
	#			  for( kin in 1:(nlen-1) )
	#			  {
	#				  iix = matrixindex[1,kin]
	#				  jix = matrixindex[2,kin]
	#				  xx = tmpdata[iix,]
	#				  yy = tmpdata[jix,]
	#				  lengxx = length(strsplit(as.character(xx[4]),'')[[1]])
	#				  lengxx1 = length(strsplit(as.character(xx[5]),'')[[1]])
	#				  if( lengxx>=1 && lengxx1==1 && lengxx == length(strsplit(as.character(yy[4]),'')[[1]]) )
	#				  {
	#					if(as.numeric(as.character(xx[3]))+lengxx == as.numeric(as.character(yy[3])))
	#					{				    				  
	#						if(as.character(xx[5]) == '-' && as.character(yy[5]) == '-')
	#						{
	#		                    seq1 = strseq(as.character(xx[2]),as.numeric(as.character(xx[3])),as.numeric(as.character(xx[3]))+lengxx-1,genome.build)
	#							seq2 = strseq(as.character(yy[2]),as.numeric(as.character(yy[3])),as.numeric(as.character(yy[3]))+lengxx-1,genome.build)
	#						    if(seq1 == seq2)
	#							{
	#							   indexix = c(indexix, jix)							
	#							}
	#						}else{
	#						    lengxx = lengxx - 1
	#						    posnew1 = as.numeric(as.character(xx[3])) + 1
	#						    posnew2 = as.numeric(as.character(yy[3])) + 1
	#		                    seq1 = strseq(as.character(xx[2]),posnew1,posnew1+lengxx-1,genome.build)
	#							seq2 = strseq(as.character(yy[2]),posnew2,posnew2+lengxx-1,genome.build)
	#						    if(seq1 == seq2)
	#							{
	#							   indexix = c(indexix, jix)							
	#							}
	#						}
	#					}
	#				  }	
	#			  }
	#			  if(length(indexix)>0)
	#			  {
	#				newmut = rbind(newmut,tmpdata[-indexix,])
	#			  }else{
	#				newmut = rbind(newmut,tmpdata)
	#			  }
	#		}else{
	#		  newmut = rbind(newmut,tmpdata)
	#		}   
	#	  } 
	#  }
	  samplesset = unique(tmpmmut[,sample.id])	  
      newmut = NULL
	  if(length(samplesset) > 0)
	  {
	  x = 1:length(samplesset)
	  fistersample = function(i, samplesset, tmpmmut)
	  {
		  samc = samplesset[i]
		  ix = tmpmmut[,sample.id] == samc
		  mutsam = tmpmmut[ix,]	
		  newmut = NULL
		  for( ch in unique(mutsam[,chr]))
		  {
			cat("Remove the same deletions from sample",samc,"on",ch,"\n")
			ix = mutsam[,chr] == ch
			tmpdata = mutsam[ix,]
			nlen = nrow(tmpdata)
			if(nlen>1)
			{
				  # make a index matrix
				  matrixindex = matrix(0,2,nlen-1)
				  matrixindex[1,] = 1:(nlen-1)
				  matrixindex[2,] = 2:nlen
				  indexix = NULL
				  for( kin in 1:(nlen-1) )
				  {
					  iix = matrixindex[1,kin]
					  jix = matrixindex[2,kin]
					  xx = tmpdata[iix,]
					  yy = tmpdata[jix,]
					  lengxx = length(strsplit(as.character(xx[4]),'')[[1]])
					  lengxx1 = length(strsplit(as.character(xx[5]),'')[[1]])
					  if( lengxx>=1 && lengxx1==1 && lengxx == length(strsplit(as.character(yy[4]),'')[[1]]) )
					  {
						if(as.numeric(as.character(xx[3]))+lengxx == as.numeric(as.character(yy[3])))
						{				    				  
							if(as.character(xx[5]) == '-' && as.character(yy[5]) == '-')
							{
								seq1 = strseq(as.character(xx[2]),as.numeric(as.character(xx[3])),as.numeric(as.character(xx[3]))+lengxx-1,genome.build)
								seq2 = strseq(as.character(yy[2]),as.numeric(as.character(yy[3])),as.numeric(as.character(yy[3]))+lengxx-1,genome.build)
								if(seq1 == seq2)
								{
								   indexix = c(indexix, jix)							
								}
							}else{
								lengxx = lengxx - 1
								posnew1 = as.numeric(as.character(xx[3])) + 1
								posnew2 = as.numeric(as.character(yy[3])) + 1
								seq1 = strseq(as.character(xx[2]),posnew1,posnew1+lengxx-1,genome.build)
								seq2 = strseq(as.character(yy[2]),posnew2,posnew2+lengxx-1,genome.build)
								if(seq1 == seq2)
								{
								   indexix = c(indexix, jix)							
								}
							}
						}
					  }	
				  }
				  if(length(indexix)>0)
				  {
					newmut = rbind(newmut,tmpdata[-indexix,])
				  }else{
					newmut = rbind(newmut,tmpdata)
				  }
			}else{
			  newmut = rbind(newmut,tmpdata)
			}   
		  }
		 return(newmut)
	  }
	  res <- mclapply(x, fistersample, samplesset, tmpmmut, mc.cores = ID.mc.cores)
	  if(is.list(res))
	  {
		 for(i in 1:length(samplesset))
		 {
			newmut = as.data.frame(rbind(newmut,res[[i]]))
		 } 	
	  }
	  if(is.array(res) | is.matrix(res))
	  {
		 newmut = as.data.frame(res)
	  }	   
	  colnames(newmut) = colnames(tmpmmutkepp)
	  }
	  newmut = rbind(newmut, tmpmmutkepp)
	  rownames(newmut) = NULL
	  mut = as.data.frame(newmut)

	  xindex = function(i, mut, genome.build, chrom)
	  {
		x = mut[i,]
		x = as.matrix(x)		
		a = strsplit(as.character(x[4]),'')[[1]]
		b = strsplit(as.character(x[5]),'')[[1]]
		if( !(length(a) == 2 && length(b) == 2) && length(a) < 100 && length(b) < 100 )
		{
		  if(as.character(x[4]) != '-' && as.character(x[5]) != '-')
		  {
			# Ins		
			if(length(a)==1 && as.character(x[4]) == b[1])
			{
			  x[4] = '-'
			  x[5] = paste(b[2:length(b)],sep="")			
			}
			
			# Del				
			if(length(b)==1 && as.character(x[5]) == a[length(a)])
			{
			  x[3] = as.numeric(as.character(x[3])) + 1
			  x[4] = paste(a[1:(length(a)-1)],sep="")
			  x[5] = '-'					
			}		
		  }
		  #Ins or Del
		  if(as.character(x[4]) == '-' || as.character(x[5]) == '-') 
		  {
			chrlength = chrom$LEN[chrom$CHR == as.character(x[2])]
			x[6] = judgeType(x, genome.build, chrlength)
		  }else
			# complex
			if(!'N' %in% as.character(x[4]) && !'N' %in% as.character(x[5]))
			{
			  x[6] = 'complex'
			}else{
			  # N?
			  x[6] = 'non_matching'
			}	    
		}
		return(x)	
	  }
	  mut$tricontext = NA
	  x = 1:nrow(mut)
	  tmpmut = NULL
	  res <- mclapply(x, xindex, mut, genome.build, chrom, mc.cores = ID.mc.cores)
	  if(is.list(res))
	  {
		 for(i in 1:length(res))
		 {
			tmpmut = as.data.frame(rbind(tmpmut,res[[i]]))
		 } 	
	  }
	  if(is.array(res) | is.matrix(res))
	  {
		 tmpmut = as.data.frame(res)
	  }	 
	  
	  colnames(tmpmut) = c(sample.id, chr, pos, ref, alt, 'tricontext') 
	  ix = is.na(tmpmut$tricontext)
	  mut = tmpmut[!ix,]
	    
	  if(!is.null(mut))
	  {
		#final.matrix = matrix(0, ncol = length(indel_types), nrow = length(unique(mut[, sample.id])))
		#colnames(final.matrix) = indel_types
		#rownames(final.matrix) = unique(mut[, sample.id]) 
		for(i in unique(mut[, sample.id]))
		{
		  tmp = subset(mut, mut[, sample.id] == i)
		  beep = table(tmp$tricontext)
		  for(l in 1:length(beep))
		  {
			trimer = names(beep[l])
			if(trimer %in% indel_types)
			{
			  final.matrix[i, trimer] = beep[trimer]
			}
		  }
		}
	  }
  }
  final.df = data.frame(final.matrix, check.names = F)
  final.df = t(final.df)
  
  #mut.full[,c(sample.id, chr, pos, ref, alt)]
  #mut <- mut.full[,c(sample.id, chr, pos, ref, alt)]  
  mut.full[,chr] = paste('chr', mut.full[,chr], sep='')  
  mut = mut[,-c(4,5)]
  tmpmut = merge(mut.full,mut,by=c(sample.id, chr, pos),all=F,sort=F)
  return(list('mut'=tmpmut,'matrix'=final.df)) 
}
