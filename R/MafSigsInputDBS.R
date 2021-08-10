MafSigsInputDBS = function(mut.ref, sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2')
{

  if(exists("mut.ref", mode = "list")){
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
  mut[,chr] = sub('chr', '', mut[,chr])
  mut[,chr] = sub('Chr', '', mut[,chr])
  mut[,chr] = sub('CHR', '', mut[,chr])
  mut[,chr] = paste('chr', mut[,chr], sep='')
  mut[,ref] = toupper(mut[,ref])
  mut[,alt] = toupper(mut[,alt])   
  mutation_types = c('CC>AA','CC>AG','CC>AT','CC>GA','CC>GG','CC>GT','CC>TA','CC>TG','CC>TT',
						  'CT>AA','CT>AC','CT>AG','CT>GA','CT>GC','CT>GG','CT>TA','CT>TC','CT>TG',
						  'TC>AA','TC>AG','TC>AT','TC>CA','TC>CG','TC>CT','TC>GA','TC>GG','TC>GT',
						  'TT>AA','TT>AC','TT>AG','TT>CA','TT>CC','TT>CG','TT>GA','TT>GC','TT>GG')
  mutation_types_non_tsb = c('AC>CA','AC>CG','AC>CT','AC>GA','AC>GG','AC>GT','AC>TA','AC>TG','AC>TT',
						  'AT>CA','AT>CC','AT>CG','AT>GA','AT>GC','AT>TA',
						  'CG>AT','CG>GC','CG>GT','CG>TA','CG>TC','CG>TT',
						  'GC>AA','GC>AG','GC>AT','GC>CA','GC>CG','GC>TA',
						  'TA>AT','TA>CG','TA>CT','TA>GC','TA>GG','TA>GT',
						  'TG>AA','TG>AC','TG>AT','TG>CA','TG>CC','TG>CT','TG>GA','TG>GC','TG>GT')
  mut78type = sort(c(mutation_types_non_tsb,mutation_types)) 

  sampleAll <- unique(mut[, sample.id])
  final.matrix = matrix(0, ncol = 78, nrow = length(sampleAll))
  colnames(final.matrix) = mut78type
  rownames(final.matrix) = sampleAll  

  # remove those length more than 2 and less than 2
  xindex = function(x, mut78type)
  {
	x = as.character(x)
	a = strsplit(x[4],'')[[1]]
	b = strsplit(x[5],'')[[1]]
	if(length(a) == 2 && length(b) == 2)
	{
		x[6] = paste(x[6],'>',sep="")
		if(! x[6] %in% mut78type)
		{
			x[6] = paste(reverse(x[4]),'>',sep="")
			x[7] = reverse(x[5])
		}
		return(t(x))
	}
  }
  mut$refcat = mut[, ref]
  mut$altcat = mut[, alt]
  tmpmut = apply(mut, 1, xindex, mut78type)
  
  if(!is.null(tmpmut))
  {
	  mut = NULL
	  for(i in 1:length(tmpmut))
	  {
		mut = as.data.frame(rbind(mut, tmpmut[[i]]))
	  }
	  colnames(mut) = c(sample.id, chr, pos, ref, alt, 'refcat', 'altcat')  
	  mut$tricontext = paste(mut$refcat,mut$altcat,sep="")
	  
	  #final.matrix = matrix(0, ncol = 78, nrow = length(unique(mut[, sample.id])))
	  #colnames(final.matrix) = mut78type
	  #rownames(final.matrix) = unique(mut[, sample.id]) 
	  for(i in unique(mut[, sample.id]))
	  {
		tmp = subset(mut, mut[, sample.id] == i)
		beep = table(tmp$tricontext)
		for(l in 1:length(beep))
		{
		  trimer = names(beep[l])
		  if(trimer %in% mut78type)
		  {
			final.matrix[i, trimer] = beep[trimer]
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
