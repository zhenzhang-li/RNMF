MafSigsInputSBS = function(mut.ref, genome.build = c("Ch37","Ch38"), sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2')
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
  }else{
    library( BSgenome.Hsapiens.UCSC.hg38 )
  }
  
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
  
  # Generate all possible trinucleotide contexts
  all.tri = c()
  for(i in c("A", "C", "G", "T"))
  {
    for(j in c("C", "T"))
	{
      for(k in c("A", "C", "G", "T"))
	  {
        if(j != k){
          for(l in c("A", "C", "G", "T"))
		  {
            tmp = paste(i, "[", j, ">", k, "]", l, sep = "")
            all.tri = c(all.tri, tmp)
          }
        }
      }
    }
  }
  all.tri <- all.tri[order(substr(all.tri, 3, 5))]  
  sampleAll <- unique(mut[, sample.id])
  final.matrix = matrix(0, ncol = 96, nrow = length(sampleAll))
  colnames(final.matrix) = all.tri
  rownames(final.matrix) = sampleAll  
  mut <- mut[which(mut[, ref] %in% c('A', 'T', 'C', 'G') & mut[, alt] %in% c('A', 'T', 'C', 'G')),] 
  if(nrow(mut) > 0)
  {
	  if(genome.build %in% c('Ch37','hg19'))
	  {
		 mut <- mut[mut[, chr] %in% GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens), ]
		 mut$context = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, mut[, chr], mut[, pos]-1, mut[, pos]+1, as.character = T)  
	  }else{
		mut <- mut[mut[, chr] %in% GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg38::Hsapiens), ]
		mut$context = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, mut[, chr], mut[, pos]-1, mut[, pos]+1, as.character = T)    
	  } 
	  mut$mutcat = paste(mut[, ref], ">", mut[, alt], sep = "")

	  # Reverse complement the G's and A's
	  gind = grep("G",substr(mut$mutcat,1,1))
	  tind = grep("A",substr(mut$mutcat,1,1))  
	  mut$std.mutcat = mut$mutcat
	  mut$std.mutcat[c(gind, tind)] <- gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", mut$std.mutcat[c(gind, tind)])))) # to lowercase
	  mut$std.mutcat[c(gind, tind)] <- gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", mut$std.mutcat[c(gind, tind)])))) # complement
	  mut$std.context = mut$context
	  mut$std.context[c(gind, tind)] <- gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", mut$std.context[c(gind, tind)])))) # to lowercase
	  mut$std.context[c(gind, tind)] <- gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", mut$std.context[c(gind, tind)])))) # complement
	  mut$std.context[c(gind, tind)] <- sapply(strsplit(mut$std.context[c(gind, tind)], split = ""), function(str) {paste(rev(str), collapse = "")}) # reverse
	  # Make the tricontext
	  mut$tricontext = paste(substr(mut$std.context, 1, 1), "[", mut$std.mutcat, "]", substr(mut$std.context, 3, 3), sep = "")


  
	  #final.matrix = matrix(0, ncol = 96, nrow = length(unique(mut[, sample.id])))
	  #colnames(final.matrix) = all.tri
	  #rownames(final.matrix) = unique(mut[, sample.id]) 
	  for(i in unique(mut[, sample.id]))
	  {
		tmp = subset(mut, mut[, sample.id] == i)
		beep = table(tmp$tricontext)
		for(l in 1:length(beep))
		{
		  trimer = names(beep[l])
		  if(trimer %in% all.tri)
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
