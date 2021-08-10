SigsInput = function(file = NULL, Filetype = c('MAF','VCF'), AnalCOSMICSigType = c('SBS','DBS','ID'), genome.build = c("Ch37","Ch38"), sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', ID.mc.cores = 2, ID.row83 = TRUE )
{
  if(length(Filetype)>1)
  {
    Filetype = 'MAF'
  }
  if(length(genome.build)>1)
  {
    genome.build = 'Ch37'
  }  
  if(length(AnalCOSMICSigType)>1)
  {
    AnalCOSMICSigType = 'SBS'
  }  
  
  ## Load all dependencies ##
  if("data.table" %in% rownames(installed.packages()) == FALSE) 
  {   
    print("install [data.table]")
    install.packages("data.table")   
  }
  if("BiocManager" %in% rownames(installed.packages()) == FALSE)
  {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  }
  if("Biostrings" %in% rownames(installed.packages()) == FALSE)
  {  
    print("install [Biostrings]")
    BiocManager::install("Biostrings")
  }
  if("GenomicRanges" %in% rownames(installed.packages()) == FALSE) 
  {         
    print("install [GenomicRanges]")
    BiocManager::install("GenomicRanges")        
  }   
  if ("BSgenome.Hsapiens.UCSC.hg19" %in% rownames(installed.packages()) == FALSE & genome.build == "Ch37") 
  {         
    print("install [BSgenome.Hsapiens.UCSC.hg19]")
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")       
  }
  if ("BSgenome.Hsapiens.UCSC.hg38" %in% rownames(installed.packages()) == FALSE & genome.build == "Ch38") 
  {         
    print("install [BSgenome.Hsapiens.UCSC.hg38]")
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")    
  }
  
  res = NULL
  try(load("sigs.input.Rdata"),silent=TRUE)
  if(is.null(res))
  {
    if(is.null(file))
    {
      stop("InputFile is NULL!\n")
    }
    if(Filetype %in% 'MAF')
    {
      if(AnalCOSMICSigType == 'SBS')
      {
        res <- MafSigsInputSBS(file, genome.build, sample.id=sample.id, chr=chr, pos=pos, ref=ref, alt=alt)    
      }else
        if(AnalCOSMICSigType == 'ID')
        {
          res <- MafSigsInputID(file, genome.build, sample.id=sample.id, chr=chr, pos=pos, ref=ref, alt=alt, ID.mc.cores = ID.mc.cores)   
        }else{
          res <- MafSigsInputDBS(file, sample.id=sample.id, chr=chr, pos=pos, ref=ref, alt=alt)   	
        } 	
    }else{
      res <- VcfSigsInput(file, AnalCOSMICSigType, genome.build, sample.id=sample.id, chr=chr, pos=pos, ref=ref, alt=alt, ID.mc.cores = ID.mc.cores)
    }
    save(res, file="sigs.input.Rdata")
  }
  sigs.input = res$matrix
  mut = res$mut
  if( isTRUE(ID.row83) && AnalCOSMICSigType == 'ID' )
  {
			if(ncol(sigs.input) == 1)
			{
			    tmpcolN = colnames(sigs.input)
				sigs.input = t(t(sigs.input[1:83,]))
				 colnames(sigs.input) = tmpcolN
			}else{
				sigs.input = sigs.input[1:83,]			
			}
  }  
  sampleNames = colnames(sigs.input)
  subtypes = rownames(sigs.input)	
  write.table(sigs.input, file='originalGenomes', quote=F, sep="\t", col.names=F, row.names=F)	
  write.table(sampleNames, file='sampleNames', quote=F, sep="\t", col.names=F, row.names=F)
  write.table(subtypes, file='subtypes', quote=F, sep="\t", col.names=F, row.names=F)	
  write.table(mut, file='MutationInputSigsTypes.txt', quote=F, sep="\t", col.names=T, row.names=F)	 

}
