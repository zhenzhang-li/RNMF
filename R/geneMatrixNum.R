# geneListSortFile
# Hugo_Symbol  Chromosome


geneMatrixNum <- function(file = NULL, geneListSortFile = NULL, Filetype = c('MAF','VCF'), AnalCOSMICSigType = c('SBS','DBS','ID'), genome.build = c("Ch37","Ch38"), sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', Hugo = 'Hugo_Symbol', ID.mc.cores = 2, ID.row83 = TRUE )
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
  library(data.table)  
  
  # geneListSortFile
  #  Hugo_Symbol	Chromosome	Length
  # KDSR	chr18	39784	
  if(is.null(geneListSortFile))
  {
    data(geneList)
    #load("geneList.RData")
    geneList = geneList[, c('Hugo_Symbol','Chromosome','Length')]
  }else{
    geneList = as.data.frame(fread(geneListSortFile))
	data(knownGene)
	#load("knownGene.RData")
	geneList = merge(geneList, knownGene, by=1:2, all.x = TRUE, sort = FALSE)
  }
  colnames(geneList) = c(Hugo, chr, 'Length')
 
  if(exists("file", mode = "list")){
	  mut <- file
  } else {
	  mutdata = as.data.frame(fread(file))[, c(Hugo, chr, pos, ref, alt, sample.id)]
	  ix = grepl('_', mutdata[,chr]) | grepl('M', mutdata[,chr])
	  mut = mutdata[!ix,]
	  mut[,chr] = sub('chr', '', mut[,chr])
	  mut[,chr] = sub('Chr', '', mut[,chr])
	  mut[,chr] = sub('CHR', '', mut[,chr])
	  mut[,chr] = paste('chr', mut[,chr], sep='')  
  }  
   
  #mut = merge(geneList, mut, by=1:2, all=F)
  #AllgenesNames = geneList$Hugo_Symbol   
  #if(is.null(nrow(mut)))
  #{
  #  stop("There is no list of overlapping genes!\n")
  #}else
  #  if(nrow(mut) == 0)
  #  {
  #    stop("There is no list of overlapping genes!\n")
  #  }
  AllgenesNames = unique(sort(as.character(mut[, Hugo])))  
  giveGlist = as.character(geneList[, Hugo])
  intersectN = intersect(AllgenesNames, giveGlist)

  if(length(intersectN) == 0)
  {
    stop("There is no list of overlapping genes!\n")
  }
  
  SampleIDs = unique(sort(mut[, sample.id]))
  dir.create('samplesResults', showWarnings = FALSE, recursive = TRUE, mode = "0777")

  geneMatrixNumData = NULL
  try(load(paste("samplesResults/", AnalCOSMICSigType, ".geneMatrixNumData.Rdata", sep="")), silent=TRUE)
  if(is.null(geneMatrixNumData))
  {  
    nl = length(SampleIDs)
    for( i in 1:nl )
    {      
      dir.create(paste('samplesResults/', SampleIDs[i], "/", AnalCOSMICSigType, sep=""), showWarnings = FALSE, recursive = TRUE, mode = "0777")
      ix = mut[, sample.id] == SampleIDs[i]
      geneMatrixNumData$SampleID[[i]] = SampleIDs[i]
      muttmp = mut[ix,]
      
#      if(AnalCOSMICSigType == 'SBS')
#      {
#        sigs.input <- MafSigsInputSBS(muttmp, genome.build, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt)$matrix    
#      }else
#        if(AnalCOSMICSigType == 'ID')
#        {
#          sigs.input <- MafSigsInputID(muttmp, genome.build, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt, ID.mc.cores = ID.mc.cores)$matrix
#		  if( isTRUE(ID.row83))
#		  {
#			if(ncol(sigs.input) == 1)
#			{
#			    tmpcolN = colnames(sigs.input)
#				sigs.input = t(t(sigs.input[1:83,]))
#				 colnames(sigs.input) = tmpcolN
#			}else{
#				sigs.input = sigs.input[1:83,]			
#			}
#		  } 
#        }else{
#          sigs.input <- MafSigsInputDBS(muttmp, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt)$matrix 	
#        }  

    if(Filetype %in% 'MAF')
    {
      if(AnalCOSMICSigType == 'SBS')
      {
        sigs.input <- MafSigsInputSBS(muttmp, genome.build, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt)$matrix    
      }else
        if(AnalCOSMICSigType == 'ID')
        {
          sigs.input <- MafSigsInputID(muttmp, genome.build, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt, ID.mc.cores = ID.mc.cores)$matrix	  
        }else{
          sigs.input <- MafSigsInputDBS(muttmp, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt)$matrix 
        } 	
    }else{
      sigs.input <- VcfSigsInput(muttmp, AnalCOSMICSigType, genome.build, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt, ID.mc.cores = ID.mc.cores)$matrix
    }
	
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
      samGeneMatrix = matrix(0, dim(sigs.input)[1], length(AllgenesNames)+1)
 
      colnames(samGeneMatrix) = c('Subtype', AllgenesNames)
      geneCols = colnames(sigs.input)
      samGeneMatrix[, 1] = rownames(sigs.input)
      
      for(j in 1:length(AllgenesNames))
      {
        index = which(geneCols == AllgenesNames[j])
        if(length(index)>0)
        {
          samGeneMatrix[, j+1] = sigs.input[,index]
        }
      }
      geneMatrixNumData$samGeneMatrix[[i]] = samGeneMatrix
      write.table(samGeneMatrix, file=paste('samplesResults/', SampleIDs[i], "/", AnalCOSMICSigType, "/", SampleIDs[i], ".G.txt",sep=""), quote=F, col.names=T, row.names=F, se="\t")
    }
	geneMatrixNumData$geneList = geneList
    save(geneMatrixNumData, file=paste("samplesResults/", AnalCOSMICSigType, ".geneMatrixNumData.Rdata", sep=""))
  }  
  
  return(geneMatrixNumData)	
}
