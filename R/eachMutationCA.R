SigngeneMatrix <- function(mut, significantGenesList, AnalCOSMICSigType, genome.build, sample.id, chr, pos, ref, alt, Hugo, ID.mc.cores, ID.row83)
{
  data(knownGene)
  geneList = merge(significantGenesList, knownGene, by=1:2, all.x = TRUE, sort = FALSE)
  colnames(geneList) = c(Hugo, chr, 'Length')
  geneList = geneList[!is.na(geneList$Length), ]
  rownames(geneList) = geneList[, Hugo]
  
  AllgenesNames = unique(sort(as.character(mut[, Hugo])))  
  giveGlist = as.character(geneList[, Hugo])
  intersectN = intersect(AllgenesNames, giveGlist)

  if(length(intersectN) == 0)
  {
    stop("There is no list of overlapping genes!\n")
  }
  
  SampleIDs = unique(sort(mut[, sample.id]))
  geneMatrixNumData = NULL
  nl = length(SampleIDs)
  for( i in 1:nl )
  {      
      ix = mut[, sample.id] == SampleIDs[i]
      geneMatrixNumData$SampleID[[i]] = SampleIDs[i]
      muttmp = mut[ix,]
      
      if(AnalCOSMICSigType == 'SBS')
      {
        sigs.input <- MafSigsInputSBS(muttmp, genome.build, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt)$matrix    
      }else
        if(AnalCOSMICSigType == 'ID')
        {
          sigs.input <- MafSigsInputID(muttmp, genome.build, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt, ID.mc.cores = ID.mc.cores)$matrix
		  if( isTRUE(ID.row83))
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
        }else{
          sigs.input <- MafSigsInputDBS(muttmp, sample.id=Hugo, chr=chr, pos=pos, ref=ref, alt=alt)$matrix 	
        }   

      samGeneMatrix = matrix(0, dim(sigs.input)[1], length(intersectN)+1)
      colnames(samGeneMatrix) = c('Subtype', intersectN)
      geneCols = colnames(sigs.input)
      samGeneMatrix[, 1] = rownames(sigs.input)
      
      for(j in 1:length(intersectN))
      {
        index = which(geneCols == intersectN[j])
        if(length(index)>0)
        {
          samGeneMatrix[, j+1] = sigs.input[,index]
        }
      }
      geneMatrixNumData$samGeneMatrix[[i]] = samGeneMatrix
  }
  geneMatrixNumData$geneList = geneList[intersectN, ]  
  return(geneMatrixNumData)	
}

eachMutationCA <- function(file = NULL, Pfile = NULL, Sfile = NULL, sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', Hugo = 'Hugo_Symbol', Subtype = 'tricontext', AnalCOSMICSigType = c('SBS','DBS','ID'), genome.build = c("Ch37","Ch38"), significantGenesList = NULL, ID.mc.cores = 2, ID.row83 = TRUE) 
{
  if(is.null(file))
  {
    stop("Please enter a mutation dataset, such as the annotation result file in MAF or VCF format!\n")
	#Tumor_Sample_Barcode	Chromosome	Start_position	Hugo_Symbol	End_position	Variant_Classification	Reference_Allele	Tumor_Seq_Allele2	Protein_Change	context	mutcat	std.mutcat	std.context	tricontext
	#FP1705100059DN01	chr1	160000534	PIGM	160000534	Silent	G	A	p.T332T	AGG	G>A	C>T	CCT	C[C>T]T
  }
  if(is.null(Pfile))
  {
    stop("Please enter a mutational signature matrix file P with a format like the output result of RNMF software!\n")
  }  
  if(is.null(Sfile))
  {
    stop("Please enter a abundance fractions matrix S with a format like the output result of RNMF software!\n")
  }  
  
  library(data.table)  
  mutdata = as.data.frame(fread(file))
  ix = grepl('_', mutdata[,chr]) | grepl('M', mutdata[,chr])
  mut = mutdata[!ix,]
  mut[,chr] = sub('chr', '', mut[,chr])
  mut[,chr] = sub('Chr', '', mut[,chr])
  mut[,chr] = sub('CHR', '', mut[,chr])
  mut[,chr] = paste('chr', mut[,chr], sep='') 

  Nduplicate = mut[, c(Hugo, chr)]
  Nduplicate = unique(Nduplicate[order(Nduplicate[, Hugo]), ])
  Nduplicatenum = t(table(Nduplicate[, Hugo]))
  Nduplicate = colnames(Nduplicatenum)[which(Nduplicatenum>1)]
  index_ix = NULL
  for( i in 1:length(Nduplicate) )
  {
    ix = which(mut[, Hugo] == Nduplicate[i])
    if(length(ix)>0)
    {
      index_ix = c(index_ix, ix)
    }
  }
  if(!is.null(index_ix))
  {
    mut = mut[-index_ix, ]
  }

  samMatrixsigsData = samMatrixSigs(Pfile, Sfile, AnalCOSMICSigType)
  # significantGenesList
  # Hugo_Symbol  Chromosome
  if(!is.null(significantGenesList))
  {
	# Hugo_Symbol  Chromosome
    significantGenesListda = as.data.frame(fread(significantGenesList))
	colnames(significantGenesListda) = c(Hugo, chr)
	significantGenesList = merge(significantGenesListda, unique(mut[, c(Hugo, chr)]), by=c(Hugo, chr), all=F, sort=F)
  }else{
	significantGenesList = unique(mut[, c(Hugo, chr)])  
  }  
  colnames(significantGenesList) = c(Hugo, chr)
  geneMatrixNumData = SigngeneMatrix(mut, significantGenesList, AnalCOSMICSigType, genome.build, sample.id, chr, pos, ref, alt, Hugo, ID.mc.cores, ID.row83)
  geneList = geneMatrixNumData$geneList
  #  Hugo_Symbol	Chromosome	Length
  # KDSR	chr18	39784
  mut = merge(geneList[, c(Hugo, chr)], mut, by=c(Hugo, chr), all=F, sort=F)
 
  dir.create('samplesResults', showWarnings = FALSE, recursive = TRUE, mode = "0777")
  Samples = intersect(unlist(samMatrixsigsData$SampleID),unlist(geneMatrixNumData$SampleID))
  n = length(Samples)
  samMatrixsigsDatatmp = list()
  geneMatrixNumDatatmp = list()
  for(i in 1:n)
  {
	ix = which(unlist(samMatrixsigsData$SampleID) == Samples[i])
	iy = which(unlist(geneMatrixNumData$SampleID) == Samples[i])
	if(length(ix)>0 && length(iy)>0)
	{
		samMatrixsigsDatatmp$SampleID[[i]] = samMatrixsigsData$SampleID[[ix]]
		samMatrixsigsDatatmp$P[[i]] = samMatrixsigsData$P[[ix]]
		geneMatrixNumDatatmp$SampleID[[i]] = geneMatrixNumData$SampleID[[iy]]
		geneMatrixNumDatatmp$samGeneMatrix[[i]] = geneMatrixNumData$samGeneMatrix[[iy]]		
	}
  }
  samMatrixsigsData = samMatrixsigsDatatmp
  geneMatrixNumData = geneMatrixNumDatatmp
  rm("geneMatrixNumDatatmp","samMatrixsigsDatatmp")

  asnumeric = function(x)
  {
    x = as.numeric(as.character(x))
    return(x)
  }
  xy = function(x)
  {
    a = sum(x)
    if(a == 0)
    {
      a = 1
    }
    x = x/a
    return(x)
  }
  
  geneMutationsContributionAbundance = NULL
  geneRelativeMutationsContributionAbundance = NULL
  #  Hugo_Symbol	Chromosome	Length
  # KDSR	chr18	39784
  typesR = samMatrixsigsData$P[[1]][,1]
  ngSamples = as.character(unlist(samMatrixsigsData$SampleID))
  for(j in 1:length(ngSamples))
  {
	  tmpMutda = NULL
	  ixmut = mutdata[, sample.id] == ngSamples[j]
	  if(sum(ixmut)>0)
	  {
		  tmpMutda = mutdata[ixmut, ]
		  dir.create(paste('samplesResults/', ngSamples[j], "/", AnalCOSMICSigType, sep=""), showWarnings = FALSE, recursive = TRUE, mode = "0777")
		  index = which(samMatrixsigsData$SampleID == ngSamples[j])
		  rho = samMatrixsigsData$P[[index]]
		  rho = rho[, -1]
		  rho = apply(rho, 2, asnumeric)
		  index = which(geneMatrixNumData$SampleID == ngSamples[j])
		  if(length(index)>0)
		  {
			G = geneMatrixNumData$samGeneMatrix[[index]]
			G = G[, -1]
			G = apply(G, 2, asnumeric)	  
			G = t(G)  # gene * type		  
		  }else{
			G =  matrix(0, nrow(geneList), nrow(samMatrixsigsData$P[[1]]))
			rownames(G) = geneList[, Hugo]
		  }
		  
		  # 1
		  GG = apply(G, 2, xy) 
		  GG = t(apply(GG, 1, xy))	  
		  theta = GG%*%rho	
		  tmp = t(apply(theta, 1, xy))
		  Mttmp = t(GG)%*%tmp   # type * Sigs
		  thetatable = cbind(t(t(typesR)), Mttmp)
		  colnames(thetatable)[1] = Subtype
		  Out.mut = merge(tmpMutda, thetatable, by = c(Subtype), all.x = TRUE, all.y = FALSE, sort = FALSE)
		  write.table(Out.mut, file=paste('samplesResults/', ngSamples[j], "/", AnalCOSMICSigType, "/", ngSamples[j], ".SubtypeSigsResult.txt",sep=""), quote=F, col.names=T, row.names=F, sep="\t")
		  geneMutationsContributionAbundance = rbind(geneMutationsContributionAbundance, Out.mut)		
		 
		  # 2	  
		  G1 = G/geneList$Length
		  GG1 = apply(G1, 2, xy)
		  GG1 = t(apply(GG1, 1, xy))
		  gamma = GG1%*%rho
		  tmp = t(apply(gamma, 1, xy))	  
		  Mttmp = t(GG1)%*%tmp   # type * Sigs
		  gammatable = cbind(t(t(typesR)), Mttmp)	   
		  colnames(gammatable)[1] = Subtype
		  Out.mut = merge(tmpMutda, gammatable, by = c(Subtype), all.x = TRUE, all.y = FALSE, sort = FALSE)		  
		  write.table(Out.mut, file=paste('samplesResults/', ngSamples[j], "/", AnalCOSMICSigType, "/", ngSamples[j], ".SubtypeRelativeSigsResult.txt",sep=""), quote=F, col.names=T, row.names=F, sep="\t")
		  geneRelativeMutationsContributionAbundance = rbind(geneRelativeMutationsContributionAbundance, Out.mut)		
	  }
  }
  
  write.table(geneMutationsContributionAbundance, file=paste('samplesResults/', AnalCOSMICSigType, ".geneMutationsContributionAbundance.txt",sep=""), quote=F, col.names=T, row.names=F, sep="\t")
  write.table(geneRelativeMutationsContributionAbundance, file=paste('samplesResults/', AnalCOSMICSigType, ".geneRelativeMutationsContributionAbundance.txt",sep=""), quote=F, col.names=T, row.names=F, sep="\t") 
}
