\name{eachMutationCA}
\alias{eachMutationCA}
\alias{eachMutationca}

\title{
	Contribution abundance of each mutation.
}
\description{
	We present a tool to calculate the contribution abundance of each mutation in cancer research.
}
\usage{
	eachMutationCA(file = NULL, Pfile = NULL, Sfile = NULL, sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', Hugo = 'Hugo_Symbol', Subtype = 'tricontext', AnalCOSMICSigType = c('SBS','DBS','ID'), genome.build = c("Ch37","Ch38"), significantGenesList = NULL, ID.mc.cores = 2, ID.row83 = TRUE) 	
}
\arguments{
  \item{file}{       
     A file with non-silent mutations obtained by annotation of mutation locations by Oncotator or other annotation software. And then can be generated through the SigsInput function. The format of the input dataset, contains MAF and VCF. File must have at least seven columns with the colnames of "Hugo_Symbol	CHROM	POS	REF	ALT	SAMPLENAME	SUBTYPE", 
	 for SBS  as: 1    1432333       T       C	T1	C[C>T]T, 
	 for DBS  as: 1    1432333       GC       AA	T1	GC>AA, 
	 for ID  as: 1    1432333       T       -	T1	1:Del:T:3.
  }
  \item{Pfile}{       
     A mutational signature matrix file P with a format like the output result of RNMF software.
  }
   \item{Sfile}{       
     A abundance fractions matrix S with a format like the output result of RNMF software.
  }
  \item{sample.id}{       
     Column name in the mutation file corresponding to the Sample ID.	
  }
  \item{chr}{       
     Column name in the mutation file corresponding to the chromosome.
  }  
  \item{pos}{       
     Column name in the mutation file corresponding to the mutation position.
  }  
  \item{ref}{       
     Column name in the mutation file corresponding to the reference base.
  }
  \item{alt}{       
     Column name in the mutation file corresponding to the alternate base.
  }
  \item{Hugo}{       
     Column name in the mutation file corresponding to the gene.
  }    
  \item{Subtype}{       
     Column name in the mutation file corresponding to the Subtype.
  }
  \item{AnalCOSMICSigType}{       
     Type of mutation based on PCAWG project, canbe SBS, DBS or ID.
  }
  \item{genome.build}{       
     The version of the reference genome used in the mutation set.
  }
  \item{significantGenesList}{ 
	Significant gene list with two columns: 'Hugo_Symbol', 'Chromosome', such as "KDSR      chr18".
  }   
  \item{ID.mc.cores}{       
     The number of cores to use in ID type, i.e. at most how many child processes will be run simultaneously. Must be exactly 1 on Windows (which uses the master process).
  }
  \item{ID.row83}{       
    Whether to keep 83 properties for ID type. 
  }
}
\author{
	Wen Luo (design), Zhenzhang Li (implementation, testing)
}
\examples{
	eachMutationCA(file = 'MutationInputSigsTypes.txt', Pfile = 'SignatureComposition.12.Normalized.txt', Sfile = 'SampleContribution.12.Normalized.txt', sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', Hugo = 'Hugo_Symbol', Subtype = 'tricontext',  AnalCOSMICSigType = 'SBS')
}
