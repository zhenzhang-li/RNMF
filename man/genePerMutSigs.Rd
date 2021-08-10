\name{genePerMutSigs}
\alias{genePerMutSigs}
\alias{genepermutsigs}

\title{
	Mutation enrichment analysis.
}
\description{
	We present a method based on R for association analysis of genes and mutational signatures. Mutation enrichment analysis identifies an association between somatic mutations and activity of signature in a discovery cohort.
}
\usage{
	genePerMutSigs(Sfile = NULL, choose.Sigs = 'SBS1', Gfile = NULL, file = NULL, sample.id = 'Tumor_Sample_Barcode', Hugo = 'Hugo_Symbol', Mutfreq = 0.04, threshold = 0.06, Qvalue = 0.1, plot = TRUE)	
}
\arguments{
   \item{Sfile}{       
     A abundance fractions matrix S with a format like the output result of RNMF software.
  }
   \item{choose.Sigs}{       
     Select a feature for subsequent analysis, and the feature name should match the feature header name in Sfile.
  }
   \item{Gfile}{       
     A abundance fractions matrix of genes for all samples with a format like the output result of RNMF software.
  }  
   \item{file}{       
     A a non-silent mutation dataset, such as the annotation result file in MAF format.
  }    
  \item{sample.id}{       
     Column name in the mutation file corresponding to the Sample ID.	
  }
  \item{Hugo}{       
     Column name in the mutation file corresponding to the gene.
  }     
  \item{Mutfreq}{       
     A minimum sample mutation rate threshold is used to determine the genes for subsequent analysis.
  } 
  \item{threshold}{       
     A threshold is used to define the minimum contribution abundance of genes that can be displayed.
  }  
  \item{Qvalue}{       
     A FDR threshold of genes are used to highlight in red.
  }  
  \item{plot}{       
     Whether to plot the image. 
  } 
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	genePerMutSigs(Sfile = 'SampleContribution.12.Normalized.txt', choose.Sigs = 'SBS1', Gfile = 'SBS1.SBS.geneCumulativeContributionAbundance.txt', file = 'Non-silent.mutation.snv.txt', sample.id = 'Tumor_Sample_Barcode', Hugo = 'Hugo_Symbol', Mutfreq = 0.04, threshold = 0.06, Qvalue = 0.1, plot = TRUE)
}