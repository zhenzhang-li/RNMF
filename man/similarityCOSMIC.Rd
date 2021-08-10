\name{similarityCOSMIC}
\alias{similarityCOSMIC}
\alias{similaritycosmic}
\alias{simcosmic}

\title{
	Similarity calculation based on COSMIC Mutational Signatures.
}
\description{
	Calculate the similarities between features and draw heat maps.
}
\usage{
	similarityCOSMIC(Pfile, AnalCOSMICSigType = c('SBS','DBS','ID'), genome.build = c('Ch37', 'Ch38'), SBS.version = c("V2","V3"), plot = TRUE, fontsize = 6, filename = NULL) 	
}
\arguments{
  \item{Pfile}{       
     A mutational signature matrix file P with a format like the output result of RNMF software.
  } 
  \item{AnalCOSMICSigType}{       
     Type of mutation based on PCAWG project, canbe SBS, DBS or ID.
  }
  \item{genome.build}{       
     The version of the reference genome used in the mutation set.
  }
  \item{SBS.version}{       
     Mutational signatures version 3 and version 2 signatures in COSMIC.
  } 
  \item{plot}{       
    Whether to plot the image. 
  }  
  \item{fontsize}{       
    base fontsize for the plot. 
  }    
  \item{filename}{       
    file path where to save the picture. Filetype is decided by the extension in the path. Currently following formats are supported: png, pdf, tiff, bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is calculated so that the plot would fit there, unless specified otherwise.
  }    
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	similarityCOSMIC(Pfile = 'SignatureComposition.12.Normalized.txt', AnalCOSMICSigType = 'SBS', genome.build = 'Ch37', SBS.version = "V2", plot = TRUE, fontsize = 6, filename = NULL)
}
