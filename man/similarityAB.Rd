\name{similarityAB}
\alias{similarityAB}
\alias{similarityab}
\alias{simab}

\title{
	Similarity calculation.
}
\description{
	Calculate the similarities between features and draw heat maps.
}
\usage{
	similarityAB(Pfile1, Pfile2, AnalCOSMICSigType = c('SBS','DBS','ID'), plot = TRUE, fontsize = 6, filename = NULL) 	
}
\arguments{
  \item{Pfile1}{       
     A mutational signature matrix file A with a format like the output result of RNMF software.
  } 
  \item{Pfile2}{       
     A mutational signature matrix file B with a format like the output result of RNMF software.
  }   
  \item{AnalCOSMICSigType}{       
     Type of mutation based on PCAWG project, canbe SBS, DBS or ID.
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
	similarityAB(Pfile1 = 'A.12.Normalized.txt', Pfile2 = 'B.12.Normalized.txt', AnalCOSMICSigType = 'SBS', plot = TRUE, fontsize = 6, filename = NULL)
}
