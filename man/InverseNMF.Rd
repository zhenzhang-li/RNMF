\name{InverseNMF}
\alias{InverseNMF}
\alias{inversenmf}
\alias{INMF}
\alias{inmf}

\title{
	Determine the contribution of known mutational processes.
}
\description{
	We present an optimized tool based on R to determine the weights of each mutational signature contributing to an individual tumor sample.
}
\usage{
	InverseNMF(originalGenomes, sampleNames, subtypes, AnalCOSMICSigType = c('SBS','DBS','ID'), genome.build = c("Ch37","Ch38"), cutoff = 0.06, SBS.version = c("V2","V3"), Pmatrix = NULL, steptol = 10^-10, plot = TRUE)	
}
\arguments{
  \item{originalGenomes}{
	a file, data matrix for NMF, [m x n, m stands for feature counts and n stands for sample counts].
  }
  \item{sampleNames}{
	a list file of sample names , data matrix for NMF, [n x 1, n stands for sample counts].
  }
  \item{subtypes}{
	a list of types of mutation, data matrix for NMF, [m x 1, m stands for feature counts]. 
  }
  \item{AnalCOSMICSigType}{       
     Type of mutation based on PCAWG project, canbe SBS, DBS or ID.
  }
  \item{genome.build}{       
     The version of the reference genome used in the mutation set.
  }
  \item{cutoff}{       
     A threshold between 0 and 1, and was chosen to correct the weight that are most likely caused by errors.
  }
  \item{SBS.version}{       
     Mutational signatures version 3 and version 2 signatures in COSMIC.
  }
  \item{Pmatrix}{       
    A mutational signature matrix file P with a format like the output result of RNMF software.
  }
  \item{steptol}{
	A positive scalar providing the minimum allowable relative step length.
  }
  \item{plot}{       
    Whether to plot the image. 
  } 
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	InverseNMF('originalGenomes', 'sampleNames', 'subtypes', AnalCOSMICSigType = 'SBS', genome.build = "Ch37", cutoff = 0.06, SBS.version = "V2", Pmatrix = NULL, steptol = 10^-10, plot = TRUE)
}