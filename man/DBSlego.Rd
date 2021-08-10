\name{DBSlego}
\alias{DBSlego}
\alias{dbslego}

\title{
	The tool is used to plot lego image for mutation spectrum.
}

\description{
   The tool is used to plot lego image for mutation spectrum, and its method references DBS Signatures.
}

\usage{
	DBSlego(file='originalGenomes', subtype='subtypes', genome.build = c("Ch37","Ch38"), SequenceType = c("WGS","WES"))
}

\arguments{
  \item{originalGenomes}{
	a file, data matrix for NMF, [m x n, m stands for feature counts and n stands for sample counts].
}
  \item{subtypes}{
	a list of types of mutation, data matrix for NMF, [m x 1, m stands for feature counts]. 
}
  \item{genome.build}{       
     The version of the reference genome used in the mutation set.
  }  
  \item{SequenceType}{       
     Type of Sequence, canbe WGS, WES.
  }  
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	DBSlego(file='originalGenomes', subtype='subtypes', genome.build = "Ch37", SequenceType = "WGS")
}