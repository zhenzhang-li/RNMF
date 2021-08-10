\name{SBSlego}
\alias{SBSlego}
\alias{sbslego}
\alias{lego}
\alias{Lego}

\title{
	The tool is used to plot lego image for mutation spectrum.
}

\description{
   The tool is used to plot lego image for mutation spectrum, and its method comes from lwlegopt(https://github.com/BGI-LuoWen/lwlegopt).
}

\usage{
	SBSlego(file='originalGenomes', subtype='subtypes', scale=10, name=NULL, title='GC', sort=NULL, top=TRUE, color=NULL, border=NULL, genome.build = c("Ch37","Ch38"), SequenceType = c("WGS","WES"), RegionLength = NULL)
}

\arguments{
  \item{originalGenomes}{
	a file, data matrix for NMF, [m x n, m stands for feature counts and n stands for sample counts].
}
  \item{subtypes}{
	a list of types of mutation, data matrix for NMF, [m x 1, m stands for feature counts]. 
}
  \item{scale}{Scale interval of Y axis.}
  \item{name}{Prefix_Name of output file.}
  \item{title}{Cancer type.}
  \item{sort}{Sort for trinucleotide percentage value: TRUE, FALSE and NULL.}
  \item{top}{Location of icon displaying: TRUE, FALSE.}
  \item{color}{Colours for lego figure.}
  \item{border}{Colours for border.}
  \item{genome.build}{       
     The version of the reference genome used in the mutation set.
  }  
  \item{SequenceType}{       
     Type of Sequence, canbe WGS, WES.
  }  
  \item{RegionLength}{       
     The length of sequencing region.
  }    
}
\author{
	Wen Luo(design, implementation, testing)
}
\examples{
	SBSlego(file='originalGenomes', subtype='subtypes', scale=10, name=NULL, title='GC', sort=NULL, top=TRUE, color=NULL, border=NULL, genome.build = "Ch37", SequenceType = "WGS")
}