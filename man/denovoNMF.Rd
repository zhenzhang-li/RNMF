\name{denovoNMF}
\alias{denovoNMF}
\alias{denovonmf}
\alias{NMF}
\alias{nmf}

\title{
	Non-Negative Matrix Factorization in cancer research.
}
\description{
	We present an optimized non-Negative Matrix Factorisation (NMF) based on R to extract the mutational singature in cancer research.
}
\usage{
	denovoNMF(originalGenomes, sampleNames, subtypes, kmin = 1, kmax = 10, AnalCOSMICSigType = 'SBS', steptol = 10^-9, totalIterations = 20, spacetime = 100, nlm = FALSE, mc.cores = 1)	
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
  \item{kmin}{
	Minimum classification of mutational signatures.
}
  \item{kmax}{
	Maximum classification of mutational signatures.
}
  \item{AnalCOSMICSigType}{       
     Type of mutation based on PCAWG project, canbe SBS, DBS or ID.
}
  \item{steptol}{
	A positive scalar providing the minimum allowable relative step length.
}
  \item{totalIterations}{
	The iterations number, must be at least 5. 
}
  \item{spacetime}{
	The length of definite solution space of NMF, must be at least 20.
}
  \item{nlm}{
	Straighten the result to calculate the optimal solution. Default FALSE.
}
  \item{mc.cores}{
	The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be exactly 1 on Windows (which uses the master process).
}
}
\author{
	Wen Luo(design, implementation, testing), Zhenzhang Li(design)
}
\examples{
	denovoNMF('originalGenomes', 'sampleNames', 'subtypes', kmin = 1, kmax = 10, AnalCOSMICSigType = 'SBS', steptol = 10^-9, totalIterations = 20, spacetime = 100, nlm = TRUE, mc.cores = 20)
}
