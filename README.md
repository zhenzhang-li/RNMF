# RNMF
Mutational Signatures Analysis Tool

# Example
library(RNMF)

file='ESCC.exon.ALL.input.maf'

SigsInput(file, Filetype = 'MAF', AnalCOSMICSigType = 'SBS', genome.build = "Ch37", sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2')

originalGenomes = 'originalGenomes'

sampleNames = 'sampleNames'

subtypes = 'subtypes'

denovoNMF(originalGenomes, sampleNames, subtypes, kmin = 1, kmax = 15, steptol = 10^-9, mc.cores = 20)

Pfile = 'SBS.Pmatrix.txt'

Sfile = 'SBS.Smatrix.txt'

maffile = 'Regular.ESCC.exon.WES.maf-input-nonsilent.txt'

Groupfile = 'Groupfile.txt'

AnalCOSMICSigType = 'SBS'

cumulativeCA(file = maffile, Filetype = 'MAF', 
             Pfile = Pfile, Sfile = Sfile, 
             sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', 
             pos = 'Start_position', ref = 'Reference_Allele', 
             alt = 'Tumor_Seq_Allele2', Hugo = 'Hugo_Symbol', 
             AnalCOSMICSigType = AnalCOSMICSigType, 
             genome.build = "Ch37", 
             groupFile = Groupfile, 
             geneListSortFile = NULL, 
             ID.mc.cores = 1, 
             ID.row83 = TRUE, 
             plot = TRUE) 

Sigstype = c('SBS3','SBS16','SBS18','New','SBS33','SBS13','SBS2','SBS5','SBS1','SBS17b','SBS22','SBS15')

for(isig in Sigstype)
{

  Gfile = paste('samplesResults/',isig,'.',AnalCOSMICSigType,'.geneCumulativeContributionAbundance.txt',sep="")
  
  samFisherSigs(Sfile = Sfile, choose.Sigs = isig, file = maffile,
                sample.id = 'Tumor_Sample_Barcode', Hugo = 'Hugo_Symbol', 
                Mutfreq = 0.04, threshold = 0.06, Qvalue = 0.1, plot = TRUE) 
				
  samPerMutSigs(Sfile = Sfile, choose.Sigs = isig,
                file = maffile, sample.id = 'Tumor_Sample_Barcode',
                Hugo = 'Hugo_Symbol', Mutfreq = 0.04, 
                threshold = 0.06, Qvalue = 0.1, plot = TRUE) 
				
  genePerMutSigs(Sfile = Sfile, choose.Sigs = isig, 
                 Gfile = Gfile, file = maffile,
                 sample.id = 'Tumor_Sample_Barcode', Hugo = 'Hugo_Symbol', 
                 Mutfreq = 0.04, threshold = 0.06, Qvalue = 0.1, plot = TRUE) 
}
