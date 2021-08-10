VcfSigsInput <- function(mut, AnalCOSMICSigType, genome.build, sample.id = 'Tumor_Sample_Barcode', chr = 'Chromosome', pos = 'Start_position', ref = 'Reference_Allele', alt = 'Tumor_Seq_Allele2', ID.mc.cores = 2) 
{
  if(AnalCOSMICSigType == 'SBS')
  {
	sigs.input <- MafSigsInputSBS(mut, genome.build, sample.id=sample.id, chr=chr, pos=pos, ref=ref, alt=alt)    
  }else
  if(AnalCOSMICSigType == 'ID')
  {
	sigs.input <- MafSigsInputID(mut, genome.build, sample.id=sample.id, chr=chr, pos=pos, ref=ref, alt=alt, ID.mc.cores = ID.mc.cores)   
  }else{
	sigs.input <- MafSigsInputDBS(mut, sample.id=sample.id, chr=chr, pos=pos, ref=ref, alt=alt) 	
  }   
	
  return(sigs.input)
}
