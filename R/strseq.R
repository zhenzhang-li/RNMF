strseq = function(chr,start,end,genome.build)
{
  if(genome.build %in% c('Ch37','hg19'))
  {
    posREF = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chr, as.numeric(start), as.numeric(end), as.character = T)
  }else{
    posREF = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, chr, as.numeric(start), as.numeric(end), as.character = T)	
  }
  return(as.character(posREF))
}