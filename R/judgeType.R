insf = function(x,num1, num2)
{
  fseq = NULL
  for( i in num1:num2 )
  {
	fseq = paste(fseq,x[i],sep="")
  }
  return(fseq)
}

judgeType = function(x, genome.build, chrlength)
{
  repnum = 0
  type1 = 'Del'
  type2 = 'R'
  nlength = 0
  a = strsplit(as.character(x[4]),'')[[1]]
  b = strsplit(as.character(x[5]),'')[[1]]	
  # Ins
  if(as.character(x[4]) == '-')
  {
    nlength = length(b)
    seqreverse = as.character(x[5])
    type1 = 'Ins'
    if(nlength==1)
    {
      if(as.character(x[5]) %in% c('C','G'))
      {
        type2 = 'C'
      }else{
        type2 = 'T'
      }			
    }
  }
  # Del
  if(as.character(x[5]) == '-')
  {
    nlength = length(a)
    type1 = 'Del'
    seqreverse = as.character(x[4])
    if(nlength==1)
    {
      if(as.character(x[4]) %in% c('C','G'))
      {
        type2 = 'C'
      }else{
        type2 = 'T'
      }			
    }
    refseq = strseq(as.character(x[2]), as.numeric(as.character(x[3])), as.numeric(as.character(x[3]))+nlength-1, genome.build)
    if(seqreverse != refseq)
    {
      seqreverse = seqs(seqreverse)
    }		
  }
  
  pos = as.numeric(as.character(x[3]))
  # judge repetition 	
  if( type1 == 'Del' )
  {
    # >=1bp Del
    ## del on left	
    start = pos + nlength
    end = start + nlength - 1
    repseq = strseq(as.character(x[2]), start, end, genome.build)
    while(end<=chrlength && repseq == seqreverse)
    {
      repnum = repnum + 1
      start = end + 1
      end = start + nlength - 1
      if(end>chrlength)
      {
        end = chrlength
      }
      repseq = strseq(as.character(x[2]), start, end, genome.build)
    }
    ## del on rigth
    end = pos - 1
    start = end - nlength + 1
    repseq = strseq(as.character(x[2]), start, end, genome.build)
    while(start>0 && repseq == seqreverse)
    {
      repnum = repnum + 1
      end = start - 1
      start = end - nlength + 1
      if(start<=0)
      {
        start = 1
      }
      repseq = strseq(as.character(x[2]), start, end, genome.build)
    }

    # MicroHomology Del 5' ---> 3'
	# C[ATCTCAT]CTCAT: According to the article, this type does not take into account
	# ACCCA[TTC]pTpAGCGGC or ACCCpCp[TTC]AAGCGGC 
    # >1bp Del
    if( nlength > 1 && repnum == 0 )
    {	
	  ## located at the 3‘ end of the deletion
	  repnuml = 0  
      start = pos + nlength
      end = start + nlength - 1
      mhl = mhl1 = 1
      delSseq = strseq(as.character(x[2]), pos, pos, genome.build)
	  NSseq = strseq(as.character(x[2]), start, start, genome.build)
      if( delSseq != NSseq )
      {
        mhl1 = 0
      }	  
      while(mhl1 && mhl && end>=start)
      {
        end = end - 1
        repnuml = end - start + 1
        repseq = strseq(as.character(x[2]), start, end, genome.build)
        if(grepl(repseq, seqreverse))
        {
			seqreversetmp = strseq(as.character(x[2]), pos, pos+repnuml-1, genome.build)
			if(repseq == seqreversetmp)
			{
			  mhl = 0
			}
		}
      }
	  
	  ## located at the 5’ end of the deletion
      repnumr = 0			
      end = pos - 1
      start = end - nlength + 1
      mhr = mhr1 = 1
      delSseq = strseq(as.character(x[2]), pos+nlength-1, pos+nlength-1, genome.build)
	  NSseq = strseq(as.character(x[2]), end, end, genome.build)
      if( delSseq != NSseq )
      {
        mhr1 = 0
      }	  	    
      while(mhr1 && mhr && end>=start)
      {
        start = start + 1
        repnumr = end - start + 1
        repseq = strseq(as.character(x[2]), start, end, genome.build)
        if(grepl(repseq, seqreverse))
		{
			seqreversetmp = strseq(as.character(x[2]), pos+nlength-1-repnumr+1, pos+nlength-1, genome.build)
			if(repseq == seqreversetmp)
			{
			  mhr = 0
			}
		}
      }	  
      if( repnuml>=repnumr )
      {
        if(mhl==0)
        {
          repnum = repnuml				
        }
      }else{
        if(mhr==0)
        {
          repnum = repnumr				
        }
      }
      if(repnum>0)
      {
        type2 = 'M'						
      }	  
      #repnuml = 0
      #start = pos + nlength
      #end = start + nlength - 1
      #mhl = 1
      #while(mhl && end>=start)
      #{
      #  end = end - 1
      #  repnuml = end - start + 1
      #  repseq = strseq(as.character(x[2]), start, end, genome.build)
      #  if(repseq %in% seqreverse)
      #  {
      #    mhl = 0
      #  }			
      #}
      ### del on rigth	
      #repnumr = 0			
      #end = pos - 1
      #start = end - nlength + 1
      #mhr = 1
      #while(mhr && end>=start)
      #{
      #  start = start + 1
      #  repnumr = end - start + 1
      #  repseq = strseq(as.character(x[2]), start, end, genome.build)
      #  if(repseq %in% seqreverse)
      #  {
      #    mhr = 0
      #  }			
      #}
      #if( repnuml>=repnumr )
      #{
      #  if(mhl==0)
      #  {
      #    repnum = repnuml				
      #  }
      #}else{
      #  if(mhr==0)
      #  {
      #    repnum = repnumr				
      #  }
      #}
      #if(repnum>0)
      #{
      #  type2 = 'M'						
      #}	  
    }
	#else{
	#  repnum = repnum + 1
	#}	
  }
  
  # judge repetition 	
  if( type1 == 'Ins' )
  {	
    # >=1bp Ins
    ## Ins on left	
    start = pos
    end = start + nlength - 1
    repseq = strseq(as.character(x[2]), start, end, genome.build)
    while(end<=chrlength && repseq == as.character(x[5]))
    {
      repnum = repnum + 1
      start = end + 1
      end = start + nlength - 1
      if(end>chrlength)
      {
        end = chrlength
      }
      repseq = strseq(as.character(x[2]), start, end, genome.build)
    }
    ## Ins on rigth
    end = pos - 1
    start = end - nlength + 1
    repseq = strseq(as.character(x[2]), start, end, genome.build)
	repseq1 = strseq(as.character(x[2]), as.numeric(as.character(x[3])), as.numeric(as.character(x[3])), genome.build)
    while(start>0 && repseq1 == as.character(x[5]) && repseq == as.character(x[5]))
    {
      repnum = repnum + 1
      end = start - 1
      start = end - nlength + 1
      if(start<=0)
      {
        start = 1
      }
      repseq = strseq(as.character(x[2]), start, end, genome.build)
    }		
    #if(repnum == 0)
    #{
      #seqreverse = seqs(seqreverse)	
      ## Ins on left	
      #start = pos
      #end = start + nlength - 1
      #repseq = strseq(as.character(x[2]), start, end, genome.build)
      #while(end<=chrlength && repseq == seqreverse)
      #{
      #  repnum = repnum + 1
      #  start = end + 1
      #  end = start + nlength - 1	
      #  if(end>chrlength)
      #  {
      #    end = chrlength
      #  }			
      #  repseq = strseq(as.character(x[2]), start, end, genome.build)
      #}
      ## Ins on rigth
      #end = pos - 1
      #start = end - nlength + 1
      #repseq = strseq(as.character(x[2]), start, end, genome.build)
	  #repseq1 = strseq(as.character(x[2]), as.numeric(as.character(x[3])), as.numeric(as.character(x[3])), genome.build)
      #while(start>0 && repseq1 == seqreverse && repseq == seqreverse)
      #{
      #  repnum = repnum + 1
      #  end = start - 1
      #  start = end - nlength + 1
      #  if(start<=0)
      #  {
      #    start = 1
      #  }
      #  repseq = strseq(as.character(x[2]), start, end, genome.build)
      #}		
    #}
 
    #MicroHomology Ins 5' --> 3'
    # >1bp Ins
    #if( nlength > 1 && repnum == 0 )
    #{
      ## Ins on left
      #repnuml = 0
      #start = pos + 1
      #end = start + nlength - 1
      #mhl = 1
      #while(mhl && end>=start)
      #{
      #  end = end - 1
      #  repnuml = end - start + 1
      #  repseq = strseq(as.character(x[2]), start, end, genome.build)
      #  if(repseq %in% as.character(x[5]))
      #  {
      #    mhl = 0
      #  }			
      #}
      ## Ins on rigth	
      #repnumr = 0			
      #end = pos - 1
      #start = end - nlength + 1
      #mhr = 1
      #while(mhr && end>=start)
      #{
      #  start = start + 1
      #  repnumr = end - start + 1
      #  repseq = strseq(as.character(x[2]), start, end, genome.build)
      #  if(repseq %in% as.character(x[5]))
      #  {
      #    mhr = 0
      #  }			
      #}
      #if( repnuml>=repnumr )
      #{
      #  if(mhl==0)
      #  {
      #    repnum = repnuml				
      #  }
      #}else{
      #  if(mhr==0)
      #  {
      #    repnum = repnumr				
      #  }
      #}
      #if(repnum == 0)
      #{
      #  ## Ins on left
      #  repnuml = 0
      #  start = pos + 1
      #  end = start + nlength - 1
      #  mhl = 1
      #  while(mhl && end>=start)
      #  {
      #   end = end - 1
      #    repnuml = end - start + 1
      #    repseq = strseq(as.character(x[2]), start, end, genome.build)
      #    if(repseq %in% seqreverse)
      #    {
      #      mhl = 0
      #    }			
      #  }
      #  ## Ins on rigth	
      #  repnumr = 0			
      #  end = pos - 1
      #  start = end - nlength + 1
      #  mhr = 1
      #  while(mhr && end>=start)
      #  {
      #    start = start + 1
      #    repnumr = end - start + 1
      #    repseq = strseq(as.character(x[2]), start, end, genome.build)
      #    if(repseq %in% seqreverse)
      #    {
      #     mhr = 0
      #    }			
      #  }
      #  if( repnuml>=repnumr )
      #  {
      #    if(mhl==0)
      #    {
      #      repnum = repnuml				
      #    }
      #  }else{
      #    if(mhr==0)
      #    {
      #      repnum = repnumr				
      #    }
      #  }			
      #}
      #if(repnum>0)
      #{
      #  type2 = 'M'						
      #}		
    #}
    #if( nlength > 1 && repnum == 0 )
    #{
		#type2 = 'NA'
	#}
	
    #MicroHomology Ins 5' --> 3'
    # >1bp Ins
    if( nlength > 1 && repnum == 0 )
    {	
	  ## located at the 3‘ end of the deletion
	  repnuml = 0  
      start = pos
      end = start + nlength - 1
      mhl = mhl1 = 1
      InsSseq = strseq(as.character(x[2]), pos, pos, genome.build)
	  NSseq = strsplit(as.character(x[5]),'')[[1]][1]
      if( InsSseq != NSseq )
      {
        mhl1 = 0
      }	  
      while(mhl1 && mhl && end>=start)
      {
        end = end - 1
        repnuml = end - start + 1
        repseq = strseq(as.character(x[2]), start, end, genome.build)
        if(grepl(repseq, as.character(x[5])))
        {
			seqreversetmp = insf(strsplit(as.character(x[5]),'')[[1]], 1, repnuml)
			if(repseq == seqreversetmp)
			{
			  mhl = 0
			}
		}
      }
	  
	  ## located at the 5’ end of the deletion
      repnumr = 0			
      end = pos - 1
      start = end - nlength + 1
      mhr = mhr1 = 1
      InsSseq = strseq(as.character(x[2]), end, end, genome.build)
	  AAseq = strsplit(as.character(x[5]),'')[[1]]
	  NSseq = AAseq[length(AAseq)]
      if( InsSseq != NSseq )
      {
        mhr1 = 0
      }	  	    
      while(mhr1 && mhr && end>=start)
      {
        start = start + 1
        repnumr = end - start + 1
        repseq = strseq(as.character(x[2]), start, end, genome.build)
        if(grepl(repseq, as.character(x[5])))
		{
			seqreversetmp = insf(AAseq, length(AAseq)-repnumr+1, length(AAseq))
			if(repseq == seqreversetmp)
			{
			  mhr = 0
			}
		}
      }	  
      if( repnuml>=repnumr )
      {
        if(mhl==0)
        {
          repnum = repnuml				
        }
      }else{
        if(mhr==0)
        {
          repnum = repnumr				
        }
      }
      if(repnum>0)
      {
        type2 = 'M'						
      }	  		
	}
  }
  if(nlength>5)
  {
    nlength = 5
  }
  if(repnum>5)
  {	
    repnum = 5
  }			
  ty = paste(nlength,":",type1,":",type2,":",repnum,sep="")	
  if( type1 == 'Del' )
  {
    if(seqreverse != refseq )
    {
      ty = 'non_matching'
    }	
  }else{
	if(type2 %in% 'NA')
	{
		ty = NA
	}
  }
  return(ty)
}
