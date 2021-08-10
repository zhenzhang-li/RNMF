DBSlego = function(file='originalGenomes', subtype='subtypes', genome.build = c("Ch37","Ch38"), SequenceType = c("WGS","WES"))
{
  if(length(genome.build)>1)
  {
    genome.build = 'Ch37'
  }
  if( genome.build %in% 'Ch37' )
  {
    data(Ch37)
  }else
  {
    data(Ch38)
  }  
  if( 'WGS' %in% SequenceType )
  {
      gelength = sum(chrom$LEN)
  }else{
	  gelength = 28*10^6  
  }
  
  library( data.table )
  data = as.data.frame(fread(file))
  subtypecode = as.character(as.data.frame(fread(subtype, header=F, sep="\t"))$V1)
  index = NULL
  data(DBsRaw)
  for(i in 1:length(DBs[,1]))
  {
	ix = which(subtypecode == DBs[,1][i])
	index = c(index, ix)  
  }
  data = data[index, ]
  num = apply(data, 1, sum)/gelength
  data = cbind(DBs, t(t(num)))
  data = as.data.frame(data)
  data = as.data.frame(data)
  data$V2 = as.numeric(data$V2)
  signum = dim(data)[2] - 1
  data = data[1:78,]
  colnames(data) = c('Subtype','DBsMutRate')
  Name = data[,1]
  Name = strsplit(Name,'>')
  Ncolnames = colnames(data)[-1]
  color = c('#09C7F7','#1A42E6','#A9C43C','#367517','#F5A89A','#DF0029','#F5B16D','#BD6B09','#A095C4','#3A2885') 
  seqnum = c(9,6,9,6,9,6,6,9,9,9)
  seqnum1 = c(1, cumsum(seqnum))
  typeseq = c('AC>NN','AT>NN','CC>NN','CG>NN','CT>NN','GC>NN','TA>NN','TC>NN','TG>NN','TT>NN')
  colorseq = NULL
  for( i in 1:10 )
  {
    colorseq = c(colorseq,rep(color[i], seqnum[i]))
  }  
  for( i in 1:signum )
  {
    pdf(paste(Ncolnames[i], ".pdf", sep=""),18,6)
    layout(matrix(data=c(1,2,3), nrow=3, ncol=1), widths=c(10), heights=c(1,4,1))
    #layout.show(3)
    xstep = 1
    ########### 1 ################
    par(mar=c(0,6,3,1))
    plot(0,0,type="l",axes=F,xlab="",ylab="",xlim=c(0,xstep*78+79),ylim=c(0,0.07))
    for(j in 1:(length(seqnum1)-1))
    {
      if(j==1)
      {
        rect(seqnum1[j]+(seqnum1[j]-1)*xstep,0,seqnum1[j+1]+seqnum1[j+1]*xstep,0.03,col=color[j],border=NA)	  
        text((seqnum1[j]+(seqnum1[j]-1)*xstep+seqnum1[j+1]+seqnum1[j+1]*xstep)/2,0.045,labels=typeseq[j],font=2,cex=2)        
      }else{
        rect(seqnum1[j]+(seqnum1[j]+1)*xstep,0,seqnum1[j+1]+seqnum1[j+1]*xstep,0.03,col=color[j],border=NA)	
        text((seqnum1[j]+(seqnum1[j]+1)*xstep+seqnum1[j+1]+seqnum1[j+1]*xstep)/2,0.045,labels=typeseq[j],font=2,cex=2)        
      }
    }  	
    
    ########### 2 ################    
    par(mar=c(0,6,0,1))
    x = as.numeric(as.character(data[,i+1]))
	xlog = floor(abs(log10(max(x))))
    #ymax = 10^(-xlog-1)
	xcn = floor(log(1/max(x)))+1
	ymax = xcn*10^(-xlog-1)	
	
    # barplot(x, col=colorseq, ylim=c(0,ymax))
    plot(0,0,type="l",axes=F,xlab="",ylab="",xlim=c(0,xstep*78+79),ylim=c(0,ymax))
    for(j in 1:78)
    {
      rect(j+(j-1)*xstep,0,j+1+(j-1)*xstep,x[j],col=colorseq[j],border=NA)
    }
    segments(0,0,xstep*78+79,0,col='#272727')
    segments(0,0,0,ymax,col='#272727')
    segments(0,ymax,xstep*78+79,ymax,col='#272727')
    segments(xstep*78+79,0,xstep*78+79,ymax,col='#272727')
    text(7,ymax-ymax/10,labels=Ncolnames[i],cex=2,font=2,col='black')
    axis(2,at=seq(0,ymax,ymax/5),labels=paste(seq(0,ymax,ymax/5)*100,"%",sep=""),col="black",col.axis="black",lwd=1,lty=1,las=1,cex.axis=1.5,line=-3.5)
    mtext("Mutation Rate of Double Base Substitutions",side=2,font=1,cex=1.3, padj = -2.5, adj = NA)
    seqseq = NULL
    for(j in 1:78)
    {
      seqseq = c(seqseq, j+(j-1)*xstep+0.5)
    } 
    axis(1,at=seqseq,labels=NA,col="#272727",col.axis="#272727",lwd=1,lty=1,las=2,cex.axis=0.8,line=-0.75)
    
    ########### 3 ################    
    par(mar=c(2,6,0,1))   
    ymax = 0.05
    grap = 0.02    
    plot(0,0,type="l",axes=F,xlab="",ylab="",xlim=c(0,xstep*78+79),ylim=c(0,ymax+grap))
    for(ij in 1:2)
    {
      for(j in 1:78)
      {
        text(seqseq[j],ymax-grap*(ij-1),labels=strsplit(Name[[j]][2],'')[[1]][ij],col='#826858',font=1,cex=1.5)         
      }  	
    }
    dev.off()
  }
}
