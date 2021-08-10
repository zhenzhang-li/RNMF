PlotSBS96 = function(file)
{
  library( data.table )
  data = as.data.frame(fread(file))
  data(SBsRaw)
  rownames(data) = as.character(data[,1])
  data = data[SBS[,1],]  
  signum = dim(data)[2] - 1
  Name = data[,1]
  sss = strsplit(Name,'')
  Name = NULL
  for(i in 1:96)
  {
    Name = c(Name, paste(sss[[i]][1],sss[[i]][3],sss[[i]][7],sep=''))    
  }
  Name = strsplit(Name,'')
  Ncolnames = colnames(data)[-1]
  #color = c('#426ab3','#87481f','#000000','#f47920','#b2d235','#f58f98')
  #color = c('#8B0016','#367517','#000000','#f47920','#33a3dc','#f58f98')
  color = c('#33a3dc','#000000','#d71345','#d3c6a6','#489620','#F6B297') 
  typeseq = c('C>A','C>G','C>T','T>A','T>C','T>G')
  colorseq = NULL
  for( i in 1:6 )
  {
    colorseq = c(colorseq,rep(color[i],16))
  }  
  #
  for( i in 1:signum )
  {
    pdf(paste(Ncolnames[i], ".pdf", sep=""),18,6)
    layout(matrix(data=c(1,2,3), nrow=3, ncol=1), widths=c(10), heights=c(1,4,1))
    #layout.show(3)
    xstep = 1
    ########### 1 ################
    par(mar=c(0,6,3,1))
    plot(0,0,type="l",axes=F,xlab="",ylab="",xlim=c(0,xstep*96+97),ylim=c(0,0.07))
    for(j in 1:6)
    {
      rect(1+16*(j-1)+16*(j-1)*xstep,0,16*(j+1)+16*(j-1)*xstep,0.03,col=color[j],border=NA)
      text((1+16*(j-1)+16*(j-1)*xstep+16*(j+1)+16*(j-1)*xstep)/2,0.05,labels=typeseq[j],font=2,cex=3)
    }  
      
    ########### 2 ################    
    par(mar=c(0,6,0,1))
    x = as.numeric(as.character(data[,i+1]))
    ymax = round(max(x),3)*7/5
    ymax = min(1,5*round(ymax/5,3))
    ymax = max(ymax, 0.01)	
    # barplot(x, col=colorseq, ylim=c(0,ymax))
    plot(0,0,type="l",axes=F,xlab="",ylab="",xlim=c(0,xstep*96+97),ylim=c(0,ymax))
    for(j in 1:96)
    {
      rect(j+(j-1)*xstep,0,j+1+(j-1)*xstep,x[j],col=colorseq[j],border=NA)
    }  
    segments(0,0,xstep*96+97,0,col='#272727')
    segments(0,0,0,ymax,col='#272727')
    segments(0,ymax,xstep*96+97,ymax,col='#272727')
    segments(xstep*96+97,0,xstep*96+97,ymax,col='#272727')
    text(5,ymax-ymax/10,labels=Ncolnames[i],cex=2,font=2,col='black')
    axis(2,at=seq(0,ymax,ymax/5),labels=paste(seq(0,ymax,ymax/5)*100,"%",sep=""),col="black",col.axis="black",lwd=1,lty=1,las=1,cex.axis=1.5,line=-3.5)
    mtext("Percentage of Single Base Substitutions",side=2,font=1,cex=1.3, padj = -1.5, adj = NA)
    seqseq = NULL
    for(j in 1:96)
    {
      seqseq = c(seqseq, j+(j-1)*xstep+0.5)
    }        
    axis(1,at=seqseq,labels=NA,col="#272727",col.axis="#272727",lwd=1,lty=1,las=2,cex.axis=0.8,line=-0.75)
    
    ########### 3 ################    
    par(mar=c(2,6,0,1))   
    ymax = 0.05
    grap = 0.02    
    plot(0,0,type="l",axes=F,xlab="",ylab="",xlim=c(0,xstep*96+97),ylim=c(0,ymax+grap))
    for(ij in 1:3 )
    {
      for(j in 1:96)
      {
        if(ij==2)
        {
          text(seqseq[j],ymax-grap*(ij-1),labels=Name[[j]][ij],col=colorseq[j],font=1.5,cex=2)          
        }else{
          text(seqseq[j],ymax-grap*(ij-1),labels=Name[[j]][ij],col='grey80',font=1,cex=2)         
        }
      }        
    }
    dev.off()
  }
}
