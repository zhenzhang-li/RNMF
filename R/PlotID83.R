PlotID83 = function(file)
{
  library( data.table )
  data = as.data.frame(fread(file))
  data(IDsRaw)
  rownames(data) = as.character(data[,1])
  data = data[ID[,1],]
  signum = dim(data)[2] - 1
  data = data[1:83,]
  Name = c('1','2','3','4','5','6+',
		   '1','2','3','4','5','6+',
		   '0','1','2','3','4','5+',
		   '0','1','2','3','4','5+',
		   '1','2','3','4','5','6+',
		   '1','2','3','4','5','6+',
		   '1','2','3','4','5','6+',
		   '1','2','3','4','5','6+',
		   '0','1','2','3','4','5+',
		   '0','1','2','3','4','5+',
		   '0','1','2','3','4','5+',
		   '0','1','2','3','4','5+',
		   '1',
		   '1','2',
		   '1','2','3',
		   '1','2','3','4','5+')
  Name_bottom = c('Homopolymer Length','Homopolymer Length','Number of Repeat Units','Number of Repeat Units','Mircohomology Length')
  Name_above = c('>1bp Deletion at Repeats','>1bp Insertions at Repeats','Mircohomology')
  Name_middle = c('1bp Deletion','1bp Insertion','(Deletion Length)','(Insertion Length)','(Deletion Length)')
  typeseq = c('C','T','C','T','2','3','4','5+','2','3','4','5+','2','3','4','5+')
  Ncolnames = colnames(data)[-1]
  
  color = c(
	  '#FAA755',
	  '#F48420',
	  '#A3CF62',
	  '#1D953F',
	  '#F8ABA6',
	  '#F05B72',
	  '#D71345',
	  '#AA363D',
	  '#AFB4DB',
	  '#6F60AA',
	  '#694D9F',
	  '#494E8F',
	  '#D2A6C7',
	  '#AF4A92',
	  '#8F006D',
	  '#64004B') 
	  
  seqnum = c(6,6,6,6,6,6,6,6,6,6,6,6,1,2,3,5)
  seqnum1 = c(1, cumsum(seqnum))
  colorseq = NULL
  for( i in 1:16 )
  {
    colorseq = c(colorseq,rep(color[i], seqnum[i]))
  }  
  colorseq1 = c('black','white','black','white','black','white','white','white','black','white','white','white','black','white','white','white')
  seqnum2 = c(1,12,24,48,72,83)
  seqnum3 = c(24,48,72,83)  
  for( i in 1:signum )
  {
    pdf(paste(Ncolnames[i], ".pdf", sep=""),20,6)
    layout(matrix(data=c(1,2,3), nrow=3, ncol=1), widths=c(10), heights=c(1,4,1))
    #layout.show(3)
    xstep = 1
    ########### 1 ################
    par(mar=c(0,6,3,1))
    plot(0,0,type="l",axes=F,xlab="",ylab="",xlim=c(0,xstep*83+84),ylim=c(0,0.08))
    for(j in 1:(length(seqnum1)-1))
    {
      if(j==1)
      {
        rect(seqnum1[j]+(seqnum1[j]-1)*xstep,0,seqnum1[j+1]+seqnum1[j+1]*xstep,0.03,col=color[j],border=NA)	  
        text((seqnum1[j]+(seqnum1[j]-1)*xstep+seqnum1[j+1]+seqnum1[j+1]*xstep)/2,0.016,labels=typeseq[j],font=1,cex=1.5,col=colorseq1[1])        
      }else{
        rect(seqnum1[j]+(seqnum1[j]+1)*xstep,0,seqnum1[j+1]+seqnum1[j+1]*xstep,0.03,col=color[j],border=NA)	
        text((seqnum1[j]+(seqnum1[j]+1)*xstep+seqnum1[j+1]+seqnum1[j+1]*xstep)/2,0.016,labels=typeseq[j],font=1,cex=1.5,col=colorseq1[j])        
      }
    }  	
    for(j in 1:(length(seqnum2)-1))
    {
      if(j==1)
      {  
        text((seqnum2[j]+(seqnum2[j]-1)*xstep+seqnum2[j+1]+seqnum2[j+1]*xstep)/2,0.045,labels=Name_middle[j],font=1,cex=1.5)        
      }else{
        text((seqnum2[j]+(seqnum2[j]+1)*xstep+seqnum2[j+1]+seqnum2[j+1]*xstep)/2,0.045,labels=Name_middle[j],font=1,cex=1.5)        
      }
    }  	
    for(j in 1:(length(seqnum3)-1))
    {
      if(j==1)
      {  
        text((seqnum3[j]+(seqnum3[j]-1)*xstep+seqnum3[j+1]+seqnum3[j+1]*xstep)/2,0.065,labels=Name_above[j],font=1,cex=1.5)        
      }else{
        text((seqnum3[j]+(seqnum3[j]+1)*xstep+seqnum3[j+1]+seqnum3[j+1]*xstep)/2,0.065,labels=Name_above[j],font=1,cex=1.5)        
      }
    }  	

    ########### 2 ################    
    par(mar=c(0,6,0,1))
    x = as.numeric(as.character(data[,i+1]))
    ymax = round(max(x),3)*7/5
    ymax = min(1,5*round(ymax/5,3))
    ymax = max(ymax, 0.01)		
    # barplot(x, col=colorseq, ylim=c(0,ymax))
    plot(0,0,type="l",axes=F,xlab="",ylab="",xlim=c(0,xstep*83+84),ylim=c(0,ymax))
    for(j in 1:83)
    {
      rect(j+(j-1)*xstep,0,j+1+(j-1)*xstep,x[j],col=colorseq[j],border=NA)
    }
    segments(0,0,xstep*83+84,0,col='#272727')
    segments(0,0,0,ymax,col='#272727')
    segments(0,ymax,xstep*83+84,ymax,col='#272727')
    segments(xstep*83+84,0,xstep*83+84,ymax,col='#272727')
    text(5,ymax-ymax/10,labels=Ncolnames[i],cex=2,font=2,col='black')
    axis(2,at=seq(0,ymax,ymax/5),labels=paste(seq(0,ymax,ymax/5)*100,"%",sep=""),col="black",col.axis="black",lwd=1,lty=1,las=1,cex.axis=1.5,line=-3.5)
    mtext("Percentage of Indels",side=2,font=1,cex=1.3, padj = -1.5, adj = NA)
   
    ########### 3 ################    
    par(mar=c(3,6,0,1))   
    ymax = 0.06
    grap = 0.02  
    plot(0,0,type="l",axes=F,xlab="",ylab="",xlim=c(0,xstep*83+84),ylim=c(0,ymax))
    for(j in 1:(length(seqnum1)-1))
    {
      if(j==1)
      {
        rect(seqnum1[j]+(seqnum1[j]-1)*xstep,ymax,seqnum1[j+1]+seqnum1[j+1]*xstep,ymax-grap,col=color[j],border=NA)	      
      }else{
        rect(seqnum1[j]+(seqnum1[j]+1)*xstep,ymax,seqnum1[j+1]+seqnum1[j+1]*xstep,ymax-grap,col=color[j],border=NA)	      
      }
    }  	
    for(j in 1:83)
    {
	  text(j+(j-1)*xstep+0.5,ymax-1.4*grap,labels=Name[j],cex=1.1,font=1,col='black')
    }
    for(j in 1:(length(seqnum2)-1))
    {
      if(j==1)
      {  
        text((seqnum2[j]+(seqnum2[j]-1)*xstep+seqnum2[j+1]+seqnum2[j+1]*xstep)/2,ymax-2.3*grap,labels=Name_bottom[j],font=1,cex=1.5)        
      }else{
        text((seqnum2[j]+(seqnum2[j]+1)*xstep+seqnum2[j+1]+seqnum2[j+1]*xstep)/2,ymax-2.3*grap,labels=Name_bottom[j],font=1,cex=1.5)        
      }
    }  
    dev.off()
  }
}
