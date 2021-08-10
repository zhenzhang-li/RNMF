SBSlego = function(file='originalGenomes', subtype='subtypes', scale=10, name=NULL, title='GC', sort=NULL, top=TRUE, color=NULL, border=NULL, genome.build = c("Ch37","Ch38"), SequenceType = c("WGS","WES"), RegionLength = NULL)
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
  if(!is.null(RegionLength))
  {
     gelength = RegionLength
  }
  h = read.table(file, header=F, sep="\t")
  hh = read.table(subtype, header=F, sep="\t")
  hhtmp1 = strsplit(as.character(hh$V1),'')
  Mut96type = NULL
  Mut6type = NULL 
  for(i in 1:96)
  {
	Mut96type = c(Mut96type, paste(hhtmp1[[i]][1],hhtmp1[[i]][3],hhtmp1[[i]][7],sep=""))
	Mut6type = c(Mut6type, paste(hhtmp1[[i]][3],hhtmp1[[i]][4],hhtmp1[[i]][5],sep=""))
  } 
  if(is.null(name))
  {
    name='cancer'
  }
  if(is.null(title))
  {
    title_name=name
  }else{
    title_name=title
  } 
  if(is.null(color))
  {
	  color = c('#FAA755','#F48420','#A3CF62','#1D953F','#F8ABA6','#F05B72')  	  
  }
 if(isTRUE(sort))
  {
    h1=NULL
    h2=NULL
    if(ncol(h)==1)
	{
      h1=rbind(h1, data.frame(Mut96type, Mut6type, h$V1/sum(h$V1)))
    }else{	
      h1=rbind(h1, data.frame(Mut96type, Mut6type, apply(h, 1, sum)/sum(apply(h, 1, sum))))
    }
    colnames(h1) = c("M96","M6","High")
    Mut6type1 = c("C>T","C>A","C>G","T>C","T>G","T>A")
    lable1 = NULL
    lable = NULL
    for(i in 1:6)
	{
      h3=NULL
      ix=h1$M6==Mut6type1[i]
      h3=h1[ix,]
      if(i==1)
	  {
        kh3 = h3[order(h3$High,decreasing = T),]
        colnames(kh3) = c("M96","M6","High")
        Mtype = as.character(kh3$M96)
        for(ij in 1:16)
		{
          lable1[ij] = paste(strsplit(Mtype[ij], "")[[1]][1],"_",strsplit(Mtype[ij], "")[[1]][3],sep="")
        }
        lable = c(lable1[13],lable1[9],lable1[5],lable1[1],lable1[14],lable1[10],lable1[6],lable1[2],lable1[15],lable1[11],lable1[7],lable1[3],lable1[16],lable1[12],lable1[8],lable1[4])
      }
      for(j in 1:16)
	  {
        st = sub("_",strsplit(Mut6type1[i],">")[[1]][1],lable1[j])
        ix = grepl(st,h3$M96)
        h2 = rbind(h2,h3[ix,])
      }
    }
    h1 = h2
    Mut6type = Mut6type1  
    name1 = paste(name,".sort",sep="")
  }else if(is.null(sort))
  {
    h1=NULL
    h2=NULL
    if(ncol(h)==1)
	{
      h1=rbind(h1, data.frame(Mut96type, Mut6type, h$V1/sum(h$V1)))
    }else{	
      h1=rbind(h1, data.frame(Mut96type, Mut6type, apply(h, 1, sum)/sum(apply(h, 1, sum))))
    }			
    colnames(h1) = c("M96","M6","High")
    Mut6type1 = c("C>T","C>A","C>G","T>C","T>G","T>A")
    lable1 = c("T_G","C_G","A_G","G_G","T_A","C_A","A_A","G_A","T_C","C_C","A_C","G_C","T_T","C_T","A_T","G_T")
    lable = c(lable1[13],lable1[9],lable1[5],lable1[1],lable1[14],lable1[10],lable1[6],lable1[2],lable1[15],lable1[11],lable1[7],lable1[3],lable1[16],lable1[12],lable1[8],lable1[4])
    for(i in 1:6)
	{
      h3=NULL
      ix=h1$M6==Mut6type1[i]
      h3=h1[ix,]
      for(j in 1:16){
        st=sub("_",strsplit(Mut6type1[i],">")[[1]][1],lable1[j])
        ix=grepl(st,h3$M96)
        h2=rbind(h2,h3[ix,])
      }
    }
    h1 = h2
    Mut6type=Mut6type1  
    name1=paste(name,".Nonsort",sep="")
  }else{
    h1=NULL	
    if(ncol(h)==1)
	{
      h1=rbind(h1, data.frame(Mut96type, Mut6type, h$V1/sum(h$V1)))
    }else{	
      h1=rbind(h1, data.frame(Mut96type, Mut6type, apply(h, 1, sum)/sum(apply(h, 1, sum))))
    }		
    colnames(h1)=c("M96","M6","High")
    Mut6type=unique(Mut6type)
    lable=c("T_A","G_A","C_A","A_A","T_C","G_C","C_C","A_C","T_G","G_G","C_G","A_G","T_T","G_T","C_T","A_T")
    name1=name
  }
  data_6type=NULL
  for(i in Mut6type)
  {
    tem=NULL
    ix = grepl(i, as.character(hh$V1))
    tem=h[ix,]
    if(ncol(h)==1)
	{
      data_6type[i] = round(sum(tem$V1))
    }else{  
      data_6type[i] = round(sum(apply(tem,1,sum)))
    }
  }
  allmutations = sum(data_6type)
  mbnum = scale
  scale = mbnum*100/allmutations*gelength/10^6
  z_high=h1$High*100
  lable_max=(max(z_high)+scale-max(z_high)%%scale)/scale
  y_high1=lable_max*scale
  y_high=10
  scale1 = scale*y_high/y_high1
  z_high = y_high*z_high/y_high1
  ##########
  data=NULL
  theta=0
  alpha=83
  beta=83
  scale2=4
  jap=0.3
  line_length1=0.7
  line_length2=1.3
  data=Box_pos(8,12,line_length=line_length1,jap=jap,alpha=alpha,beta=beta,type=FALSE)
  data1=Box_pos(4,4,line_length=line_length2,jap=jap,alpha=alpha,beta=beta,data2=data,type=TRUE,top=top)
  x1=max(apply(data1[,,1],1,max))
  x1=max(x1,max(apply(data[,,1],1,max)))
  x2=min(apply(data[,,1],1,min))
  y2=min(apply(data1[,,2],1,min))
  data_mid=Box_pos(8,12,line_length=line_length1/2,jap=jap+line_length1/2,alpha=alpha,beta=beta,type=FALSE,top=TRUE)
  data1_mid=Box_pos(4,4,line_length=line_length2/2,jap=jap+line_length2/2,alpha=alpha,beta=beta,data2=data,type=TRUE,top=top)
  pos=pos_get(data[1,1,],data[1,24,],data[2,24,],data[2,1,],y_high+scale2,theta)
  y_max=pos$B1[2]
  pdf(file=paste(getwd(),"/",name1,".lego.pdf",sep=""),width=15,height=8) 
  layout(matrix(data=c(1,2), nrow=2, ncol=1), widths=c(20), heights=c(6,3))
  # 1 Draw lego  
  par(mar = c(2,1,2,1))
  plot(0,0,type='l',xlim=c(-ceiling(x2)-1,ceiling(x1)+1),ylim=c(y2,y_max),axes=F,xlab='',ylab='')
  text(x2,y_max,title_name,font=2,cex=1.8)
  text(x2+0.55,y_max-8.5,paste("n=",allmutations,sep=""),font=1,cex=1.2)
  segments(data[1,1,1],data[1,1,2],data[1,24,1],data[1,24,2],lwd=2)
  segments(data[1,1,1],data[1,1,2],data[16,1,1],data[16,1,2],lwd=2)
  segments(data[16,1,1],data[16,1,2],data[16,24,1],data[16,24,2],lwd=2)
  segments(data[1,24,1],data[1,24,2],data[16,24,1],data[16,24,2],lwd=2)
  num_x=c(1,seq(2,24,by=2),25)
  num_y=seq(1,lable_max+1,by=1)
  data_ge1=array(NA,dim=c(length(num_x),length(num_y),2))
  for(i in 1:num_y[length(num_y)])
  {
    for(j in num_x){
      if(i==1){
        if(j==1){
          data_ge1[1,1,]=data_mid[1,j,]
        }else if(j==25){
          data_ge1[length(num_x),1,]=data[1,j-1,]
        }else{
          data_ge1[j/2+1,1,]=data_mid[1,j,]
        }
      }else{
        if(j==1){
          kk=pos_get(A=data_ge1[1,i-1,],z_high=scale1,theta=theta)
          data_ge1[1,i,]=kk$A1        
        }else if(j==25){
          kk=pos_get(A=data_ge1[length(num_x),i-1,],z_high=scale1,theta=theta) 
          data_ge1[length(num_x),i,]=kk$A1
        }else{
          kk=pos_get(A=data_ge1[j/2+1,i-1,],z_high=scale1,theta=theta) 
          data_ge1[j/2+1,i,]=kk$A1
        }         
      }
    }
  }
  for(i in 1:13)
  {
    if(i==1){
      lty=1
      lwd=2
    }else{
      lty=2
      lwd=1
    }
    segments(data_ge1[i,1,1],data_ge1[i,1,2],data_ge1[i,length(num_y),1],data_ge1[i,length(num_y),2],lwd=lwd,lty=lty,col="grey65")
  }
  for(i in 1:1:length(num_y))
  {
    segments(data_ge1[1,i,1],data_ge1[1,i,2],data_ge1[14,i,1],data_ge1[14,i,2],lwd=1,lty=2,col="grey65")
  }
  k1=tan(-pi*alpha/180)
  k2=tan(pi*beta/180)
  segments(data_ge1[1,1,1],data_ge1[1,1,2],data_ge1[1,length(num_y),1],data_ge1[1,length(num_y),2],lty=1,lwd=2)
  for(i in 1:length(num_y))
  {
    pos_z=posfun2(data_ge1[1,i,],length=0.2,k=k1,axis='z')
    segments(pos_z[1],pos_z[2],data_ge1[1,i,1],data_ge1[1,i,2],lty=1,lwd=2)
    text(pos_z[1]-0.1,pos_z[2],(i-1)*mbnum,font=2,cex=1.2)
  }
  mtext("Mutations/Mb",side=2,line=0.5,adj=0.38,padj=14.5,font=2,cex=1.3)
  for(i in seq(2,16,by=2)){
    pos_z=posfun2(data_mid[i,1,],length=0.2,k=k2,axis='x')
    segments(pos_z[1],pos_z[2],data_mid[i,1,1],data_mid[i,1,2],lty=1,lwd=2)
  }
  for(i in seq(2,24,by=2)){
    pos_z=posfun2(data_mid[16,i,],length=0.5,k=k1,axis='y')
    segments(pos_z[1],pos_z[2],data_mid[16,i,1],data_mid[16,i,2],lty=1,lwd=2)
  }
  num_x=c(1,seq(2,16,by=2),17)
  num_y=seq(1,lable_max+1,by=1)
  data_ge1=array(NA,dim=c(length(num_x),length(num_y),2))
  for(i in 1:num_y[length(num_y)])
  {
    for(j in num_x){
      if(i==1){
        if(j==1){
          data_ge1[1,1,]=data[1,24,]
        }else if(j==17){
          data_ge1[length(num_x),1,]=data[16,24,]
        }else{
          data_ge1[j/2+1,1,]=posfun(data_mid[j,24,],data_ge1[j/2,1,],k1,k2) 
        }
      }else{
        if(j==1){
          kk=pos_get(A=data_ge1[1,i-1,],z_high=scale1,theta=theta)
          data_ge1[1,i,]=kk$A1        
        }else if(j==17){
          kk=pos_get(A=data_ge1[length(num_x),i-1,],z_high=scale1,theta=theta) 
          data_ge1[length(num_x),i,]=kk$A1
        }else{
          kk=pos_get(A=data_ge1[j/2+1,i-1,],z_high=scale1,theta=theta) 
          data_ge1[j/2+1,i,]=kk$A1
        }         
      }
    }
  }
  for(i in 2:9){
    if(i==1){
      lty=1
      lwd=2
    }else{
      lty=2
      lwd=1
    }
    segments(data_ge1[i,1,1],data_ge1[i,1,2],data_ge1[i,length(num_y),1],data_ge1[i,length(num_y),2],lwd=lwd,lty=lty,col="grey65")
  }
  for(i in 1:length(num_y)){
    segments(data_ge1[1,i,1],data_ge1[1,i,2],data_ge1[10,i,1],data_ge1[10,i,2],lwd=1,lty=2,col="grey65")
  }
  
  Box_plot(data1,z_high=NULL,col=NULL,type=T,method=NULL,theta=theta,border=border)
  Box_plot(data,z_high=z_high,col=color,type=F,method="polygon",theta=theta,border=border)
  dim_text=dim(data1_mid)[1:2]
  for(i in seq(2,dim_text[1],by=2)){
    if(isTRUE(top)){
      num_top=seq(1,dim_text[1]-1,by=2)
    }else{
      num_top=seq(2,dim_text[1],by=2)
    }
    for(j in num_top){
      KK=pos_get(C=data1_mid[i,j,],z_high=0.2,theta=theta)
      if(isTRUE(top)){
        t=2*(i-2)+(j+1)/2
      }else{
        t=2*(i-2)+j/2
      }
      text(KK$C1[1],KK$C1[2],lable[t],cex=1.1,font=1)
    }
  }
  x_0=max(apply(data[,,1],1,max))
  y_0=y_max-5
  geg=1.9
  for(i in 1:6){
    points(x_0+0.05,y_0-i*geg-0.5,pch=22,cex=3,col=color[i],bg=color[i])
    text(x_0+0.12,y_0-i*geg-0.5,labels=Mut6type[i],font=2,cex=1.2,adj=0,srt=0)
  }
 
  # 2 Draw the pie
  par(mar = c(2,2,2,2))
  total <-sum(data_6type)
  value <-c(data_6type[1]/total,data_6type[2]/total,data_6type[3]/total,data_6type[4]/total,data_6type[5]/total,data_6type[6]/total)
  cols <- color
  pie(value,col=cols,labels=paste(colnames(t(value))," ",round(value*100,2),"%",sep=""))
  dev.off()
}
