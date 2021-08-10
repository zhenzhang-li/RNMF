Box_pos = function(x_n=4,y_n=4,line_length=0.9,jap=0.1,alpha=60,beta=30,data2=NULL,type=FALSE,top=FALSE)
{
  data=array(NA,dim=c(2*x_n,2*y_n,2))
  k1=tan(-pi*alpha/180)
  k2=tan(pi*beta/180)
  if(isTRUE(type)){
    if(!is.null(data2)){
      dimm1=dim(data2)
      if(isTRUE(top)){
        pos_start=data2[16,24,]
        for(i in 1:(2*x_n)){
          if(i%%2==1){
            sum=0.6+(i-1)/2*line_length+(i-1)/2*jap
            data[i,2*y_n,]=c(pos_start[1]+sum*cos(pi*alpha/180),pos_start[2]-sum*sin(pi*alpha/180))
          }else{
            sum1=0.6+i/2*line_length+(i/2-1)*jap
            data[i,2*y_n,]=c(pos_start[1]+sum1*cos(pi*alpha/180),pos_start[2]-sum1*sin(pi*alpha/180))        
          }
        }      
        data3=array(NA,dim=c(1,2*y_n,2))
        pos_start1=data2[1,24,]
        for(i in 1:(2*y_n)){
          if(i%%2==1){
            sum=(i-1)/2*(line_length+jap)
            data3[1,2*y_n-i+1,]=c(pos_start1[1]-sum*cos(pi*beta/180),pos_start1[2]-sum*sin(pi*beta/180))
          }else{
            sum1=i/2*line_length+(i/2-1)*jap
            data3[1,2*y_n-i+1,]=c(pos_start1[1]-sum1*cos(pi*beta/180),pos_start1[2]-sum1*sin(pi*beta/180))        
          }
        } 
        for(j in (2*y_n-1):1){    
          data[1,j,]=posfun(data[1,2*y_n,],data3[1,j,],k1,k2) 
        }          
        for(j in (2*y_n-1):1){
          for(i in 2:(2*x_n)){
            data[i,j,]=posfun(data[i,j+1,],data[i-1,j,],k1,k2) 
          }
        }
      }else{   
        pos_start=data2[16,1,] 
        for(i in 1:(2*x_n))
		{
          if(i%%2==1)
		  {
            #        sum=pos_start[1]/cos(pi*alpha/180)+0.6+(i-1)/2*line_length+(i+1)/2*jap
            sum=pos_start[1]/cos(pi*alpha/180)+0.6+(i-1)/2*line_length+(i-1)/2*jap
            data[i,1,]=c(sum*cos(pi*alpha/180),-sum*sin(pi*alpha/180))
          }else{
            #        sum1=pos_start[1]/cos(pi*alpha/180)+0.6+i/2*line_length+i/2*jap
            sum1=pos_start[1]/cos(pi*alpha/180)+0.6+i/2*line_length+(i/2-1)*jap
            data[i,1,]=c(sum1*cos(pi*alpha/180),-sum1*sin(pi*alpha/180))        
          }
        }
        data3=array(NA,dim=c(1,2*y_n,2))
        for(i in 1:(2*y_n))
		{
          if(i%%2==1){
            sum=(i-1)/2*(line_length+jap)
            data3[1,i,]=c(sum*cos(pi*beta/180),sum*sin(pi*beta/180))
          }else{
            sum1=i/2*line_length+(i/2-1)*jap
            data3[1,i,]=c(sum1*cos(pi*beta/180),sum1*sin(pi*beta/180))        
          }
        } 
        for(j in 2:(2*y_n))
		{    
          data[1,j,]=posfun(data[1,1,],data3[1,j,],k1,k2) 
        }          
        for(j in 2:(2*y_n))
		{
          for(i in 2:(2*x_n))
		  {
            data[i,j,]=posfun(data[i,j-1,],data[i-1,j,],k1,k2) 
          }
        } 
      }
    }else{
      for(i in 1:(2*x_n))
	  {
        if(i%%2==1){
          sum=(i-1)/2*(line_length+jap)
          data[i,1,]=c(sum*cos(pi*alpha/180),-sum*sin(pi*alpha/180))
        }else{
          sum1=i/2*line_length+(i/2-1)*jap
          data[i,1,]=c(sum1*cos(pi*alpha/180),-sum1*sin(pi*alpha/180))        
        }
      }  
      for(i in 1:(2*y_n))
	  {
        if(i%%2==1){
          sum=(i-1)/2*(line_length+jap)
          data[1,i,]=c(sum*cos(pi*beta/180),sum*sin(pi*beta/180))
        }else{
          sum1=i/2*line_length+(i/2-1)*jap
          data[1,i,]=c(sum1*cos(pi*beta/180),sum1*sin(pi*beta/180))        
        }
      } 
      for(j in 2:(2*y_n))
	  {
        for(i in 2:(2*x_n)){
          data[i,j,]=posfun(data[i,j-1,],data[i-1,j,],k1,k2) 
        }
      }
    }      
  }else{
    for(i in 1:(2*x_n))
	{
      if(i%%2==1){
        sum=(i-1)/2*(line_length+jap)
        data[i,1,]=c(sum*cos(pi*alpha/180),-sum*sin(pi*alpha/180))
      }else{
        sum1=i/2*line_length+(i/2-1)*jap
        data[i,1,]=c(sum1*cos(pi*alpha/180),-sum1*sin(pi*alpha/180))        
      }
    }  
    for(i in 1:(2*y_n))
	{
      if(i%%2==1){
        sum=(i-1)/2*(line_length+jap)
        data[1,i,]=c(sum*cos(pi*beta/180),sum*sin(pi*beta/180))
      }else{
        sum1=i/2*line_length+(i/2-1)*jap
        data[1,i,]=c(sum1*cos(pi*beta/180),sum1*sin(pi*beta/180))        
      }
    } 
    for(j in 2:(2*y_n)){
      for(i in 2:(2*x_n)){
        data[i,j,]=posfun(data[i,j-1,],data[i-1,j,],k1,k2) 
      }
    }
  }
  return(data)
}