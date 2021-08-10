Box_plot = function(data=NULL,z_high=NULL,col=NULL,type=FALSE,method=NULL,theta=0,border=border){
  if(!is.null(data)){
    dimm=dim(data)
    if(!is.null(z_high)){ 
      num=seq(1,7,by=2)
      num2=seq(9,15,by=2)
      num1=seq(23,1,by=-2)
      if(is.null(method)){
        polyplot(num1,num,data,z_high,col,theta,border=border) 
        polyplot(num1,num2,data,z_high,col,theta,border=border)        
      }else{
        polyplot(num1,num,data,z_high,col,method,theta,border=border) 
        polyplot(num1,num2,data,z_high,col,method,theta,border=border)
      }
    }else{
      if(!isTRUE(type)){
        num=seq(1,15,by=2)
        num1=seq(1,23,by=2)
        for(j in num1){
          for(i in num){
            A=data[i,j,]
            B=data[i+1,j,]
            D=data[i,j+1,]
            C=data[i+1,j+1,]
            segments(A[1],A[2],B[1],B[2])
            segments(A[1],A[2],D[1],D[2])
            segments(C[1],C[2],B[1],B[2])
            segments(C[1],C[2],D[1],D[2])
          }
        } 
      }else{
        num=c(1,3,5,7)
        for(j in num){
          for(i in num){
            A=data[i,j,]
            B=data[i+1,j,]
            D=data[i,j+1,]
            C=data[i+1,j+1,]
            pos_p=pos_get(A,B,C,D,z_high=0.3,theta)
            A1=pos_p$A1
            B1=pos_p$B1
            C1=pos_p$C1
            D1=pos_p$D1          
            segments(A[1],A[2],B[1],B[2])
            segments(C[1],C[2],B[1],B[2])
            segments(A[1],A[2],A1[1],A1[2])
            segments(B[1],B[2],B1[1],B1[2])
            segments(C[1],C[2],C1[1],C1[2])
            segments(A1[1],A1[2],D1[1],D1[2])
            segments(A1[1],A1[2],B1[1],B1[2])
            segments(B1[1],B1[2],C1[1],C1[2])
            segments(C1[1],C1[2],D1[1],D1[2])
          }
        } 
      }
    }
  }
}
