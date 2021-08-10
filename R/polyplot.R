polyplot = function(num1,num,data,z_high,color,method="polygon",theta,border='grep80'){
  if(all(num>8)){
    k=44
  }else{
    k=0
  }
  for(j in num1){
    for(i in num){
      A=data[i,j,]
      D=data[i+1,j,]
      B=data[i,j+1,]
      C=data[i+1,j+1,]
      pos_p=pos_get(A,B,C,D,z_high[k+2*(23-j)+(i+1)/2],theta)
      A1=pos_p$A1
      B1=pos_p$B1
      C1=pos_p$C1
      D1=pos_p$D1
      m=k+2*(23-j)+(i+1)/2      
      #     m=floor(m/17)+1
      if(m < 17){
        m=1
      }else if(m < 33){
        m=2
      }else if(m < 49){
        m=3
      }else if(m < 65){
        m=4
      }else if(m < 81){
        m=5
      }else{
        m=6
      }      
      if(method=="polygon"){
        polygon_color(A,A1,D1,D,color[m])
        polygon_color(A1,B1,C1,D1,col=color[m])
        polygon_color(D,D1,C1,C,col=color[m])    
      }else{
        fitcolor(A,A1,D1,D,color[m])
        fitcolor(A1,B1,C1,D1,col=color[m])
        fitcolor(D,D1,C1,C,col=color[m])
      }
      segments(A[1],A[2],D[1],D[2],col=border,lwd=0.1)
      segments(C[1],C[2],D[1],D[2],col=border,lwd=0.1)
      segments(A[1],A[2],A1[1],A1[2],col=border,lwd=0.1)
      segments(D[1],D[2],D1[1],D1[2],col=border,lwd=0.1)
      segments(C[1],C[2],C1[1],C1[2],col=border,lwd=0.1)
      segments(A1[1],A1[2],D1[1],D1[2],col=border,lwd=0.1)
      segments(A1[1],A1[2],B1[1],B1[2],col=border,lwd=0.1)
      segments(B1[1],B1[2],C1[1],C1[2],col=border,lwd=0.1)
      segments(C1[1],C1[2],D1[1],D1[2],col=border,lwd=0.1)
    }
  }
}
