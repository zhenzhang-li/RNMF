fitcolor = function(A=NULL,B=NULL,C=NULL,D=NULL,col=NULL)
{
  if(!is.null(col)){
    x=c(A[1],B[1],C[1],D[1])
    y=c(A[2],B[2],C[2],D[2])
    x_lim=seq(min(x),max(x),by=0.01)
    y_lim=seq(min(y),max(y),by=0.01)
    for(i in x_lim){
      for(j in y_lim){
        P=c(i,j)
        ###for A,B ###
        d1=cal_d(P,P1=A,P2=B)
        ###for B,C ###
        d2=cal_d(P,P1=B,P2=C)
        ###for C,D ###
        d3=cal_d(P,P1=C,P2=D)
        ###for D,A ###
        d4=cal_d(P,P1=D,P2=A)
        ## A,B  C,D ##
        dd1=cal_line_d(A,B,C,D)
        ## B,C  A,D ##
        dd2=cal_line_d(B,C,D,A)
        if(d1+d3<=dd1 && d2+d4<=dd2){
          points(i,j,col=col,pch=16,cex=0.01)
        }
      } 
    }
  }
}
