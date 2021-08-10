polygon_color = function(A=NULL,B=NULL,C=NULL,D=NULL,col=NULL)
{
  if(!is.null(col)){
    x=c(A[1],B[1],C[1],D[1])
    y=c(A[2],B[2],C[2],D[2])
    polygon(x,y,col=col,border=NULL)
  }
}
