pos_get = function(A=NULL,B=NULL,C=NULL,D=NULL,z_high,theta)
{
  if(!is.null(A)){A1=c(A[1]-z_high*tan(pi*theta/180),A[2]+z_high)}else{A1=NULL}
  if(!is.null(B)){B1=c(B[1]-z_high*tan(pi*theta/180),B[2]+z_high)}else{B1=NULL}
  if(!is.null(C)){C1=c(C[1]-z_high*tan(pi*theta/180),C[2]+z_high)}else{C1=NULL}
  if(!is.null(D)){D1=c(D[1]-z_high*tan(pi*theta/180),D[2]+z_high)}else{D1=NULL}
  return(list('A1'=A1,'B1'=B1,'C1'=C1,'D1'=D1))
}