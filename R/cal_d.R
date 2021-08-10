cal_d = function(P,P1=NULL,P2=NULL)
{
  if(is.null(P1) || is.null(P2)){
    d=0
  }else{
    A=P2[2]-P1[2]
    B=P1[1]-P2[1]
    C=(P2[1]-P1[1])*P1[2]-P1[1]*(P2[2]-P1[2])
    d=abs(A*P[1]+B*P[2]+C)/sqrt(A^2+B^2)
  } 
  return(d)
}