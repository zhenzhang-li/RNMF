cal_line_d = function(P1=NULL,P2=NULL,P3=NULL,P4=NULL)
{
  if(is.null(P1) || is.null(P3)){
    d=0
  }else{
    if(P2[1]-P1[1]==0 || P3[1]-P4[1]==0){
      d=abs(P4[1]-P1[1])
    }else{
      k=(P2[2]-P1[2])/(P2[1]-P1[1])
      b1=P1[2]-k*P1[1]
      b2=P3[2]-k*P3[1]
      d=abs(b1-b2)/sqrt(1+k^2)
    }
  } 
  return(d)
}