posfun=function(p1,p2,k1,k2)
{
  kpos=NULL
  x0=(p2[2]-p1[2]+k2*p1[1]-k1*p2[1])/(k2-k1)
  y0=p2[2]-k1*p2[1]+k1*x0
  kpos=c(x0,y0)
  # cat("kpos",kpos,"\n")
  return(kpos)
}