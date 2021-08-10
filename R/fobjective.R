fobjective <- function(xin, k, n, m, X, beta, alpha)
{
  xin = abs(xin)
  P0 = matrix(xin[1:(n*k)], n, k)  # n*k
  S0 = matrix(xin[-c(1:(n*k))], k, m)  # k*m
  res = 1/2 * sum( (X - P0%*%S0)^2 ) + alpha * ( sum(crossprod(P0)) - sum(P0^2) ) + beta * sum(S0)
  res
}
