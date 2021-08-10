DirichletCancerGenomes = function(genomes)
{
  m = ncol(genomes)
  n = nrow(genomes)
  summ = apply(genomes, 2, sum) 
  dirichletGenomes = NULL
  for( i in 1:m )
  {
    dirichletGenomes = cbind(dirichletGenomes, t(round(rdirichlet(1, genomes[, i])*summ[i], 0)))    
  }
  return(dirichletGenomes)
}