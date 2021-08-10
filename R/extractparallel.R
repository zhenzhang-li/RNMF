extractparallel = function(i, genomes, totalIterations, totalProcesses, alpha, beta, eps)
{
  removeWeakMutationTypes = 0.001   
  index = apply(genomes, 1, sum)/sum(genomes) < removeWeakMutationTypes
  genomes[index,] = 0  
  totalMutationTypes = nrow(genomes)
  totalGenomes = ncol(genomes)
  dirichletGenomes = DirichletCancerGenomes( genomes ) 
  dirichletGenomes[is.nan(dirichletGenomes)] = 0
  dirichletGenomes[is.na(dirichletGenomes)] = 0  
  dirichletGenomes[is.infinite(dirichletGenomes)] = 10^abs(log10(eps))
  dirichletGenomes[dirichletGenomes<eps] = eps
  dirichletGenomes = normalize(dirichletGenomes)    
  res = nmf(dirichletGenomes, totalProcesses, eps, alpha, beta) 
  E = dirichletGenomes - res$P%*%res$S
  return(list('P'=res$P,'S'=res$S,'E'=E))
}