denovoNMF = function(originalGenomes, sampleNames, subtypes, AnalCOSMICSigType = 'SBS', kmin = 1, kmax = 10, steptol = 10^-9, totalIterations = 20, spacetime = 100, mc.cores = 1)
{
  library( data.table )
  library( MCMCpack )
  library( parallel )  
  originalGenomes = as.data.frame(fread(originalGenomes, header=F))
  sampleNames = as.data.frame(fread(sampleNames, header=F))
  subtypes = as.data.frame(fread(subtypes, header=F))
  Y = as.matrix(originalGenomes)
  X = normalize(Y)  
  error = steptol
  myexp = 10^-9
  alpha = 17.6 
  beta = 0.001 
  totalIterations = totalIterations
  if(totalIterations < 5)
  {
    totalIterations = 5
	cat("The parameter totalIterations setting is less than 5, so reset to automatically select the minimum threshold of 5 ... \n")
  }     
  if(spacetime < 20)
  {
    spacetime = 20
	cat("The parameter spacetime setting is less than 20, so reset to automatically select the minimum threshold of 20 ... \n")
  }    
  spacetime = spacetime  
  eps = .Machine$double.eps
  kseqs = kmin:kmax
  allres = NULL
  try(load("RESULT.RData"),silent=TRUE)
  if(is.null(allres))
  {	
	  for( kst in 1:length(kseqs) )
	  {
		totalProcesses = kseqs[kst]
		for( ist in 1:spacetime )
		{
		  Pall = NULL
		  Sall = NULL
		  Eall = NULL
		  xtmp = 1:totalIterations
		  Es = rep(0, totalIterations)	
		  istepNUM = 1
		  while(istepNUM)
		  {
			res = mclapply(xtmp, extractparallel, Y, totalIterations, totalProcesses, alpha, beta, eps, mc.cores = mc.cores) 
			for( i in xtmp )
			{ 
				Pall[[i]] = res[[i]]$P
				Sall[[i]] = res[[i]]$S
				Es[i] = sum((res[[i]]$E)^2)   
			}
			index = sort.list(Es)[1:5]
			Palltmp = NULL
			Salltmp = NULL
			for(i in 1:length(index))
			{
				Palltmp = cbind(Palltmp, Pall[[index[[i]]]])
				Salltmp = rbind(Salltmp, Sall[[index[[i]]]])  
			}		  
			if(is.matrix(Palltmp) && is.matrix(Salltmp))
			{
				istepNUM = 0
			}
		  }
		  
		  P = t(kmeans(t(Palltmp), totalProcesses)$centers)
		  S = kmeans(Salltmp, totalProcesses)$centers
		  m = ncol(S)
		  n = nrow(P)
		  xin = c(as.vector(P), as.vector(S))
		  res = nlm(fobjective, xin, totalProcesses, n, m, X, beta, alpha, gradtol = 1e-4, stepmax = 1000, steptol = 1e-4, iterlim = 100)
		  xin = abs(res$estimate)
		  xin[xin<eps] = eps
		  P = matrix(xin[1:(n*totalProcesses)], n, totalProcesses)  
		  S = matrix(xin[-c(1:(n*totalProcesses))], totalProcesses, m)  
		  P = normalize(P) 
		  S = normalize(S) 
		  res = IterateInitialvalue( X, totalProcesses, S, P, error, myexp)		  
		  allres$P[[ist + spacetime*(kst-1)]] = res$P
		  allres$S[[ist + spacetime*(kst-1)]] = res$S
		  allres$E[[ist + spacetime*(kst-1)]] = sum((X - res$P%*%res$S)^2)
		  allres$k[ist + spacetime*(kst-1)] = totalProcesses   
		  cat("Category:", kst, "; Imitative space:", ist, " ...\n")
		}
	  }
	  save(allres, file = "RESULT.RData")
  }  
  k = allres$k 
  if(spacetime < 100)
  {
     seqs = 5:(spacetime-5)
  }else{
     seqs = 5:95 
  }  
  datmp = array(0, dim=c(length(seqs), length(c(min(k):max(k)))))
  for( i in 1:length(seqs) )
  {
     restmp = S0choose( X, allres, tail = seqs[i] )
     datmp[i, ] = SilhouetteCoefficient(restmp)		
  }	
  cfin95dataset = apply(datmp, 2, confidencexy)
  processStabAvg <- unlist(sapply(cfin95dataset, "[[", "mean"))
  cd <- unlist(sapply(cfin95dataset, "[[", "lower"))
  cu <- unlist(sapply(cfin95dataset, "[[", "upper"))  
  y = c(cd, rev(cu))
  aseq = min(k):max(k)
  xseq = c(aseq, rev(aseq))
  colorset =  rgb(186, 84, 255, 50, maxColorValue=255)  
  pdf("StabilityEvaluation.pdf", 6, 6)
  plot(processStabAvg, type='l', ylim=c(0, 1.12), xlim=c(min(k),max(k)), ylab='Stability Evaluation', xlab='Rank', lwd=1, col='black')
  polygon(xseq, y, col=colorset, border=NA)
  legend("topright", legend=c("Mean","95% confidence interval"), col=c("black",colorset), lty=1, lwd=3 ) 
  points(min(k):max(k), processStabAvg, pch = 20, col='black', cex=0.6)
  text(min(k):max(k), processStabAvg+0.05, round(processStabAvg,2), col='red')  
  dev.off()	  
  errorstabledata = as.data.frame(cbind(aseq, processStabAvg, cd, cu))
  colnames(errorstabledata) = c('Rank', 'Mean', 'Left of 95%-CI', 'Right of 95%-CI')
  write.table(errorstabledata, file='StabilityEvaluation.txt', quote=F, col.names=TRUE, row.names=FALSE, sep="\t")
  Berrorplo( X, allres, sampleNames, subtypes, AnalCOSMICSigType = AnalCOSMICSigType )
}
