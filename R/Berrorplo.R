Berrorplo = function( X, data, sampleNames, subtypes, AnalCOSMICSigType )
{
	P = data$P
	S = data$S
	k = data$k
	ks = sort(unique(k))
	datmp = array(0, dim=c(length(P), length(ks)))	
	for(i in 1:length(ks))
	{	
		ix = which( k == ks[i] )
		if(length(ix)>0)
		{
			ctmp = NULL
			for(j in 1:length(ix))
			{
				ctmp = c(ctmp, sum((X-P[[ix[j]]]%*%S[[ix[j]]])^2))		
			}
			datmp[, i] = ctmp
			indexmink = which.min(ctmp)
			ptxt = P[[ix[indexmink]]]
			colN = c(paste(AnalCOSMICSigType, 1:ks[i], sep=""))
			colnames(ptxt) = NULL
			ptxt = as.data.frame(rbind(colN, ptxt))
			rownames(ptxt) = c('Subtype', subtypes$V1)
			write.table(ptxt, file=paste('SignatureComposition.', ks[i], '.Normalized.txt', sep=""), quote=F, col.names=FALSE, row.names=TRUE, sep="\t")
			stxt = S[[ix[indexmink]]]
			stxt = as.data.frame(rbind(sampleNames$V1, stxt))
			rowN = c('SampleID', paste(AnalCOSMICSigType, 1:ks[i], sep=""))
			rownames(stxt) = rowN			
			write.table(stxt, file=paste('SampleContribution.', ks[i], '.Normalized.txt', sep=""), quote=F, col.names=FALSE, row.names=TRUE, sep="\t")
		}
	}
	cfin95dataset = apply(datmp, 2, confidencexy)
    	sda <- unlist(sapply(cfin95dataset, "[[", "mean"))
    	colorset =  rgb(186, 84, 255, 50, maxColorValue=255)  	
	pdf("NormalizedError.pdf", 6, 6)	
	plot(sda, type='l', ylim=c(0, 1.12*max(sda)), xlim=c(min(k),max(k)), ylab='Normalized Error', xlab='Rank', lwd=0.2, col='black') 
	points(min(k):max(k), sda, pch = 20, col='black', cex=0.8)
	dev.off()	
	cfin95dataset = apply(diff(t(datmp)), 1, confidencexy)
    	dsda <- unlist(sapply(cfin95dataset, "[[", "mean"))
	pdf("StabilityofNormalizedError.pdf", 6, 6)
	plot(ks[-1], dsda, type='l', ylab='Stability of Normalized Error', xlab='Rank', lwd=0.2, col='black')
	points(ks[-1], dsda, pch = 20, col='black', cex=0.8)
	abline(h=0, col='red', lwd = 1.5)		
	dev.off()	
}
