ScanCBSSimPlot <-
function(cases, controls, CBSObj, trueTau, SpikeMat, filename, mainTitle, length.out=1000, smoothF=25, xlabScale=10^6, width=12, height=9) {
	maxCase = max(cases)
	maxControl = max(controls)
	maxVal = max(c(maxCase, maxControl))
	cpts = matrix(CBSObj$statHat[,c(1,2,5)], nrow=nrow(CBSObj$statHat))
	tauHat = CBSObj$tauHat
	tauHat = tauHat[c(-1, -length(tauHat))]
	relCN = CBSObj$relCN
	ylims = c(min(c(0, cpts[,3])), max(c(0, cpts[,3])))
	relCNlims = c(min(relCN), max(relCN))
	
	grid.fix = seq(1, maxVal, length.out=length.out)
	gridSize = grid.fix[2]-grid.fix[1]
	casesCountInGrid = getCountsInWindow(cases, 0, maxVal, gridSize, sorted=FALSE)
	casesCountInGridSmooth = lowess(x=grid.fix, y=casesCountInGrid, smoothF/length.out)
	casesCountInGridSmooth$y[casesCountInGridSmooth$y<0] = 0
	controlCountInGrid = getCountsInWindow(controls, 0, maxVal, gridSize, sorted=FALSE)
	controlCountInGridSmooth = lowess(x=grid.fix, y=controlCountInGrid, smoothF/length.out)
	controlCountInGridSmooth$y[controlCountInGridSmooth$y<0] = 0
	PInGrid = casesCountInGrid/(casesCountInGrid+controlCountInGrid)
	PInGrid[is.nan(PInGrid)]=0
	PInGridSmooth = lowess(x=grid.fix, y=PInGrid, smoothF/length.out)
	tauHatInGrid = grid.fix[tauHat %/% gridSize]/xlabScale
	trueTauInGrid = grid.fix[trueTau %/% gridSize]/xlabScale

	pdf(paste(filename, ".pdf", sep=""), width=width, height=height)
	par(mfrow=c(2,1))
	plot(x=grid.fix/xlabScale, y=rep(0, length(grid.fix)), type="n", ylim=ylims, main=mainTitle, ylab="Statistic", xlab=paste("Base Pairs", xlabScale))
	for(i in 1:nrow(cpts)) {
		plotX = c(grid.fix[max(floor(cpts[i,1]/gridSize), 1)]/xlabScale, grid.fix[ceiling(cpts[i,2]/gridSize)]/xlabScale)
		lines(x=plotX, y=rep(cpts[i,3],2), lwd=2)
	}
	for(i in 1:nrow(SpikeMat)) {
		plotX = c(grid.fix[max(floor(SpikeMat[i,3]/gridSize), 1)]/xlabScale, grid.fix[ceiling(SpikeMat[i,4]/gridSize)]/xlabScale)
		lines(x=plotX, y=rep(SpikeMat[i,5],2), lwd=2, col=2)
	}
	abline(v=tauHatInGrid, lty=3, col=2)
	abline(v=trueTauInGrid, lty=3, col=3)
	
	matplot(x=grid.fix/xlabScale, y=log(cbind(casesCountInGridSmooth$y, controlCountInGridSmooth$y)+1), type="l", lty=c(1,1), main="Smoothed Log Read Intensity", ylab="Smoothed Intensity, P", xlab=paste("Base Pairs", xlabScale))
	abline(v=tauHatInGrid, lty=3, col=2)
	abline(v=trueTauInGrid, lty=3, col=3)
	legend("topright", c("case","control"), lty=c(1,1), col=1:2)
	dev.off()
}

