ScanCBSPlot <-
function(cases, controls, CBSObj, filename, mainTitle, CIObj=NULL, length.out=10000, localWindow=0.5*10^5, localSeperatePlot=TRUE, smoothF=0.025, xlabScale=10^6, width=12, height=18) {
	p = length(cases)/(length(cases)+length(controls))
	maxCase = max(cases)
	maxControl = max(controls)
	maxVal = max(c(maxCase, maxControl))
	cpts = matrix(CBSObj$statHat[,c(1,2,5)], nrow=nrow(CBSObj$statHat))
	tauHatFull = CBSObj$tauHat
	tauHat = tauHatFull[c(-1, -length(tauHatFull))]
	relCN = CBSObj$relCN
	ylims = c(min(c(0, cpts[,3])), max(c(0, cpts[,3])))
	relCNlims = c(min(relCN), max(relCN))
	if(!is.null(CIObj)) {
		CIBounds = CIObj$CIRes[3:4,]
		CIL = as.numeric(CIObj$CIRes[5,])
		CIU = as.numeric(CIObj$CIRes[6,])
		relCNCIL = CIL/(1-CIL+10^(-5))/(p/(1-p))
		relCNCIU = CIU/(1-CIU+10^(-5))/(p/(1-p))
	}
	
	grid.fix = seq(1, maxVal, length.out=length.out)
	gridSize = grid.fix[2]-grid.fix[1]
	casesCountInGrid = getCountsInWindow(cases, 0, maxVal, gridSize, sorted=FALSE)
	casesCountInGridSmooth = lowess(x=grid.fix, y=casesCountInGrid, smoothF)
	casesCountInGridSmooth$y[casesCountInGridSmooth$y<0] = 0
	controlCountInGrid = getCountsInWindow(controls, 0, maxVal, gridSize, sorted=FALSE)
	controlCountInGridSmooth = lowess(x=grid.fix, y=controlCountInGrid, smoothF)
	controlCountInGridSmooth$y[controlCountInGridSmooth$y<0] = 0
	PInGrid = casesCountInGrid/(casesCountInGrid+controlCountInGrid)
	PInGrid[is.nan(PInGrid)]=0
	PInGridSmooth = lowess(x=grid.fix, y=PInGrid, smoothF)
	tauHatInGrid = grid.fix[tauHat %/% gridSize]/xlabScale
	gridYLims = c(min(log(casesCountInGrid+1) - log(controlCountInGrid+1)), log(max(casesCountInGrid, controlCountInGrid)))
	
	## 1. Plot the chromosome global view
	pdf(paste(filename, ".pdf", sep=""), width=width, height=height)
	par(mfrow=c(3,1))
	
	plot(x=grid.fix/xlabScale, y=rep(0, length(grid.fix)), type="n", ylim=ylims, main=mainTitle, ylab="Statistic", xlab=paste("Base Pairs", xlabScale))
	for(i in 1:nrow(cpts)) {
		plotX = c(grid.fix[max(floor(cpts[i,1]/gridSize), 1)]/xlabScale, grid.fix[ceiling(cpts[i,2]/gridSize)]/xlabScale)
		lines(x=plotX, y=rep(cpts[i,3],2), lwd=3)
	}
	abline(v=tauHatInGrid, lty=3, col=4)
	
	matplot(x=grid.fix/xlabScale, y=log(cbind(casesCountInGridSmooth$y, controlCountInGridSmooth$y)+1), type="l", lty=c(1,1), col=c(2,1), main="Log Read Intensity", ylab="Read Intensity", xlab=paste("Base Pairs", xlabScale), ylim=gridYLims)
	points(x=grid.fix/xlabScale, y=log(casesCountInGrid+1) - log(controlCountInGrid+1), pch=".", col=1)
	abline(v=tauHatInGrid, lty=3, col=4)
	legend("topright", c("case","control", "case-control"), pch=".", lty=c(1,1,0), col=c(2,1,1))
	
	plotTauHat=cbind(c(1,tauHat), c(tauHat, maxVal))
	plot(x=grid.fix/xlabScale, y=rep(0, length(grid.fix)), type="n", ylim=relCNlims, main="Relative Copy Number", ylab="Relative CN", xlab=paste("Base Pairs", xlabScale))
	for(i in 1:nrow(plotTauHat)) {
		plotx = c(grid.fix[max(floor(plotTauHat[i,1]/gridSize), 1)]/xlabScale, grid.fix[ceiling(plotTauHat[i,2]/gridSize)]/xlabScale)
		lines(x=plotx, y=rep(relCN[i], 2), lwd=3)
	}
	abline(v=tauHatInGrid, lty=3, col=4)
	dev.off()
	
	## 2. Plot Local View of each Change Point
	nTauHat = length(tauHat)
	if(localSeperatePlot == FALSE) {
		nPlotCol = as.integer(sqrt(nTauHat/(height/width)))
		nPlotRow = ceiling(nTauHat/nPlotCol)
		pdf(paste(filename, "_localDetails.pdf", sep=""), width=width*2, height=height*2)
		par(mfrow=c(nPlotRow, nPlotCol))
	}
	for(i in 1:nTauHat) {
		if(localSeperatePlot) {
			pdf(paste(filename, "_local_", i, "_", tauHat[i], ".pdf", sep=""), width=width, height=height/2)
		}
		lBound = max(0, tauHat[i]-localWindow)
		rBound = min(maxVal, tauHat[i]+localWindow)
		localCas = cases[cases >= lBound & cases < rBound]
		localCon = controls[controls >= lBound & controls < rBound]
		grid.fix = seq(lBound, rBound, length.out=length.out/100)
		gridSize = grid.fix[2]-grid.fix[1]
		grid.mpt = grid.fix + gridSize/2
		CasCountInGrid = getCountsInWindow(localCas, lBound, rBound, gridSize, sorted=FALSE)
		ConCountInGrid = getCountsInWindow(localCon, lBound, rBound, gridSize, sorted=FALSE)
		pInGrid = CasCountInGrid/(CasCountInGrid+ConCountInGrid)
		pInGrid[is.nan(pInGrid)] = 0.0
		
		combLocalCasCon = CombineCaseControlC(localCas, localCon)
		plotReadX = combLocalCasCon$combL
		plotReadY = combLocalCasCon$combZ/combLocalCasCon$combX
		
		plotPX = cbind(tauHatFull[-length(tauHatFull)], tauHatFull[-1])
		pSegment = relCN*p/(1-p)/(1+relCN*p/(1-p))
		plotPY = cbind(pSegment, pSegment)
		
		plot(x=plotReadX, y=plotReadY, pch=".", xlim=c(lBound, rBound), ylim=c(0,1), main=paste("Reads and Inference around", tauHat[i]), xlab="Base Pair Locations", ylab="p(case read)", cex.main=0.75, cex.lab=0.5, cex.axis=0.5)
		if(length(grid.mpt) != length(pInGrid)) {
			grid.mpt = grid.mpt[1:max(length(grid.mpt), length(pInGrid))]
			pInGrid = pInGrid[1:max(length(grid.mpt), length(pInGrid))]
		}
		points(x=grid.mpt, y=pInGrid, pch=20, col=3)
		if(!is.null(CIObj)) {
			for(j in 1:ncol(CIBounds)) {
				lines(x=CIBounds[,j], y=rep(CIL[j],2), col=4, lwd=2)
				lines(x=CIBounds[,j], y=rep(CIU[j],2), col=2, lwd=2)
			}
		}
		for(j in 1:nrow(plotPY)) {
			lines(x=plotPX[j,], y=plotPY[j,], lwd=2.5)
		}
		abline(v=tauHat, lty=3, col=4)
		if(localSeperatePlot) {
			dev.off()
		}
	}
	if(localSeperatePlot == FALSE)	dev.off()
}
