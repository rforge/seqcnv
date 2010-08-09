relCNComp <-
function(combX, combZ, tauHatInd, p) {
	baselineRelCN = p/(1-p)
	nTau = length(tauHatInd)-1
	relCN = rep(0, nTau)
	for(i in 1:nTau) {
		tolRead = sum(combX[tauHatInd[i]:tauHatInd[i+1]])
		caseRead = sum(combZ[tauHatInd[i]:tauHatInd[i+1]])
		relCN[i] = (caseRead/(tolRead-caseRead+1))/baselineRelCN
	}
	return(relCN)
}

