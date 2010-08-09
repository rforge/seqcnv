
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Print.h>
#include <Rmath.h>
#include <math.h>

/*
# CombineSortedVectorC.c
# Jeremy J. Shen
# Combine two vectors sorted in increasing order
# Updated: 4/13/2010
*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Print.h>
#include <math.h>

SEXP CombineSortedVectorC(SEXP casesS, SEXP controlsS) {
	double *cases = REAL(casesS);
	double *controls = REAL(controlsS);
	long long nCas = length(casesS);
	long long nCon = length(controlsS);
	long long nTot = nCas+nCon;
	SEXP combCC;
	PROTECT(combCC = allocVector(REALSXP, nTot));
	double *combCCPtr = REAL(combCC);
	long long i, fCas, fCon;
	fCas=fCon=0;
	for(i=0; i<nTot; i++) {
		if(fCon >= nCon) {
			combCCPtr[i] = cases[fCas];
			fCas++;
		}
		else if(fCas >= nCas) {
			combCCPtr[i] = controls[fCon];
			fCon++;
		}
		else if(cases[fCas] < controls[fCon]) {
			combCCPtr[i] = cases[fCas];
			fCas++;
		}
		else {
			combCCPtr[i] = controls[fCon];
			fCon++;
		}
	}
	UNPROTECT(1);
	return(combCC);
}

SEXP CombineToUniqueValueC(SEXP casesS, SEXP controlsS, SEXP combLS) {
	double *cases = REAL(casesS);
	double *controls = REAL(controlsS);
	double *combL = REAL(combLS);
	long long nCas = length(casesS);
	long long nCon = length(controlsS);
	long long nTot = nCas+nCon;
	long long nL = length(combLS);
	SEXP combZX;
	PROTECT(combZX = allocMatrix(REALSXP, nL, 2));
	double *combZXPtr = REAL(combZX);
	long long i, fCas, fCon, fL, nCasL, nConL;
	fCas=fCon=fL=0;
	for(i=0; i<nL; i++) {
		nCasL=nConL=0;
		while(fCas<nCas && cases[fCas]==combL[i]) {
			nCasL++;
			fCas++;
		}
		while(fCon<nCon && controls[fCon]==combL[i]) {
			nConL++;
			fCon++;
		}
		combZXPtr[i] = nCasL;
		combZXPtr[i+nL] = nConL+nCasL;
	}
	UNPROTECT(1);
	return(combZX);
}

SEXP FindUniqueInSortedArrayC(SEXP combCCS) {
	double *combCC = REAL(combCCS);
	long long i, j, nUnique;
	long long nEntry = length(combCCS);
	for(i=1, nUnique=1; i<nEntry; i++) {
		if(combCC[i]!=combCC[i-1])	nUnique++;
	}
	SEXP combL;
	PROTECT(combL = allocVector(REALSXP, nUnique));
	double *combLPtr = REAL(combL);
	combLPtr[0] = combCC[0];
	for(i=1, j=1; i<nEntry; i++) {
		if(combCC[i]!=combCC[i-1]) {
			combLPtr[j] = combCC[i];
			j++; 
		}
	}
	UNPROTECT(1);
	return(combL);
}


SEXP ScanIGSGridCumSumC(SEXP combYS, SEXP gridCurS) {
	double *combY = REAL(combYS);
	double *gridCur = REAL(gridCurS);
	long long gridCurLen = length(gridCurS);
	
	SEXP combYCumSum;
	PROTECT(combYCumSum = allocVector(REALSXP, gridCurLen));
	double *combYCumSumPtr = REAL(combYCumSum);
	long long i, j;
	
	combYCumSumPtr[0] = combY[0];
	for(i=1; i<gridCurLen; i++) {
		combYCumSumPtr[i] = combYCumSumPtr[i-1];
		for(j=gridCur[i-1]; j<gridCur[i]; j++) {
			combYCumSumPtr[i] += combY[j];
		}
	}
	UNPROTECT(1);
	return(combYCumSum);
}



SEXP ScanStatNewCompBinomC(SEXP combZCumSumS, SEXP combXCumSumS, SEXP combZPointS, SEXP combXPointS, SEXP pS, SEXP nTotalS, SEXP gridCurS, SEXP maxWinS) {
	/*
	combZCumSumS = coerceVector(combZCumSumS, INTSXP);
	combXCumSumS = coerceVector(combXCumSumS, INTSXP);
	combZPointS = coerceVector(combZPointS, INTSXP);
	combXPointS = coerceVector(combXPointS, INTSXP);
	nTotalS = coerceVector(nTotalS, INTSXP);
	gridCurS = coerceVector(gridCurS, INTSXP);
	maxWinS = coerceVector(maxWinS, INTSXP);
	*/

	double *combZCumSum = REAL(combZCumSumS);
	double *combXCumSum = REAL(combXCumSumS);
	double *combZPoint = REAL(combZPointS);
	double *combXPoint = REAL(combXPointS);
	double p = REAL(pS)[0];
	double nTotal = REAL(nTotalS)[0];
	long long gridCurLen = length(gridCurS);
	long long gridCurMaxInd = gridCurLen-1;
	double *gridCur = REAL(gridCurS);
	double maxWin = REAL(maxWinS)[0];
	long long i, j, jMax, bestWinI, bestWinJ;
	double nCas, nObs, pNObs, Rij, bestWinR, bestWinRAbs;
	SEXP newRes;
	PROTECT(newRes = allocMatrix(REALSXP, gridCurMaxInd, 3));
	double *newResPtr = REAL(newRes);
	
	for(i=0.0; i<gridCurMaxInd; i++) {
		jMax = i + maxWin;
		if(jMax > gridCurMaxInd) {
			jMax=gridCurMaxInd;
		}
		bestWinI = i;
		bestWinJ = jMax;
		bestWinR = 0.0;
		bestWinRAbs = 0.0;
		for(j=i+1; j<=jMax; j++) {
			nObs = combXCumSum[j]-combXCumSum[i]+combXPoint[i];
			if (nObs == 0.0) {
				Rij = 0.0;
			}
			else {
				nCas = combZCumSum[j]-combZCumSum[i]+combZPoint[i];
				pNObs = p*nObs;
				if(nCas < pNObs){
					Rij = pbinom(nCas, nObs, p, 1,1);
					Rij = qnorm(Rij, 0,1,1,1);
					if(R_FINITE(Rij) == 0) Rij=0;
				}
				else {
					Rij = pbinom(nCas, nObs, p, 0,1);
					Rij = qnorm(Rij, 0,1,0,1);
					if(R_FINITE(Rij) == 0) Rij=0;
				}
			}
			if(fabs(Rij) > bestWinRAbs) {
				bestWinI = i;
				bestWinJ = j;
				bestWinR = Rij;
				bestWinRAbs = fabs(bestWinR);
			}
		}
		newResPtr[i] = gridCur[bestWinI];
		newResPtr[i + gridCurMaxInd] = gridCur[bestWinJ];
		newResPtr[i + 2*gridCurMaxInd] = bestWinR;
	}
	UNPROTECT(1);
	return(newRes);
}


SEXP ScanStatRefineCompBinomC(SEXP combZCumSumS, SEXP combXCumSumS, SEXP combZPointS, SEXP combXPointS, SEXP pS, SEXP nTotalS, SEXP gridCurS, SEXP gridLRS, SEXP maxWinS) {
	/*
	combZCumSumS = coerceVector(combZCumSumS, INTSXP);
	combXCumSumS = coerceVector(combXCumSumS, INTSXP);
	combZPointS = coerceVector(combZPointS, INTSXP);
	combXPointS = coerceVector(combXPointS, INTSXP);
	nTotalS = coerceVector(nTotalS, INTSXP);
	gridCurS = coerceVector(gridCurS, INTSXP);
	gridLS = coerceVector(gridLS, INTSXP);
	gridRS = coerceVector(gridRS, INTSXP);
	maxWinS = coerceVector(maxWinS, INTSXP);
	*/

	double *combZCumSum = REAL(combZCumSumS);
	double *combXCumSum = REAL(combXCumSumS);
	double *combZPoint = REAL(combZPointS);
	double *combXPoint = REAL(combXPointS);
	double p = REAL(pS)[0];
	double nTotal = REAL(nTotalS)[0];
	long long gridCurLen = length(gridCurS);
	long long gridCurMaxInd = gridCurLen-1;
	double *gridCur = REAL(gridCurS);
	double gridL = REAL(gridLRS)[0]-1;
	double gridR = REAL(gridLRS)[1]-1;
	double maxWin = REAL(maxWinS)[0];
	double jMin, gridLL, gridLR, gridRL, gridRR;
	long long i, j, nRows, bestWinI, bestWinJ, rCt;
	double nCas, nObs, pNObs, Rij, bestWinR, bestWinRAbs;

	gridLL = gridL - floor(maxWin/2);
	if(gridLL < 0) {
		gridLL = 0;
	}
	gridLR = gridL + floor(maxWin/2);
	if(gridLR > gridCurMaxInd-1) {
		gridLR = gridCurMaxInd-1;
	}
	gridRL = gridR - floor(maxWin/2);
	if(gridRL < 0) {
		gridRL = 0;
	}
	gridRR = gridR + floor(maxWin/2);
	if(gridRR > gridCurMaxInd) {
		gridRR = gridCurMaxInd;
	}
	nRows = gridLR-gridLL+1;
	
	SEXP newRes;
	PROTECT(newRes = allocMatrix(REALSXP, nRows, 3));
	double *newResPtr = REAL(newRes);
	rCt=0;
	
	//printf("gridL: %d \t gridR: %d \t maxWin: %d\n", gridL, gridR, maxWin);
	//printf("gridLL: %d \t gridLR: %d \t gridRL: %d \t gridRR: %d \t nRows: %d\n", gridLL, gridLR, gridRL, gridRR, nRows);
	
	for(i=gridLL; i<=gridLR; i++) {
		jMin = i+1;
		if(jMin > gridRL) {
			jMin=gridRL;
		}
		bestWinI = i;
		bestWinJ = gridRR;
		bestWinR = 0.0;
		bestWinRAbs = 0.0;
		for(j=jMin; j<=gridRR; j++) {
			nObs = combXCumSum[j]-combXCumSum[i]+combXPoint[i];
			if (nObs == 0.0) {
				Rij = 0.0;
			}
			else {
				nCas = combZCumSum[j]-combZCumSum[i]+combZPoint[i];
				pNObs = p*nObs;
				if(nCas < pNObs){
					Rij = pbinom(nCas, nObs, p, 1,1);
					Rij = qnorm(Rij, 0,1,1,1);
					if(R_FINITE(Rij) == 0) Rij=0;
				}
				else {
					Rij = pbinom(nCas, nObs, p, 0,1);
					Rij = qnorm(Rij, 0,1,0,1);
					if(R_FINITE(Rij) == 0) Rij=0;
				}
			}
			if(fabs(Rij) > bestWinRAbs) {
				bestWinI = i;
				bestWinJ = j;
				bestWinR = Rij;
				bestWinRAbs = fabs(bestWinR);
			}
		}
		newResPtr[rCt] = gridCur[bestWinI];
		newResPtr[rCt + nRows] = gridCur[bestWinJ];
		newResPtr[rCt + 2*nRows] = bestWinR;
		rCt++;
	}
	UNPROTECT(1);
	return(newRes);
}

SEXP ScanStatNewCompNormalC(SEXP combZCumSumS, SEXP combXCumSumS, SEXP combZPointS, SEXP combXPointS, SEXP SijFactorNS, SEXP pS, SEXP nTotalS, SEXP gridCurS, SEXP maxWinS) {
	/*
	combZCumSumS = coerceVector(combZCumSumS, INTSXP);
	combXCumSumS = coerceVector(combXCumSumS, INTSXP);
	combZPointS = coerceVector(combZPointS, INTSXP);
	combXPointS = coerceVector(combXPointS, INTSXP);
	nTotalS = coerceVector(nTotalS, INTSXP);
	gridCurS = coerceVector(gridCurS, INTSXP);
	maxWinS = coerceVector(maxWinS, INTSXP);
	*/

	double *combZCumSum = REAL(combZCumSumS);
	double *combXCumSum = REAL(combXCumSumS);
	double *combZPoint = REAL(combZPointS);
	double *combXPoint = REAL(combXPointS);
	double SijFactorN = REAL(SijFactorNS)[0];
	double p = REAL(pS)[0];
	double nTotal = REAL(nTotalS)[0];
	long long gridCurLen = length(gridCurS);
	long long gridCurMaxInd = gridCurLen-1;
	double *gridCur = REAL(gridCurS);
	double maxWin = REAL(maxWinS)[0];
	long long i, j, jMax, bestWinI, bestWinJ;
	double nCas, nObs, pNObs, SijFactor2, Rij, bestWinR, bestWinRAbs;
	SEXP newRes;
	PROTECT(newRes = allocMatrix(REALSXP, gridCurMaxInd, 3));
	double *newResPtr = REAL(newRes);
	
	for(i=0.0; i<gridCurMaxInd; i++) {
		jMax = i + maxWin;
		if(jMax > gridCurMaxInd) {
			jMax=gridCurMaxInd;
		}
		bestWinI = i;
		bestWinJ = jMax;
		bestWinR = 0.0;
		bestWinRAbs = 0.0;
		for(j=i+1; j<=jMax; j++) {
			nObs = combXCumSum[j]-combXCumSum[i]+combXPoint[i];
			SijFactor2 = nObs;
			if (SijFactor2 == 0.0) {
				Rij = 0.0;
			}
			else {
				nCas = combZCumSum[j]-combZCumSum[i]+combZPoint[i];
				pNObs = p*nObs;
				Rij = (nCas - pNObs)/(sqrt(SijFactorN*SijFactor2));
			}
			if(fabs(Rij) > bestWinRAbs) {
				bestWinI = i;
				bestWinJ = j;
				bestWinR = Rij;
				bestWinRAbs = fabs(bestWinR);
			}
		}
		newResPtr[i] = gridCur[bestWinI];
		newResPtr[i + gridCurMaxInd] = gridCur[bestWinJ];
		newResPtr[i + 2*gridCurMaxInd] = bestWinR;
	}
	UNPROTECT(1);
	return(newRes);
}


SEXP ScanStatRefineCompNormalC(SEXP combZCumSumS, SEXP combXCumSumS, SEXP combZPointS, SEXP combXPointS, SEXP SijFactorNS, SEXP pS, SEXP nTotalS, SEXP gridCurS, SEXP gridLRS, SEXP maxWinS) {
	/*
	combZCumSumS = coerceVector(combZCumSumS, INTSXP);
	combXCumSumS = coerceVector(combXCumSumS, INTSXP);
	combZPointS = coerceVector(combZPointS, INTSXP);
	combXPointS = coerceVector(combXPointS, INTSXP);
	nTotalS = coerceVector(nTotalS, INTSXP);
	gridCurS = coerceVector(gridCurS, INTSXP);
	gridLS = coerceVector(gridLS, INTSXP);
	gridRS = coerceVector(gridRS, INTSXP);
	maxWinS = coerceVector(maxWinS, INTSXP);
	*/

	double *combZCumSum = REAL(combZCumSumS);
	double *combXCumSum = REAL(combXCumSumS);
	double *combZPoint = REAL(combZPointS);
	double *combXPoint = REAL(combXPointS);
	double SijFactorN = REAL(SijFactorNS)[0];
	double p = REAL(pS)[0];
	double nTotal = REAL(nTotalS)[0];
	long long gridCurLen = length(gridCurS);
	long long gridCurMaxInd = gridCurLen-1;
	double *gridCur = REAL(gridCurS);
	double gridL = REAL(gridLRS)[0]-1;
	double gridR = REAL(gridLRS)[1]-1;
	double maxWin = REAL(maxWinS)[0];
	double jMin, gridLL, gridLR, gridRL, gridRR;
	long long i, j, nRows, bestWinI, bestWinJ, rCt;
	double nCas, nObs, pNObs, SijFactor2, Rij, bestWinR, bestWinRAbs;

	gridLL = gridL - floor(maxWin/2);
	if(gridLL < 0) {
		gridLL = 0;
	}
	gridLR = gridL + floor(maxWin/2);
	if(gridLR > gridCurMaxInd-1) {
		gridLR = gridCurMaxInd-1;
	}
	gridRL = gridR - floor(maxWin/2);
	if(gridRL < 0) {
		gridRL = 0;
	}
	gridRR = gridR + floor(maxWin/2);
	if(gridRR > gridCurMaxInd) {
		gridRR = gridCurMaxInd;
	}
	nRows = gridLR-gridLL+1;
	
	SEXP newRes;
	PROTECT(newRes = allocMatrix(REALSXP, nRows, 3));
	double *newResPtr = REAL(newRes);
	rCt=0;
	
	//printf("gridL: %d \t gridR: %d \t maxWin: %d\n", gridL, gridR, maxWin);
	//printf("gridLL: %d \t gridLR: %d \t gridRL: %d \t gridRR: %d \t nRows: %d\n", gridLL, gridLR, gridRL, gridRR, nRows);
	
	for(i=gridLL; i<=gridLR; i++) {
		jMin = i+1;
		if(jMin > gridRL) {
			jMin=gridRL;
		}
		bestWinI = i;
		bestWinJ = gridRR;
		bestWinR = 0.0;
		bestWinRAbs = 0.0;
		for(j=jMin; j<=gridRR; j++) {
			nObs = combXCumSum[j]-combXCumSum[i]+combXPoint[i];
			SijFactor2 = nObs;
			if (SijFactor2 == 0.0) {
				Rij = 0.0;
			}
			else {
				nCas = combZCumSum[j]-combZCumSum[i]+combZPoint[i];
				pNObs = p*nObs;
				Rij = (nCas - pNObs)/(sqrt(SijFactorN*SijFactor2));
			}
			if(fabs(Rij) > bestWinRAbs) {
				bestWinI = i;
				bestWinJ = j;
				bestWinR = Rij;
				bestWinRAbs = fabs(bestWinR);
			}
		}
		newResPtr[rCt] = gridCur[bestWinI];
		newResPtr[rCt + nRows] = gridCur[bestWinJ];
		newResPtr[rCt + 2*nRows] = bestWinR;
		rCt++;
	}
	UNPROTECT(1);
	return(newRes);
}


SEXP ScanStatNewCompRabinC(SEXP combZCumSumS, SEXP combXCumSumS, SEXP combZPointS, SEXP combXPointS, SEXP SijFactorRS, SEXP pS, SEXP nTotalS, SEXP gridCurS, SEXP maxWinS) {
	/*
	combZCumSumS = coerceVector(combZCumSumS, INTSXP);
	combXCumSumS = coerceVector(combXCumSumS, INTSXP);
	combZPointS = coerceVector(combZPointS, INTSXP);
	combXPointS = coerceVector(combXPointS, INTSXP);
	nTotalS = coerceVector(nTotalS, INTSXP);
	gridCurS = coerceVector(gridCurS, INTSXP);
	maxWinS = coerceVector(maxWinS, INTSXP);
	*/

	double *combZCumSum = REAL(combZCumSumS);
	double *combXCumSum = REAL(combXCumSumS);
	double *combZPoint = REAL(combZPointS);
	double *combXPoint = REAL(combXPointS);
	double SijFactorR = REAL(SijFactorRS)[0];
	double p = REAL(pS)[0];
	double nTotal = REAL(nTotalS)[0];
	long long gridCurLen = length(gridCurS);
	long long gridCurMaxInd = gridCurLen-1;
	double *gridCur = REAL(gridCurS);
	double maxWin = REAL(maxWinS)[0];
	long long i, j, jMax, bestWinI, bestWinJ;
	double nCas, nObs, pNObs, SijFactor2, Rij, bestWinR, bestWinRAbs;
	SEXP newRes;
	PROTECT(newRes = allocMatrix(REALSXP, gridCurMaxInd, 3));
	double *newResPtr = REAL(newRes);
	
	for(i=0.0; i<gridCurMaxInd; i++) {
		jMax = i + maxWin;
		if(jMax > gridCurMaxInd) {
			jMax=gridCurMaxInd;
		}
		bestWinI = i;
		bestWinJ = jMax;
		bestWinR = 0.0;
		bestWinRAbs = 0.0;
		for(j=i+1; j<=jMax; j++) {
			nObs = combXCumSum[j]-combXCumSum[i]+combXPoint[i];
			SijFactor2 = nObs - nObs*nObs/nTotal;
			if (SijFactor2 == 0.0) {
				Rij = 0.0;
			}
			else {
				nCas = combZCumSum[j]-combZCumSum[i]+combZPoint[i];
				pNObs = p*nObs;
				Rij = (nCas - pNObs)/(sqrt(SijFactorR*SijFactor2));
			}
			if(fabs(Rij) > bestWinRAbs) {
				bestWinI = i;
				bestWinJ = j;
				bestWinR = Rij;
				bestWinRAbs = fabs(bestWinR);
			}
		}
		newResPtr[i] = gridCur[bestWinI];
		newResPtr[i + gridCurMaxInd] = gridCur[bestWinJ];
		newResPtr[i + 2*gridCurMaxInd] = bestWinR;
	}
	UNPROTECT(1);
	return(newRes);
}


SEXP ScanStatRefineCompRabinC(SEXP combZCumSumS, SEXP combXCumSumS, SEXP combZPointS, SEXP combXPointS, SEXP SijFactorRS, SEXP pS, SEXP nTotalS, SEXP gridCurS, SEXP gridLRS, SEXP maxWinS) {
	/*
	combZCumSumS = coerceVector(combZCumSumS, INTSXP);
	combXCumSumS = coerceVector(combXCumSumS, INTSXP);
	combZPointS = coerceVector(combZPointS, INTSXP);
	combXPointS = coerceVector(combXPointS, INTSXP);
	nTotalS = coerceVector(nTotalS, INTSXP);
	gridCurS = coerceVector(gridCurS, INTSXP);
	gridLS = coerceVector(gridLS, INTSXP);
	gridRS = coerceVector(gridRS, INTSXP);
	maxWinS = coerceVector(maxWinS, INTSXP);
	*/

	double *combZCumSum = REAL(combZCumSumS);
	double *combXCumSum = REAL(combXCumSumS);
	double *combZPoint = REAL(combZPointS);
	double *combXPoint = REAL(combXPointS);
	double SijFactorR = REAL(SijFactorRS)[0];
	double p = REAL(pS)[0];
	double nTotal = REAL(nTotalS)[0];
	long long gridCurLen = length(gridCurS);
	long long gridCurMaxInd = gridCurLen-1;
	double *gridCur = REAL(gridCurS);
	double gridL = REAL(gridLRS)[0]-1;
	double gridR = REAL(gridLRS)[1]-1;
	double maxWin = REAL(maxWinS)[0];
	double jMin, gridLL, gridLR, gridRL, gridRR;
	long long i, j, nRows, bestWinI, bestWinJ, rCt;
	double nCas, nObs, pNObs, SijFactor2, Rij, bestWinR, bestWinRAbs;

	gridLL = gridL - floor(maxWin/2);
	if(gridLL < 0) {
		gridLL = 0;
	}
	gridLR = gridL + floor(maxWin/2);
	if(gridLR > gridCurMaxInd-1) {
		gridLR = gridCurMaxInd-1;
	}
	gridRL = gridR - floor(maxWin/2);
	if(gridRL < 0) {
		gridRL = 0;
	}
	gridRR = gridR + floor(maxWin/2);
	if(gridRR > gridCurMaxInd) {
		gridRR = gridCurMaxInd;
	}
	nRows = gridLR-gridLL+1;
	
	SEXP newRes;
	PROTECT(newRes = allocMatrix(REALSXP, nRows, 3));
	double *newResPtr = REAL(newRes);
	rCt=0;
	
	//printf("gridL: %d \t gridR: %d \t maxWin: %d\n", gridL, gridR, maxWin);
	//printf("gridLL: %d \t gridLR: %d \t gridRL: %d \t gridRR: %d \t nRows: %d\n", gridLL, gridLR, gridRL, gridRR, nRows);
	
	for(i=gridLL; i<=gridLR; i++) {
		jMin = i+1;
		if(jMin > gridRL) {
			jMin=gridRL;
		}
		bestWinI = i;
		bestWinJ = gridRR;
		bestWinR = 0.0;
		bestWinRAbs = 0.0;
		for(j=jMin; j<=gridRR; j++) {
			nObs = combXCumSum[j]-combXCumSum[i]+combXPoint[i];
			SijFactor2 = nObs - nObs*nObs/nTotal;
			if (SijFactor2 == 0.0) {
				Rij = 0.0;
			}
			else {
				nCas = combZCumSum[j]-combZCumSum[i]+combZPoint[i];
				pNObs = p*nObs;
				Rij = (nCas - pNObs)/(sqrt(SijFactorR*SijFactor2));
			}
			if(fabs(Rij) > bestWinRAbs) {
				bestWinI = i;
				bestWinJ = j;
				bestWinR = Rij;
				bestWinRAbs = fabs(bestWinR);
			}
		}
		newResPtr[rCt] = gridCur[bestWinI];
		newResPtr[rCt + nRows] = gridCur[bestWinJ];
		newResPtr[rCt + 2*nRows] = bestWinR;
		rCt++;
	}
	UNPROTECT(1);
	return(newRes);
}

double pBetaMixRootEval(double x, double pRoot, double *betaParam1, double *betaParam2, double *wks, long nMix) {
	double px = 0.0;
	for(long i=0; i<nMix; i++) {
		px = px + wks[i]*pbeta(x, betaParam1[i], betaParam2[i], 1, 0);
	}
	px = px - pRoot;
	return(px);
}

double dBetaMixEval(double x, double *betaParam1, double *betaParam2, double *wks, long nMix) {
	double dx = 0.0;
	for(long i=0; i<nMix; i++) {
		dx = dx + wks[i]*dbeta(x, betaParam1[i], betaParam2[i], 0);
	}
	return(dx);
}

double rtBetaMixCDF(double pRoot, double *betaParam1, double *betaParam2, double *wks, long nMix, double epsCDF) {
	// Newton-Raphson with Bisection Safety Mechanism
	// Concepts of Numerical Recipes 3rd edition used
	const int maxIter=1000;
	double x1 = 0;
	double x2 = 1;
	double xl = x1;
	double xh = x2;
	double fl = -pRoot;
	double fh = 1-pRoot;
	double rts, dxold, dx, f, df, temp;

	if (fl == 0.0) { 
		return x1;
	}
	if (fh == 0.0) {
		return x2;
	}

	rts = 0.5*(x1+x2); 
	dxold = fabs(x2-x1); 
	dx = dxold; 
	f = pBetaMixRootEval(rts, pRoot, betaParam1, betaParam2, wks, nMix);
	df = dBetaMixEval(rts, betaParam1, betaParam2, wks, nMix);
	for (int j=0; j<maxIter; j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) return rts;
		}
		else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < epsCDF) {
			return(rts);
		}
		f = pBetaMixRootEval(rts, pRoot, betaParam1, betaParam2, wks, nMix);
		df= dBetaMixEval(rts, betaParam1, betaParam2, wks, nMix);
		if (f < 0.0) {
			xl=rts;
		}
		else {
			xh=rts;
		}
	}	
}

SEXP BayesCptCICompC(SEXP betaParam1S, SEXP betaParam2S, SEXP wksS, SEXP alphaS, SEXP epsCDFS) {
	double *betaParam1 = REAL(betaParam1S);
	double *betaParam2 = REAL(betaParam2S);
	double *wks = REAL(wksS);
	double alpha = REAL(alphaS)[0];
	double epsCDF = REAL(epsCDFS)[0];
	long i, j;
	long nMix = length(betaParam1S);
	double pLow = alpha/2;
	double pHigh = 1.0 - pLow;
	
	SEXP CI;
	PROTECT(CI = allocVector(REALSXP, 2));
	double *CIPtr = REAL(CI);
	
	// Find Roots and Return
	CIPtr[0] = rtBetaMixCDF(pLow, betaParam1, betaParam2, wks, nMix, epsCDF);
	CIPtr[1] = rtBetaMixCDF(pHigh, betaParam1, betaParam2, wks, nMix, epsCDF);
	UNPROTECT(1);
	return(CI);
}
