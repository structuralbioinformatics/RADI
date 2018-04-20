#include "raDI.h"



void Correction_APC (int L, double **MIo, double **MI){
	int i, j ;
	double apc,aaMI,*aMI;
	
        aMI = (double*) calloc(L,sizeof(double));

        aaMI = 0;
	for (i=0; i<L; i++){
	  aMI[i] = 0;
	  for (j=0; j<L; j++){
            MI[i][j]          = 0;
            aaMI              = aaMI   + MIo[i][j];
	    if (j!=i){ aMI[i] = aMI[i] + MIo[i][j]; }
          }
          aMI[i] = aMI[i] / L ;
        }
        aaMI = aaMI / (L*L) ;
        for (i=0; i<L; i++){
        for (j=i+1; j<L; j++){
            apc      = aMI[i]*aMI[j]/aaMI;
            MI[i][j] = MIo[i][j] - apc;
            MI[j][i] = MIo[j][i] - apc;
        }}

	return;
}


