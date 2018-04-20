#include "raDI.h"



void MI_Compute (int L, int q, double **f, double ****ff, double **MI){
	int i, j, A, B;
	double daux;
	
	for (i=0; i<L; i++){
		MI[i][i] = 0;
		for (j=i+1; j<L; j++){
			for (A=0; A<q; A++){
				for (B=0; B<q; B++){
					daux = ff[i][j][A][B]/(f[i][A]*f[j][B]);
					if (daux < 1e-10){
						continue;
					}
					MI[i][j] = MI[i][j] + ff[i][j][A][B]*log(daux);
				}
			}
			MI[j][i] = MI[i][j];
		}
	}	

	return;
}


