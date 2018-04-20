#include "raDI.h"



/* Routine to compute Moore-Penrose pseudoinverse of A using SVD procedure, and store it into A */
void Pseudoinverse (int n, double **A){
	int i, j, k;
	double E, aux;
	double *D, **V, **AUX;

	E = 1e-3;

	D = (double*) calloc(n,sizeof(double)); /* Vector for the diagonal matrix of SVD descomposition */
	AUX = (double**) malloc(n*sizeof(double*));
	V = (double**) malloc(n*sizeof(double*));
	if (D == NULL || AUX == NULL || V == NULL){
		printf("It was not possible to allocate memory\n");
		exit(1);
	}
	for (i=0; i<n; i++){
		AUX[i] = (double*) calloc(n,sizeof(double));
		V[i] = (double*) calloc(n,sizeof(double));
		if (V[i] == NULL || AUX == NULL){
			printf("It was not possible to allocate memory\n");
			exit(1);
		}
	}
	/* Let's perform a SVD decomposition A=UDV^t, where U will be stored into A */
	/* warning: V is not transposed in the result! */
	svdcmp(A,n,n,D,V);
	/* Compute the pseudoinverse of the diagonal matrix D */
	for (i=0; i<n; i++){
		if (fabs(D[i])>E){
			D[i] = 1./D[i];
		}
		else{
			D[i] = 0.;
		}
	}

	/* Undo the basis change for SVD descomposition, with the pseudoinverse of D, in order to obtain the pseudoinverse of A */
	/* Remember that we have to transpose V in order to undo the basis change! And that D is a diagonal matrix! */
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			A[i][j] = A[i][j]*D[j];
		}
	}
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			AUX[i][j] = 0.;
			for (k=0; k<n; k++){
				AUX[i][j] = AUX[i][j] + A[i][k]*V[j][k];
			}
		}
	}
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			A[i][j] = AUX[i][j];
		}
	}

	/* Free memory for V and D */
	for (i=0; i<n; i++){
		free(AUX[i]);
		free(V[i]);
	}
	free(AUX);
	free(V);
	free(D);

	return;
}
