#include "raDI.h"


void MaxEnt (int L, int q, double **f, double ****ff, double **DI){

	int i, j, k, l, q1, Lq, A, B, A1, B1, count, Maxcount, auxexit;
	double c, daux, mu0, sum, sum1, sum2, dist, eps, Z;
	double *mu1, *mu2, *scra1, *scra2, *new1, *new2;
	double **C, **h, **h_new, **Pdir;
	double ****e;

	Maxcount = 1e8;
	eps = 1e-8;

	q1 = q-1;
	Lq = L*q1;

	C = (double**) malloc(Lq*sizeof(double*));
	e = (double****) malloc(L*sizeof(double***));
	if (C == NULL || e == NULL){
		printf("It was not possible to allocate memory \n");
		exit(1);
	}
	for (i=0; i<L; i++){
		e[i] = (double***) malloc(L*sizeof(double**));
		if (e[i] == NULL){
			printf("It was not possible to allocate memory (e)\n");
			exit(1);
		}
		for (j=0; j<L; j++){
			e[i][j] = (double**) malloc(q*sizeof(double*));
			if (e[i][j] == NULL){
				printf("It was not possible to allocate memory (e)\n");
				exit(1);
			}
			for (A=0; A<q; A++){
				e[i][j][A] = (double*) calloc(q,sizeof(double));
				if (e[i][j][A] == NULL){
					printf("It was not possible to allocate memory (e)\n");
					exit(1);
				}
			}
		}
	}
	for(i=0; i<Lq; i++){
		C[i] = (double*) calloc(Lq,sizeof(double));
		if (C[i] == NULL){
			printf("It was not possible to allocate memory (C)\n");
			exit(1);
		}
	}

	/* Gauge invariance method for obtaining interaction lagrange multipliers e */
	for (i=0; i<L; i++){
		for (j=0; j<L; j++){
			for (A=1; A<q; A++){
				for (B=1; B<q; B++){
					A1 = A-1;
					B1 = B-1;
					C[i*q1+A1][j*q1+B1] = ff[i][j][A][B]-f[i][A]*f[j][B];
				}
			}
		}
	}

	Pseudoinverse(Lq, C);

	for (i=0; i<Lq; i++){
		for (j=0; j<Lq; j++){
			k = i/q1;
			A = i%q1+1;
			l = j/q1;
			B = j%q1+1;
			e[k][l][A][B] = exp(-C[i][j]);
		}
	}
	for (i=0; i<L; i++){
		for (j=0; j<L; j++){
			e[i][j][0][0] = 1.;
			for (A=0; A<q; A++){
				e[i][j][A][0] = 1.;
				e[i][j][0][A] = 1.;
			}
		}
	}

	for (i=0; i<Lq; i++){
		free(C[i]);
	}
	free(C);

	h = (double**) malloc(L*sizeof(double*));
	h_new = (double**) malloc(L*sizeof(double*));
	if (h == NULL || h_new ==NULL){
		printf("It was not possible to allocate memory\n");
		exit(1);
	}
	for (i=0; i<L; i++){
		h[i] = (double*) malloc(q*sizeof(double));
		h_new[i] = (double*) malloc(q*sizeof(double));
		if (h[i] == NULL || h_new[i] ==  NULL){
			printf("It was not possible to allocate memory\n");
			exit(1);
		}
	}
	Pdir = (double**) malloc(q*sizeof(double*));
	if (Pdir == NULL){
		printf("It was not possible to allocate memory\n");
		exit(1);
	}
	for (A=0; A<q; A++){
		Pdir[A] = (double*) calloc(q,sizeof(double));
		if (Pdir == NULL){
			printf("It was not possible to allocate memory\n");
			exit(1);
		}
	}
	mu1 = (double*) malloc(q*sizeof(double));
	mu2 = (double*) malloc(q*sizeof(double));
	scra1 = (double*) malloc(q*sizeof(double));
	scra2 = (double*) malloc(q*sizeof(double));
	new1 = (double*) malloc(q*sizeof(double));
	new2 = (double*) malloc(q*sizeof(double));
	mu0 = exp(1./q);
	for (i=0; i<L-1; i++){
		DI[i][i] = 1.;
		for (j=i+1; j<L; j++){
			for (A=0; A<q; A++){
				mu1[A] = mu0;
				mu2[A] = mu0;
			}
			for (count=0; count<Maxcount; count++){
				dist = 0;
				for (A=0; A<q; A++){
					scra1[A] = 0;
					scra2[A] = 0;
					for (B=0; B<q; B++){
						scra1[A] = scra1[A] + e[i][j][A][B]*mu2[B];
						scra2[A] = scra2[A] + e[i][j][B][A]*mu1[B];
					}
				}
				sum1 = 0;
				sum2 = 0;
				for (A=0; A<q; A++){
					new1[A] = f[i][A]/scra1[A];
					new2[A] = f[j][A]/scra2[A];
					sum1 = sum1 + new1[A];
					sum2 = sum2 + new2[A];
				}
				for (A=0; A<q; A++){
					new1[A] = new1[A]/sum1;
					daux = fabs(new1[A]-mu1[A]);
					if (daux > dist){
						dist = daux;
					}
					new2[A] = new2[A]/sum2;
					daux = fabs(new2[A]-mu2[A]);
					if (daux > dist){
						dist = daux;
					}
					mu1[A] = new1[A];
					mu2[A] = new2[A];
				}
				if (dist < eps){
					break;
				}
			}
			if (count == Maxcount){
				printf("Lagrange single-site marginals could not be obtained in %d iterations\n", Maxcount);
				exit(9);
			}
			Z = 0.;		
			for (A=0; A<q; A++){
				for (B=0; B<q; B++){
					Z = Z + e[i][j][A][B]*mu1[A]*mu2[B];
				}
			}
			for (A=0; A<q; A++){
				for (B=0; B<q; B++){
					Pdir[A][B] = (e[i][j][A][B]*mu1[A]*mu2[B])/Z;
				}
			}
			for(A=0; A<q; A++){
				for (B=0; B<q; B++){
					daux = Pdir[A][B]/(f[i][A]*f[j][B]);
					if (daux < 1e-100){
						continue;
					}
					DI[i][j] = DI[i][j] + Pdir[A][B]*log(daux);
				}
			}
			DI[j][i] = DI[i][j];
		}
	}
	free(mu1);
	free(mu2);
	free(scra1);
	free(scra2);
	free(new1);
	free(new2);

	for (i=0; i<L; i++){
		for (j=0; j<L; j++){
			for (A=0; A<q; A++){
				free(e[i][j][A]);
			}
			free(e[i][j]);
		}
		free(e[i]);
		free(h[i]);
		free(h_new[i]);
	}
	free(e);
	free(h);
	free(h_new);

	/* Free memory for Pdir */
	for (A=0; A<q; A++){
		free(Pdir[A]);
	}
	free(Pdir);

	return;	
}


