#include "raDI.h"

double EffectiveCounts( MSA,mini_aa,aa,q,numseqs, L, l,pos,f,ff,verbose)
char **MSA;
char *mini_aa, *aa;
int  numseqs, L, l,q,*pos,verbose;
double **f,****ff;
{
  int i,j,k,aux,A,B,id;
  double Max_id_use, Max_lid, Max_gaps_use, Meff, lambda_psc_use, lambda_Meff;
  double *seqs;
//Request memory
	seqs = (double*) calloc(numseqs,sizeof(double));

//Initialize parameters  
        Max_id_use     = (double) Max_id;
        Max_gaps_use   = (double) Max_gaps;
        lambda_psc_use = (double) lambda_psc;

//Obtain effective number of sequences
	Max_lid = Max_id_use*l;
	Meff = 0;
        if (verbose) {printf("\t--Parameters Max id %f  \n",Max_lid);}
	for (i=0; i<numseqs; i++){
		seqs[i] = 0;
		for (j=0; j<numseqs; j++){
			id = 0;
			for (k=0; k<l; k++){
				if (MSA[i][pos[k]] == MSA[j][pos[k]]){ id++; }
			}
			if (id >= Max_lid){ seqs[i]++; }
		}
		seqs[i] = 1./seqs[i];
		Meff = Meff + seqs[i];
	}
	lambda_psc_use = Meff; //Best value for lambda_psc
	lambda_Meff    = lambda_psc_use + Meff;
        lambda_psc_use = lambda_psc_use/q;

//Single AA frequencies
	for (i=0; i<l; i++){
		for (A=0; A<q; A++){
			f[i][A] = lambda_psc_use;
			ff[i][i][A][A] = lambda_psc_use; 
			for (k=0; k<numseqs; k++){
				if (MSA[k][pos[i]] == aa[A]){ f[i][A] = f[i][A] + seqs[k]; }
			}
			f[i][A] = f[i][A]/lambda_Meff;
		}
	}
//Pair AA frequencies
        lambda_psc_use = lambda_psc_use/q;
	for (i=0; i<l; i++){
	for (j=0; j<l; j++){
		for (A=0; A<q; A++){
		for (B=0; B<q; B++){
			if (i != j){ ff[i][j][A][B] = lambda_psc_use; }
			for (k=0; k<numseqs; k++){
				if((MSA[k][pos[i]] == aa[A]) && (MSA[k][pos[j]] == aa[B])){ ff[i][j][A][B] = ff[i][j][A][B] + seqs[k]; }
			}
			ff[i][j][A][B] = ff[i][j][A][B]/lambda_Meff;
		} }
	} }
//Release memory
	free(seqs);

//Report verbose
        if (verbose){printf("\t--Effective number of sequences: %f\n",Meff); }

//Return total number of effective sequences

        return Meff;

}
