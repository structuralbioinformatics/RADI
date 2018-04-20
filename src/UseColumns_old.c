#include "raDI.h"

int  *UseColumns( MSA, mini_aa , aa, q, numseqs, L, ll, verbose )
char **MSA;
char *mini_aa, *aa;
int  numseqs, L, *ll,q,verbose;
{
  int   *pos,i,j,k,aux,A,lold,id,l,print_return;
  double Max_id_use, Max_lid, Max_gaps_use, Min_Meff_gaps_use, Meff, lambda_psc_use, lambda_Meff, daux;
  double *gaps, *seqs;

//Request memory
        if (verbose) {printf("\t--Allocate memory\n");}
  	gaps = (double*) calloc(L,sizeof(double));
	seqs = (double*) calloc(numseqs,sizeof(double));
	pos  = (int*) calloc(L,sizeof(int));

//Initialize parameters  
        Max_id_use     = (double) Max_id;
        Max_gaps_use   = (double) Max_gaps;
        lambda_psc_use = (double) lambda_psc;
        Min_Meff_gaps_use = (double) Min_Meff_gaps;
        for (i=0; i<L; i++){ pos[i] = i; }

//Skip columns with lowercase or gap (if the first sequence has one, all other sequences have it)
        j = 0;
        for (i=0; i<L; i++){
		aux = 0;
		if (MSA[0][i] == aa[0]){ continue; }
		for (A=0; A<q; A++){
		 if (MSA[0][i] == mini_aa[A]){ aux = 1; break; }
		}
		if (aux == 0){ pos[j] = i; j++; }
	}
	l = j;

        if (verbose) {printf("\t--Number of starting positions  %d\n",l);}
//Skip columns with a ratio of gaps larger than Max_gaps when the number of effective sequences is high enough (Min_Meff_gaps)
	lold = 0;
	while (l != lold){
		Max_lid = Max_id_use*l;
		Meff = 0;
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
		daux           = lambda_psc_use/q;
		for (i=0; i<l; i++){
			gaps[i] = daux;
			for (k=0; k<numseqs; k++){
				if (MSA[k][pos[i]] == aa[0]){ gaps[i] = gaps[i] + seqs[k]; }
			}
			gaps[i] = gaps[i]/lambda_Meff;
		}
		j = 0;
		for (i=0; i<l; i++){
			aux = 1;
			if (gaps[i] > Max_gaps_use && Meff > Min_Meff_gaps_use ){ aux = 0; continue; }
			if (aux == 1){ aux = pos[i]; pos[j] = aux; j++; }
		}
		lold = l;
		l = j;
	}
//Release memory
	free(gaps);          
	free(seqs);    

//Size of selected columns      
        if (verbose) {printf("\t--Number of selected positions  %d\n",l);}
        if (verbose) {printf("\t--Number of effective sequences %f\n",Meff);}
        *ll = l;
	print_return=0;
	if (verbose){
         printf("\t--Selected positions:\n\t\t");
	 for (i=0;i<l;i++){
		 print_return++;
		 printf("%3d,",pos[i]);
		 if (print_return>20){printf("\n\t\t");print_return=0; }
	 }
	 printf("\n");
	}

        return pos;

}

