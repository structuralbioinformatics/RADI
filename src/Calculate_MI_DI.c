#include "raDI.h"

/*


#          Program for Mutual Information and Direct Information between positions of a Multiple Sequence Alignment(MSA)
#          It includes the possibility to check contact-pairs using a PDB file
#          Authors:
#          Bernat Anton, Mireia Besalú, Gemma de las Cuevas, Oriol Fornes, Jaume Bonet &  Baldo Oliva
#          Structural Bioinformatics lab (GRIB)
#          Universitat Pompeu Fabra
#          Catalonia\n\
#          Contact: baldo.oliva@upf.edu  (+34 933160509)
#          Last Update: March 2018

*/


double  Calculate_MI_DI(MSA,mini_aa,aa,q,numseqs,L,l,pos,MI,DI,time,verbose)
char **MSA,*aa, *mini_aa;
int    q,numseqs,L,l,verbose,*pos;
double **DI, **MI, *time;
{
 int  i,j,A;
 double **MIo,**DIo,**f, ****ff;
 double  Meff,EffectiveCounts();
 clock_t t,t_start,t_end,t_init;

//Starting time on calculations
        t_init = clock();
//Allocate amino-acid effective frequencies per position
        if (verbose){printf("-- Allocate memory for residue-frequencies\n");}
	f = (double**) malloc(l*sizeof(double*));
	ff = (double****) malloc(l*sizeof(double***)); 
        MIo = (double**) malloc(l*sizeof(double*));
        DIo = (double**) malloc(l*sizeof(double*));
	if (f == NULL || ff == NULL || MIo == NULL || DIo == NULL ){ printf("It was not possible to allocate memory\n"); exit(1); }
	for (i=0; i<l; i++){
		f[i] = (double*) calloc(q,sizeof(double));
		ff[i] = (double***) malloc(l*sizeof(double**));
                MIo[i]  = (double*) calloc(l,sizeof(double));
                DIo[i]  = (double*) calloc(l,sizeof(double));
		if (f[i] == NULL || ff[i] == NULL || MIo[i] == NULL || DIo[i] == NULL){ printf("It was not possible to allocate memory\n"); exit(1); }
		for (j=0; j<l; j++){
			ff[i][j] = (double**) malloc(q*sizeof(double*));
			if (ff[i][j] == NULL){ printf("It was not possible to allocate memory\n"); exit(1); }
			for (A=0; A<q; A++){
				ff[i][j][A] = (double*) calloc(q,sizeof(double));
				if (ff[i][j][A] == NULL){ printf("It was not possible to allocate memory\n"); exit(1); }
			}				
		}
	}


//Extract amino-acid effective counts per position
        if (verbose){printf("-- Use effective counts and compute their frequencies\n");}
        Meff = EffectiveCounts(MSA,mini_aa,aa,q,numseqs,L,l,pos,f,ff,verbose);
        t_end   = clock();
        t = t_end - t_init;
        if (verbose) {printf("Effective count execution time t= %f seconds\n", ((float)t)/CLOCKS_PER_SEC);}

//Compute MI and correction APC
        t_start = clock();
        MI_Compute(l, q, f, ff, MIo);
        Correction_APC(l,MIo,MI);
        t_end   = clock();
        t = t_end - t_start;
        if (verbose) {printf("MI execution time t= %f seconds\n", ((float)t)/CLOCKS_PER_SEC);}
        time[0]=((double)t)/CLOCKS_PER_SEC;

//Compute Direct Information and correction APC;
        t_start = clock();
        MaxEnt (l, q, f, ff, DIo);
        Correction_APC(l,DIo,DI);
        t_end   = clock();
        t = t_end - t_start;
        if (verbose) {printf("DI execution time t= %f seconds\n", ((float)t)/CLOCKS_PER_SEC);}
        time[1]=((double)t)/CLOCKS_PER_SEC;
        t =  clock() - t_init;
        if (verbose) {printf("Total time Calculation t= %f seconds\n", ((float)t)/CLOCKS_PER_SEC);}

//Free memory
	for (i=0; i<l; i++){
          for (j=0; j<l; j++){
            for (A=0; A<q; A++){ free(ff[i][j][A]);}
            free(ff[i][j]);
            }
          free(f[i]);
          free(ff[i]);
          }
        free(f);
        free(ff);

        return Meff;

}
