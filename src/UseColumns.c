#include "raDI.h"

int  *UseColumns( MSA, mini_aa , aa, q, numseqs, L, ll, verbose )
char **MSA;
char *mini_aa, *aa;
int  numseqs, L, *ll,q,verbose;
{
  int   *pos,*pos_new,i,j,k,aux,A,lold,id,l,print_return,use_max,kk,i_meff,lj;
  double Max_id_use, Max_lid, Max_gaps_use, Min_Meff_gaps_use, Meff, lambda_psc_use, lambda_Meff, daux, max_Meff,limit_of_gaps,coverage,min_coverage_use;
  double *gaps, *gaps_check, *seqs;

//Request memory
        if (verbose) {printf("\t--Allocate memory\n");}
  	gaps = (double*) calloc(L,sizeof(double));
	seqs = (double*) calloc(numseqs,sizeof(double));
	pos  = (int*) calloc(L,sizeof(int));
	gaps_check= (double*) calloc(11,sizeof(double));

//Initialize parameters  
        Max_id_use     = (double) Max_id;
        Max_gaps_use   = fmin((double)Max_gaps,(double)Lim_gaps);
	limit_of_gaps  = fmax((double)Max_gaps,(double)Lim_gaps);
	min_coverage_use=(double)min_coverage;
        lambda_psc_use = (double) lambda_psc;
        Min_Meff_gaps_use = (double) Min_Meff_gaps;
        for (i=0; i<L; i++){ pos[i] = i; }
	if (verbose)printf("\t--Using parameters Max_id %f Max_gaps %f limit_of_gaps %f minimum_coverage %f\n",Max_id_use,Max_gaps_use,limit_of_gaps,min_coverage_use);
	for (i=0; i<=10; i++){ 
		gaps_check[i] = Max_gaps_use + (double)i*(limit_of_gaps-Max_gaps_use)/10.0;
	}

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
	lj = j;

        if (verbose) {printf("\t--Number of starting positions  %d\n",lj);}
//Skip columns with a ratio of gaps larger than Max_gaps until Meff > Min_Meff_gaps or use the Maximum value of Meff
	i_meff   = 0;
	use_max  = 0;
	max_Meff = 1.0;
	for (kk=0;kk<=10;kk++){
	    if (verbose){ printf("\t\t--Using Maximum ratio of gaps[%d]:  %f\n",kk,gaps_check[kk]); }
	    lold  = 0;
	    l=lj;
            for (i=0; i<L; i++){ pos[i] = i; }
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
			if (gaps[i] > gaps_check[kk]){ aux = 0; continue; }
			if (aux == 1){ aux = pos[i]; pos[j] = aux; j++; }
		}
		lold = l;
		l = j;
	     }
	     if (verbose){
                printf("\t\t--Effective number of sequences:%f\n",Meff); 
                printf("\t\t--Number of positions:%d\n",l); 
	     }
	     coverage=(double)l/(double)L;
	     if (Meff > Min_Meff_gaps_use && coverage > min_coverage){
	       i_meff   = kk;
	       use_max  = 0;
               if (verbose)printf("\t\t--Minimum effective number of sequences (%f) found with gaps ratio %f covering %f of sequence\n",Meff,gaps_check[i_meff],coverage); 
	       break;
	     }else{
               if (Meff<Min_Meff_gaps_use){
	         if (Meff > max_Meff){
		  max_Meff = Meff;
	          use_max  = 1;
		  i_meff   = kk;
                  if (verbose)printf("\t\t--Maximum effective number of sequences: %f with gaps ratio %f\n",max_Meff,gaps_check[i_meff]); 
	         }
	       }
             } 
	}
	if (use_max){
            if (verbose)printf("\t--Using Maximum effective number of sequences: %f with gaps ratio %f\n",max_Meff,gaps_check[i_meff]); 
	    kk    = i_meff;
            lold  = 0;
	    l     = lj;
            for (i=0; i<L; i++){ pos[i] = i; }
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
			if (gaps[i] > gaps_check[kk]){ aux = 0; continue; }
			if (aux == 1){ aux = pos[i]; pos[j] = aux; j++; }
		}
		lold = l;
		l = j;
	     }
	     coverage=(double)l/(double)L;
	}
//Release memory
	free(gaps);          
	free(seqs);
        free(gaps_check);	

//Size of selected columns      
        if (verbose) {printf("\t--Number of selected positions  %d\n",l);}
        if (verbose) {printf("\t--Number of effective sequences %f\n",Meff);}
        if (verbose) {printf("\t--Coverage of sequence %f\n",coverage);}
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

