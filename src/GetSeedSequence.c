#include "raDI.h"

char **GetSeedSequence(MSA_File,L)
FILE *MSA_File;
int  L;
{
  int  i,j,k,A,aux;
  char *line,**seed;
  char *aa,*mini_aa;

// Request memory
        line = (char*)  calloc(linel,sizeof(char));
	if (L==0 ){
		printf("MSA file is empty\n");
		exit(4);
	}
        seed = (char**) malloc(1*sizeof(char*));
        if (seed == NULL){
		printf("Not enough memory\n");
		exit(1);
	}
        for (i=0; i<1; i++){
		seed[i] = (char*) calloc(L,sizeof(char));
		if (seed[i] == NULL){
			printf("Not enough memory to allocate Seed Sequence\n");
			exit(1);
		}
	}


// Construct seed as a vector of characters
	mini_aa = ".rhkdestnqavlimfwycgp";
        aa      = "-RHKDESTNQAVLIMFWYCGP";
	aux = 0;
	k = -1;
	while (fgets(line,linel,MSA_File) != NULL){
	  if (line[0] == '>'){
		aux = 1;
		i = 0;
		k++;
		continue;
	  }
	  j = 0;
	  if (k>=1) break;
	  if (aux > 0 && k<1){
	    while (line[j] != '\n'){
	       aux = 2;
	       for (A=0; A<21; A++){ if (line[j] == mini_aa[A]){aux = 3;break;}}
	       for (A=0; A<21; A++){ if (line[j] == aa[A]){ aux = 1; break; } }
	       if (aux == 2){	line[j] = '-'; }
	       seed[k][i]=(char) line[j]; 
	       i++;
	       j++;
	    }
	  }

	}

 return seed;

}
