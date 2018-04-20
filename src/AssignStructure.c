#include "raDI.h"

secondary_structure  *AssignStructure(SSA_File,MSA,L,verbose)
FILE *SSA_File;
char **MSA;
int  L,verbose;
{
  int    i,k,n,helices,strands;
  char   *line,**ss_sequence;
  char   dummy_secondary,regular_structure,ss_seq;
  char   reduce_structure();
  //char   dummy_secondary,regular_structure,ss_seq;
  secondary_structure *SSA;

//Request memory
    line         = (char*)  calloc(linel,sizeof(char));
    SSA          = (secondary_structure*) calloc(L,sizeof(secondary_structure));
    ss_sequence  = (char**) malloc(2*sizeof(char*));
    if (ss_sequence == NULL){ printf("It was not possible to allocate memory\n"); exit(1); }
    for (i=0; i<2; i++){
        ss_sequence[i] = (char*) calloc(L,sizeof(char));
        if (ss_sequence[i] == NULL){ 
             printf("Not enough memory to allocate MSA\n"); 
             exit(1);
           }else{
	     memset(ss_sequence[i],'\0',L);
           }

        }


//Read file
    if (SSA_File != NULL){
        k=-1;
        while (!feof(SSA_File) && k<2 ) {
         if (fgets(line,linel,SSA_File) != NULL){ 
            if      (line[0] == '>'){ k++; continue;} 
            else if (line[0] !='\n'){ strncat(ss_sequence[k],line,strlen(line)-1); }
            }	
         }	
        if (verbose){
            printf("\t--Read secondary structure\n");
            printf("\t\t--Seed sequence = %s\n",ss_sequence[0]);
            printf("\t\t--Sec. structure= %s\n",ss_sequence[1]);
            }  
     }else{
        for (i=0;i<L;i++){ 
            ss_sequence[0][i]=MSA[0][i];
            ss_sequence[1][i]='Z';
            }
     }
              

//Check continuity between MSA and SSA seed sequences
    for (i=0;i<L;i++){
      if (MSA[0][i] == '-' && ss_sequence[0][i] != '-'){
         printf("#WARNING: Check differences between MSA seed and the sequence with assigned secondary structure\n");
         printf("\t#Sec. structure= %s\n",ss_sequence[1]);
         printf("\t#MSA  structure= %s\n",MSA[0]);
         printf("\n");
         //exit(1);
         }
     }

//Assign secondary structure
    n=1;
    for (k=0;k<L;k++){
           ss_seq = ss_sequence[1][k];
           regular_structure = reduce_structure(ss_seq);
           if (k==0){dummy_secondary = regular_structure;}
           SSA[k].type = regular_structure ;
           if (dummy_secondary != SSA[k].type ||  SSA[k].type == 'Z'){
               dummy_secondary = SSA[k].type;
               n++;
              }
           SSA[k].order=n;
           }
    if (verbose){
       printf("\t\t--Sec. regular  = ");
       for (k=0;k<L;k++){
         printf("%c",SSA[k].type);
       }
       printf("\n");
    }

// Free memory
    for (i=0; i<2; i++){ free(ss_sequence[i]);}
    free(ss_sequence);
    free(line);

//Report verbose
    if (verbose){
       printf("\t\t--Number of structural fragments: %d\n",SSA[L-1].order);
       helices=0;
       strands=0;
       n=0;
       for (k=0;k<L;k++){
         if (SSA[k].order != n){
           n=SSA[k].order;
           if (SSA[k].type == 'H' ) helices++;
           if (SSA[k].type == 'E' ) strands++;
         }
       }
       printf("\t\t\t--Number of helices: %d\n",helices);
       printf("\t\t\t--Number of strands: %d\n",strands);
    }

// return SSA
    return SSA;
}


char reduce_structure(s)
char s;
{
 char ss;
 if       (s=='H' || s=='h') {ss='H';}
 else if  (s=='E' || s=='e') {ss='E';} 
 else if  (s=='G' || s=='g') {ss='H';} 
 else if  (s=='D' || s=='d') {ss='D';} 
 else if  (s=='-' )          {ss='-';} 
 else if  (s=='Z' )          {ss='Z';} 
 else                        {ss='C';}
 return ss;
}

