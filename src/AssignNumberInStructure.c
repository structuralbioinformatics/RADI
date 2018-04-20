#include "raDI.h"

int  *AssignNumberInStructure(XSA_File,MSA,L,verbose)
FILE *XSA_File;
char **MSA;
int  L,verbose;
{
  int    i,k,n,print_return,*xpos;
  char   *line,**sequence;

//Request memory
    line         = (char*)  calloc(linel,sizeof(char));
    sequence     = (char**) malloc(2*sizeof(char*));
    xpos         = (int*) calloc(L,sizeof(int));

    if (sequence == NULL){ printf("It was not possible to allocate memory\n"); exit(1); }
    for (i=0; i<2; i++){
        sequence[i] = (char*) calloc(L,sizeof(char));
        if (sequence[i] == NULL){ 
             printf("Not enough memory to allocate MSA\n"); 
             exit(1);
           }else{
	     memset(sequence[i],'\0',L);
           }

        }


//Read file
    if (XSA_File != NULL){
        k=-1;
        while (!feof(XSA_File) && k<2 ) {
         if (fgets(line,linel,XSA_File) != NULL){ 
            if      (line[0] == '>'){ k++; continue;} 
            else if (line[0] !='\n'){ strncat(sequence[k],line,strlen(line)-1); }
            }	
         }	
        if (verbose){
            printf("\t--Read sequence with structure\n");
            printf("\t\t--Seed sequence = %s\n",sequence[0]);
            printf("\t\t--Seq. structure= %s\n",sequence[1]);
            }  
     }else{
	for (i=0;i<L;i++){xpos[i]=i+1;}
     }
              

//Check continuity between MSA and SSA seed sequences
    for (i=0;i<L;i++){
      if (MSA[0][i] == '-' && sequence[0][i] != '-' ){
         printf("#WARNING: Check differences between MSA seed and the sequence with assigned secondary structure\n");
         printf("\t#Seq. structure= %s\n",sequence[1]);
         printf("\t#MSA  structure= %s\n",MSA[0]);
         printf("\n");
         //exit(1);
         }
     }

//Assign  structure
    n=1;
    for (i=0;i<L;i++){
           if ( MSA[0][i] == '-'  && sequence[1][i] != '-') {n++;}
           if ( MSA[0][i] != '-'  && sequence[1][i] != '-') {xpos[i]=n;n++;}
           if ( MSA[0][i] != '-'  && sequence[1][i] == '-') {xpos[i]=0;}
           }

// Free memory
    for (i=0; i<2; i++){ free(sequence[i]);}
    free(sequence);
    free(line);

//Report verbose
    if (verbose){
         print_return=0;
         printf("\t--PDB positions:\n\t\t");
	 for (i=0;i<L;i++){
		 print_return++;
		 printf("%3d,",xpos[i]);
		 if (print_return>20){printf("\n\t\t");print_return=0; }
	 }
	 printf("\n");

    }
 

// return xpos
    return xpos;
}



