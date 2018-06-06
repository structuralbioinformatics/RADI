#include "raDI.h"

char **ReadMSA( MSA_File, numseqs, ra_cluster, L, mini_aa, aa, q, swap, verbose )
FILE *MSA_File;
char *mini_aa,*aa;
int  numseqs, L, ra_cluster, q, swap, verbose;
{
  int i, j, k, A, aux;
  char *line, **MSA;

  if (verbose) {
    printf("\t--Data: L %d q %d AA %s \n", L, q, aa);
  }

  // Request memory
  line = (char*) malloc(linel*sizeof(char));
  if (L==0 || numseqs == 0){
    printf("MSA file is empty\n");
    exit(4);
  }
  MSA = (char**) malloc(numseqs*sizeof(char*));
  if (MSA == NULL){
    printf("Not enough memory\n");
    exit(1);
  }
  for (i=0; i<numseqs; i++){
    MSA[i] = (char*) calloc(L,sizeof(char));
    if (MSA[i] == NULL){
      printf("Not enough memory to allocate MSA\n");
      exit(1);
    }
  }

  // Construct MSA as a matrix of characters
  k = -1;
  aux = 0;
  while (fgets(line,linel,MSA_File) != NULL){
    if (line[0] == '>'){
      aux = 1;
      i = 0;
      k++;
      continue;
    }
    j = 0;
    if (aux > 0) {
      while (line[j] != '\n') {
        for (A=0; A<21; A++) {
          if (line[j] == mini_aa[A]) {
            aux = 3; break;
          }
        }
        switch (ra_cluster) {
          case 0:
            aux = 2;
            for (A=0; A<q; A++) {
              if (line[j] == aa[A]) {
                aux = 1;
                break;
              }
            }
            if (aux == 2){
              line[j] = '-'; //non-std aa into gaps
            }
          break;
          case 1:
            aux = 2;
            if (line[j] == '-') { line[j] = '0'; aux = 1;}
            else if (line[j] == 'R' || line[j] == 'H' || line[j] == 'K') { line[j] = '1'; aux = 1; }
            else if (line[j] == 'D' || line[j] == 'E') { line[j] = '2'; aux = 1; }
            else if (line[j] == 'S' || line[j] == 'T' || line[j] == 'N' || line[j] == 'Q') {	line[j] = '3'; aux = 1; }
            else if (line[j] == 'A' || line[j] == 'V' || line[j] ==	'L' || line[j] == 'I' || line[j] == 'M') { line[j] = '4'; aux = 1; }
            else if (line[j] == 'F' || line[j] == 'W' || line[j] == 'Y') { line[j] = '5'; aux = 1; }
            else if (line[j] == 'C') { line[j] = '6'; aux = 1; }
            else if (line[j] == 'G') { line[j] = '7'; aux = 1; }
            else if (line[j] == 'P') { line[j] = '8'; aux = 1; }
            if (aux == 2){ //Control for non-std aa
              line[j] = '0'; //non-std aa into gaps
            }
          break;
          case 2:
            aux = 2;
            if (line[j] == '-'){ line[j] = '0'; aux = 1; }
            else if (line[j] == 'R' || line[j] == 'H' || line[j] == 'K' || line[j] == 'D' || line[j] == 'E' || line[j] == 'S' || line[j] == 'T' || line[j] == 'N' || line[j] == 'Q' || line[j] == 'C' ){line[j] = '1'; aux = 1; }
            else if (line[j] == 'A' || line[j] == 'V' || line[j] ==	'L' || line[j] == 'I' || line[j] == 'M' || line[j] == 'F' || line[j] == 'W' || line[j] == 'Y'){ line[j] = '2'; aux = 1; }
            else if (line[j] == 'G'){ line[j] = '3'; aux = 1; }
            else if (line[j] == 'P'){ line[j] = '4'; aux = 1; }
            if (aux == 2){ //Control for non-std aa
              line[j] = '0'; //non-std aa into gaps
            }
          break;
          case 3:
            aux = 2;
            if (line[j] == '-'){ line[j] = '0'; aux = 1; }
            else if (line[j] == 'R' || line[j] == 'H' || line[j] == 'K' || line[j] == 'D' || line[j] == 'E' || line[j] == 'S' || line[j] == 'T' || line[j] == 'N' || line[j] == 'Q' || line[j] == 'C' || line[j] == 'G' ){line[j] = '1'; aux = 1; }
            else if (line[j] == 'A' || line[j] == 'V' || line[j] ==	'L' || line[j] == 'I' || line[j] == 'M' || line[j] == 'F' || line[j] == 'W' || line[j] == 'Y' || line[j] == 'P' ){ line[j] = '2'; aux = 1; }
            if (aux == 2){ //Control for non-std aa
              line[j] = '0'; //non-std aa into gaps
            }
          break;
          default:
            aux = 2;
            for (A=0; A<q; A++){ if (line[j] == aa[A]){ aux = 1; break; } }
            if (aux == 2){
              line[j] = '-'; //non-std aa into gaps
            }
          break;
        }
        if (swap == 1) {
          if (line[j] == '-' || line[j] == '0' || aux == 3){
             MSA[k][i] = (char) line[j];
           }else{
             MSA[k][i] = rnd_aa(q,aa,(char) line[j]);
           }
        } else {
           MSA[k][i] = (char) line[j];
        }
        i++;
        j++;
      }
    }
  }
  free(line);

  return MSA;
}


char rnd_aa(int q, char *aa, char test) {
  int   i, k, r, qq;
  char  *aar, result;
  time_t time();

  srand(time(NULL));

  aar = (char*) malloc(q*sizeof(char));
  qq = 0;
  for (i=1; i<q; i++){
    if (aa[i] != test){
      aar[qq] = aa[i];
      qq++;
    }
  }

  r = rand();
  k = r % (qq);
  result = aar[k];
  free( aar);

  return result;
}


