#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>


// Internal parameters for execution:
// starting maximum number of gaps per column-position 15%,
// Limitting number of accepted gaps per column position 75%
// maximum %ID 80%
// pseudo-count lamda 0.5
// maximum number of characters of a line 131072
// maximum number of characters of a sequence 10000
// Minimum number of effective sequences to remove gaps Min_Meff_gaps

#define MAXS 	500	/* Maximum number of characters in a line of PDB */
#define MAXATOM 100000	/* Maximum number of atoms in a protein */
#define MAXCHAIN 1   /* Maximum number of CHAINS, must be 1 to analyze a single fold*/
#define Max_id      0.8
#define lambda_psc  0.5
#define Max_gaps    0.15
#define Lim_gaps    0.75
#define min_coverage 0.8
#define Min_Meff_gaps   1000
#define linel       131072
#define MAXSEQ      10000
#define GAP         1   /* gap to calculate contact-maps */
#define THRESHOLD_CA  12.0
#define THRESHOLD_CB   8.0
#define THRESHOLD_MIN  5.0
#define TOPRANK  40




//Structure for secondary structure
typedef struct secondary_structure
{
 char  type;
 int   order;
} secondary_structure;

//Structure for atoms of PDB
typedef struct atom
  {
    char name[5];
    char res[5];
    char res_code[2];
    char chain[5];
    int  res_number;
    char a_res_number[6];
    float x;
    float y;
    float z;
    float bfact;
  } ATOM;

//STructure of protein
typedef struct prot
  {
    int number_of_atoms;
    int number_of_res;
    ATOM atoms_of_prot[MAXATOM];
  } PROT;


// Functions

extern void svdcmp(double**, int, int, double*, double**);
void MaxEnt (int, int, double**, double****, double**);
void MI_Compute (int, int, double**, double****, double**);
void Correction_APC (int, double**,double**);
void Pseudoinverse (int, double**);
char **ReadMSA(FILE*,int, int, int, char*,char*,int,int,int );
int  *UseColumns(char**, char*,char*, int, int, int, int* ,int);
char rnd_aa(int,char*,char);
void sort2(int,double*,double*);
double  EffectiveCounts(char**, char*,char*, int, int, int, int,int*,double**,double****,int );
secondary_structure  *AssignStructure(FILE*,char**,int,int);
PROT readpdb(FILE*,int);
int  **contact_map(PROT,int,char*,double,double**,int);
int RankOrder(int, int, int, int, int*, double**, secondary_structure*, int*, int*, int*, double* ,double*);
double NearestContact(int,int,int,int,double**);
char **GetSeedSequence(FILE*,int);
double Calculate_MI_DI(char**, char*, char*, int, int, int, int, int*, double**, double**, double*, int);
void Results_MI_DI( FILE*,char*,char*,char*,char*,char*,char*,char*,int**,double**,int**,double**,int**,double**,int, int, int, int, int, int , double, int*,int*, secondary_structure*, double**, double**, double*, int,int);
void MakeFileNames(char*,int,char*,char*,char*,char*,char*,char*,char*);
