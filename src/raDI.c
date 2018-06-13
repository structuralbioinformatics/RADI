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

q = number of possible positions for each position in the alignment
L = length of the MSA
numseqs = number of homolog sequences in the MSA
ra_cluster = reduced alphabet type of clustering
minimum_gap = closest distance in sequence between residues (10 by default)
fragment_size = interval of residues for accepting correct contacts (9 by default, is a 4-turn radius)

*/


int main(int argc, char *argv[])
{

  int i, j, k, A, L, l,  q, numseqs, ra_cluster, help, verbose,swap;
  int minimum_gap, fragment_size, nr;
  int *pos, *xpos;
  int *UseColumns(), *AssignNumberInStructure();
  int **cmap_ca, **cmap_cb, **cmap_min, **contact_map();
  double  Meff, **DI, **MI, **contact_ca, **contact_cb, **contact_min,*time;
  char *aa, *mini_aa, *type, *line, *MSA_filename, *DI_filename, *SSA_filename, *XSA_filename, *PDB_filename;
  char *root,*DI_top_filename, *MI_top_filename, *cmap_filename, *cmap_off_filename, *gnu_filename, *gnu_gif;
  char **MSA, **seed_sequence, **ReadMSA(), **GetSeedSequence();
  PROT protein, readpdb();
  FILE *MSA_File, *DI_File, *SSA_File, *XSA_File, *PDB_File, *DI_top_File, *MI_top_File, *cmap_File, *gnu_File;
  secondary_structure *SSA, *AssignStructure();
  void Results_MI_DI(), MakeFileNames();
  double Calculate_MI_DI();


  //By default, the algorithm runs with all aminoacids
  ra_cluster    = 0;
  verbose       = 0;
  swap          = 0;
  minimum_gap   = 10;
  fragment_size = 9;
  SSA_File      = NULL;
  MSA_File      = NULL;
  PDB_File      = NULL;
  XSA_File      = NULL;


  //Read input
  help=0;
  for (i=0; i<argc; i++){
    if (strcmp(argv[i], "-h") == 0 )  {help=1;}            // help
    if (strcmp(argv[i], "-v") == 0)   {                    // verbose to print log file
      verbose=1;
      printf("#INPUT -v: \tVerbosity is ON\n");
    }
    if (strcmp(argv[i], "-sdis") == 0) {	                // Minimum distance between residues
      sscanf(argv[i+1], "%d", &minimum_gap);
      printf("#INPUT -sdis:  \t%d\n", minimum_gap);
    }
    if (strcmp(argv[i], "-lfrg") == 0) {	                // Fragment length to avoid redundant results
      sscanf(argv[i+1], "%d", &fragment_size);
      printf("#INPUT -lfrg:  \t%d\n", fragment_size);
    }

    if (strcmp(argv[i], "-ra") == 0) {	                  // Selection of Alphabet clustering
      sscanf(argv[i+1], "%d", &ra_cluster);
      printf("#INPUT -ra:  \t%d\n", ra_cluster);
    }
    if (strcmp(argv[i], "-h") == 0) {help=1;}             // print help
    if (strcmp(argv[i], "-msa") == 0){                    // MSA file
      MSA_filename=argv[i+1];
      MSA_File = fopen(MSA_filename, "r");
      printf("#INPUT -msa: \t%s\n", MSA_filename);
      if (MSA_File == NULL){
        printf("--Problems opening the INPUT file %s\n",MSA_filename);
        exit(2);
      }
    }
    if (strcmp(argv[i], "-ssa") == 0){                     // Secondary structure information aligned with target sequence
      SSA_filename=argv[i+1];
      SSA_File = fopen(SSA_filename, "r");
      printf("#INPUT -ssa: \t%s\n", SSA_filename);
      if (SSA_File == NULL){
        printf("--Problems opening the INPUT file %s\n", SSA_filename);
        exit(2);
      }
    }
    if (strcmp(argv[i], "-xsa") == 0){                     // File with sequence of PDB structure aligned with the seed in the MSA file
      XSA_filename=argv[i+1];
      XSA_File = fopen(XSA_filename, "r");
      printf("#INPUT -xsa: \t%s\n", XSA_filename);
      if (XSA_File == NULL){
        printf("--Problems opening the INPUT file %s\n", XSA_filename);
        exit(2);
      }
    }
    if (strcmp(argv[i], "-pdb") == 0){                     // File with PDB structure
      PDB_filename=argv[i+1];
      PDB_File = fopen(PDB_filename, "r");
      printf("#INPUT -pdb: \t%s\n", PDB_filename);
      if (PDB_File == NULL){
        printf("--Problems opening the INPUT file %s\n", PDB_filename);
        exit(2);
      }
    }
    if (strcmp(argv[i], "-o") == 0){                       // Output
      root=argv[i+1];
      printf("#INPUT -plot:\t%s\n", root);
    }
    if (strcmp(argv[i], "-swap") == 0 )  {
      swap=1;
      printf("#INPUT -swap: Randomize sequences is ON\n");
    }
  }

  if (help==1 || (MSA_File == NULL)){
    printf("\
#          Program for Mutual Information and Direct Information between positions of a Multiple Sequence Alignment(MSA)\n\
#          It includes the possibility to check contact-pairs using a PDB file\n\
#          Authors: \n\
#                 Bernat Anton, Mireia Besalú, Gemma de las Cuevas, \n\
#                 Oriol Fornes, Jaume Bonet &  Baldo Oliva\n\
#          Structural Bioinformatics lab (GRIB)\n\
#          Universitat Pompeu Fabra\n\
#          Catalonia\n\
#          Contact: baldo.oliva@upf.edu  (+34 933160509)\n\
#          Last Update: March 2018\n\n\
           \n\t -msa \t\t File with Multiple Sequence Alignments (MSA) \
           \n\t -ssa \t\t File with Secondary Structure (SS) aligned with the seed in the MSA file\
           \n\t -xsa \t\t File with sequence of PDB structure aligned with the seed in the MSA file\
           \n\t -pdb \t\t PDB File of the SINGLE chain structure of the seed\
           \n\t -sdis\t\t Minimum distance in sequence between residue-pairs (default 10) \
           \n\t -lfrg\t\t Fragment length to avoid redundant results (selecting the best pair between two fragments, default 5) \
           \n\t -ra  \t\t Reduced alphabet of AA (default is no-reduction:0)\
           \n\t      \t\t    0 {-}{A}{C}{D}{E}{F}{G}{H}{I}{K}{L}{M}{N}{P}{Q}{R}{S}{T}{V}{W}{Y} \
           \n\t      \t\t    1 {-}{RKH}{DE}{STNQ}{AVLIM}{FWY}{C}{G}{P} \
           \n\t      \t\t    2 {-}{RKHDESTNQC}{AVLIMFWY}{G}{P} \
           \n\t      \t\t    3 {-}{RKHDESTNQCG}{AVLIMFWYP} \
           \n\t -o   \t\t Root name of the outputs with MI, DI and Contact-Map data for GNU plot.\
           \n\t -swap\t\t Swap the AA (or raAA) without affecting gaps \
           \n\t -v   \t\t Verbose \
           \n\t -h   \t\t This help (if MSA file is null help is also shown)\n\n");
    exit(0);
  }


  //Allocate space for filenames
  DI_filename        =(char*)calloc(linel,sizeof(char));
  cmap_filename      =(char*)calloc(linel,sizeof(char));
  cmap_off_filename  =(char*)calloc(linel,sizeof(char));
  MI_top_filename    =(char*)calloc(linel,sizeof(char));
  DI_top_filename    =(char*)calloc(linel,sizeof(char));
  gnu_filename       =(char*)calloc(linel,sizeof(char));
  gnu_gif            =(char*)calloc(linel,sizeof(char));

  //Detect the size of the MSA matrix
  numseqs = 0;
  L = 0;
  line = (char*) calloc(linel, sizeof(char));
  while (fgets(line, MAXSEQ, MSA_File) != NULL){
    if (line[0] == '>'){
      numseqs++;
      if (verbose) printf("#Read sequence %s",line);
      continue;
    }
    if (numseqs == 1){
      i = 0;
      while (line[i] != '\n'){ i++; }
      L = L+i;
    }
  }
  free(line);
  if (verbose) printf("#Sequence length %d\n",L);

  //Get the Seed sequence to define positions with structural information
  rewind(MSA_File);
  seed_sequence=GetSeedSequence(MSA_File,L);

  //Define the aminoacid number in structure for each column-position IF there is no structure positions are pos+1
  if (verbose) printf("-- Assign terciary structure positions\n");
  xpos = AssignNumberInStructure(XSA_File, seed_sequence, L, verbose);
  if (XSA_File != NULL) fclose(XSA_File);

  //Define the secondary structure order for each column-position
  if (verbose) printf("-- Assign secondary structure\n");
  SSA = AssignStructure(SSA_File, seed_sequence, L, verbose);
  if (SSA_File != NULL) fclose(SSA_File);

  //Use the structure of the protein when this is known
  nr=1;
  if (PDB_File != NULL){
    if (verbose) printf("-- Read the PDB file\n");
    protein = readpdb(PDB_File,verbose);
    fclose(PDB_File);
    //Calculate the contact-maps for CB, CA and closest distance
    if (verbose) printf("-- Calculate contact maps\n");
    //Allocate memory
    nr=protein.number_of_res;
    contact_ca =(double**) malloc((nr+1)*sizeof(double*));
    contact_cb =(double**) malloc((nr+1)*sizeof(double*));
    contact_min=(double**) malloc((nr+1)*sizeof(double*));
    for (i=0; i<=nr; i++){ contact_ca[i]=(double*) calloc((nr+1),sizeof(double)); }
    for (i=0; i<=nr; i++){ contact_cb[i]=(double*) calloc((nr+1),sizeof(double)); }
    for (i=0; i<=nr; i++){ contact_min[i]=(double*) calloc((nr+1),sizeof(double)); }
    for (i=0; i<=nr; i++){
      for (j=0; j<=nr; j++){
        contact_ca[i][j]=-1.;
        contact_cb[i][j]=-1.;
        contact_min[i][j]=-1.;
      }
    }
    //Get CMAP
    type="CA";
    cmap_ca  = contact_map(protein, (int)GAP, type, (double)THRESHOLD_CA, contact_ca, verbose);
    type="CB";
    cmap_cb  = contact_map(protein, (int)GAP, type, (double)THRESHOLD_CB, contact_cb, verbose);
    type="MIN";
    cmap_min = contact_map(protein, (int)GAP, type, (double)THRESHOLD_MIN, contact_min, verbose);
  }else{
    nr=L;
    contact_ca =(double**) malloc((nr+1)*sizeof(double*));
    contact_cb =(double**) malloc((nr+1)*sizeof(double*));
    contact_min=(double**) malloc((nr+1)*sizeof(double*));
    for (i=0; i<=nr; i++){ contact_ca[i]=(double*) calloc((nr+1), sizeof(double)); }
    for (i=0; i<=nr; i++){ contact_cb[i]=(double*) calloc((nr+1), sizeof(double)); }
    for (i=0; i<=nr; i++){ contact_min[i]=(double*) calloc((nr+1), sizeof(double)); }
    cmap_ca =(int**) malloc((nr+1)*sizeof(int*));
    cmap_cb =(int**) malloc((nr+1)*sizeof(int*));
    cmap_min=(int**) malloc((nr+1)*sizeof(int*));
    for (i=0; i<=nr; i++){ cmap_ca[i]=(int*) calloc((nr+1), sizeof(int)); }
    for (i=0; i<=nr; i++){ cmap_cb[i]=(int*) calloc((nr+1), sizeof(int)); }
    for (i=0; i<=nr; i++){ cmap_min[i]=(int*) calloc((nr+1), sizeof(int)); }
    for (i=0; i<=nr; i++){
      for (j=0; j<=nr; j++){
        contact_ca[i][j]=-1.;
        contact_cb[i][j]=-1.;
        contact_min[i][j]=-1.;
      }
    }
    for (i=0; i<=nr; i++){
      for (j=i-1;j<=i+1; j++){
        if (j>0 && j<=nr) {
          contact_ca[i][j]=5.0;
          contact_cb[i][j]=5.0;
          contact_min[i][j]=5.0;
        }
      }
    }
    for (i=0; i<=nr; i++){
      for (j=0; j<=nr; j++){
        cmap_ca[i][j]=0;
        cmap_cb[i][j]=0;
        cmap_min[i][j]=0;
      }
    }
    for (i=0; i<=nr; i++){
      for (j=i-1;j<=i+1; j++){
        if (j>0 && j<=nr) {
          cmap_ca[i][j]=1;
          cmap_cb[i][j]=1;
          cmap_min[i][j]=1;
        }
      }
    }
  }


  //Define plot/output names
  MakeFileNames(root, ra_cluster, cmap_filename, cmap_off_filename, MI_top_filename, DI_top_filename, DI_filename, gnu_filename, gnu_gif);
  printf("#OUTPUT GNU:     \t%s\n", gnu_filename);
  printf("#OUTPUT CMAP off:\t%s\n", cmap_off_filename);
  printf("#OUTPUT CMAP:    \t%s\n", cmap_filename);
  printf("#OUTPUT MI Top:  \t%s\n", MI_top_filename);
  printf("#OUTPUT DI Top:  \t%s\n", DI_top_filename);
  printf("#OUTPUT DI :     \t%s\n", DI_filename);

  //Define the mini_aa and aa dictionaries
  mini_aa = ".rhkdestnqavlimfwycgp";
  switch (ra_cluster)
  {
    case  0: aa = "-RHKDESTNQAVLIMFWYCGP"; q=21; break;
    case  1: aa = "012345678"; q=9; break;
    case  2: aa = "01234"; q=5; break;
    case  3: aa = "012"; q=3; break;
    default: aa = "-RHKDESTNQAVLIMFWYCGP"; q=21; break;
  }

  //Construct the MSA matrix
  if (verbose) printf("-- Reading MSA file\n");
  rewind(MSA_File);
  MSA = ReadMSA(MSA_File, numseqs, ra_cluster, L,mini_aa, aa, q, swap, verbose);
  fclose(MSA_File);

  //Define the column-position to use pos[j]=i, with 0<i<L and 0<j<l, being l the real size of columns to study
  if (verbose){printf("-- Select Columns\n");}
  pos = UseColumns(MSA, mini_aa, aa, q, numseqs, L, &l, verbose);

  //Allocate memory for DI and MI
  if (verbose){printf("-- Allocate memory for MI and DI\n");}
  DI  = (double**) malloc(l*sizeof(double*));
  MI  = (double**) malloc(l*sizeof(double*));
  if (DI == NULL || MI == NULL ){
    printf("It was not possible to allocate memory\n");
    exit(1);
  }
  for (i=0; i<l; i++){
    DI[i]  = (double*) calloc(l,sizeof(double));
    MI[i]  = (double*) calloc(l,sizeof(double));
    if (DI[i] == NULL || MI[i] == NULL  ){
      printf("It was not possible to allocate memory\n");
      exit(1);
    }
  }

  //Calculate MI and DI
  time=  (double*) calloc(2,sizeof(double));
  Meff=Calculate_MI_DI(MSA, mini_aa, aa, q, numseqs, L, l, pos, MI, DI, time, verbose);

  //Rank and Print Results
  if (Meff > 1){
    if (verbose) printf("-- Rank and plot MI and DI \n");
    Results_MI_DI(PDB_File, DI_filename, cmap_filename, cmap_off_filename, MI_top_filename, DI_top_filename, gnu_filename,
                  gnu_gif, cmap_ca, contact_ca, cmap_cb, contact_cb, cmap_min, contact_min, nr, l, L, minimum_gap, numseqs,
                  fragment_size, Meff, pos, xpos, SSA, MI, DI, time, verbose);
  }else{
    printf("-- Not enough effective sequences to calculate MI\n");
  }


  //Free memory
  if (verbose) printf("Clean Memory %d\n",nr);
  for (i=0; i<=nr; i++){
    free(cmap_ca[i]);
    free(cmap_cb[i]);
    free(cmap_min[i]);
  }
  for (i=0; i<=nr;i++) free(contact_ca[i]);
  for (i=0; i<=nr;i++) free(contact_cb[i]);
  for (i=0; i<=nr;i++) free(contact_min[i]);
  free(cmap_ca);
  free(cmap_cb);
  free(cmap_min);
  free(contact_ca);
  free(contact_cb);
  free(contact_min);

  for (i=0; i<l; i++){
    free(DI[i]);
    free(MI[i]);
  }
  free(DI);
  free(MI);
  free(time);
  free(pos);
  free(xpos);
  free(cmap_filename);
  free(MI_top_filename);
  free(DI_top_filename);
  free(gnu_filename);
  free(gnu_gif);
  for (i=0; i<L; i++){
    free(MSA[i]);
    if (i==0) free(seed_sequence[i]);
  }
  free(MSA);
  free(seed_sequence);

  //Exit
  printf("Done\n");
  exit(0);

 }


