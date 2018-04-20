#include "raDI.h"

/*
q = number of possible positions for each position in the alignment
L = length of the MSA
numseqs = number of homolog sequences in the MSA
*/


int main(int argc, char *argv[])
{

	int i, j, k, A, B, id, aux, aux2, L, l, lold, q, numseqs, ra_cluster, help, verbose,n_index,ss_i,ss_j,rank,index_mi,index_di,swap;
        int minimum_gap,fragment_size,print_return,max_rank_mi,max_rank_di,nr,i_mi,j_mi,i_di,j_di;
	int tp_mi_ca,tp_mi_cb,tp_mi_min,tp_mi_sum,tp_di_ca,tp_di_cb,tp_di_min,tp_di_sum;
	int fp_mi_ca,fp_mi_cb,fp_mi_min,fp_mi_sum,fp_di_ca,fp_di_cb,fp_di_min,fp_di_sum;
	int *pos,*xpos,*di_rank,*mi_rank,*j_index,*i_index;
        int *UseColumns(),*AssignNumberInStructure(),RankOrder();
	int **cmap_ca,**cmap_cb,**cmap_min,**cmap_sum,**contact_map();
        double *mi_index,*di_index,*DI_order,*MI_order;
	double  Meff,EffectiveCounts();
	double  NearestContact();
	double **DI, **MI, **f, ****ff, **contact_ca,**contact_cb,**contact_min;
	char *aa, *mini_aa, *type, *line, *ss, *MSA_filename, *DI_filename, *SSA_filename, *XSA_filename,*PDB_filename;
	char *root,*DI_top_filename,*MI_top_filename,*cmap_filename,*gnu_filename,*gnu_gif;
	char **MSA, **ReadMSA();
	PROT protein,readpdb();
	clock_t t,t_start,t_end,t_init;
	FILE *MSA_File, *DI_File, *SSA_File, *XSA_File,*PDB_File,*DI_top_File,*MI_top_File,*cmap_File,*gnu_File;
        secondary_structure *SSA,*AssignStructure();
        void sort2();

//Starting time
        t_init = clock();

//By default, the algorithm runs with all aminoacids 	
	ra_cluster = 0;
        SSA_File=NULL;
        MSA_File=NULL;
	PDB_File=NULL;
	XSA_File=NULL;
        MI_top_File = NULL;
        DI_top_File = NULL;
        cmap_File   = NULL;
        gnu_File    = NULL;
	verbose=0;
        swap=0;
	minimum_gap=5;
	fragment_size=1;

//Read input
        help=0;
        for (i=0;i<argc;i++){
           if (strcmp(argv[i], "-h") == 0 )  {help=1;}            // help
     	   if (strcmp(argv[i], "-v") == 0)   {                    // verbose to print log file
              verbose=1; 
              printf("#INPUT -v: \tVerbosity is ON\n");
           }          
           if (strcmp(argv[i], "-sdis") == 0) {	                  // Minimum distance between residues
              sscanf(argv[i+1], "%d", &minimum_gap);
              printf("#INPUT -sdis:  \t%d\n",minimum_gap);
           }
           if (strcmp(argv[i], "-lfrg") == 0) {	                  // Fragment length to avoid redundant results
              sscanf(argv[i+1], "%d", &fragment_size);
              printf("#INPUT -lfrg:  \t%d\n",fragment_size);
           }

           if (strcmp(argv[i], "-ra") == 0) {	                  // Selection of Alphabet clustering
              sscanf(argv[i+1], "%d", &ra_cluster);
              printf("#INPUT -ra:  \t%d\n",ra_cluster);
           }
     	   if (strcmp(argv[i], "-h") == 0) {help=1;}		  // print help
	   if (strcmp(argv[i], "-msa") == 0){                     // MSA file
	      MSA_filename=argv[i+1];
	      MSA_File = fopen(MSA_filename, "r");
              printf("#INPUT -msa: \t%s\n",MSA_filename);
              if (MSA_File == NULL){
		printf("--Problems opening the INPUT file %s\n",MSA_filename);
		exit(2);
                }
	   }
	   if (strcmp(argv[i], "-ssa") == 0){                     // Secondary structure information aligned with target sequence
	      SSA_filename=argv[i+1];
	      SSA_File = fopen(SSA_filename, "r");
              printf("#INPUT -ssa: \t%s\n",SSA_filename);
              if (SSA_File == NULL){
		printf("--Problems opening the INPUT file %s\n",SSA_filename);
		exit(2);
                }
	   }
	   if (strcmp(argv[i], "-xsa") == 0){                     // File with sequence of PDB structure aligned with the seed in the MSA file
	      XSA_filename=argv[i+1];
	      XSA_File = fopen(XSA_filename, "r");
              printf("#INPUT -xsa: \t%s\n",XSA_filename);
              if (XSA_File == NULL){
		printf("--Problems opening the INPUT file %s\n",XSA_filename);
		exit(2);
                }
	   }
	   if (strcmp(argv[i], "-pdb") == 0){                     // File with PDB structure
	      PDB_filename=argv[i+1];
	      PDB_File = fopen(PDB_filename, "r");
              printf("#INPUT -pdb: \t%s\n",PDB_filename);
              if (PDB_File == NULL){
		printf("--Problems opening the INPUT file %s\n",PDB_filename);
		exit(2);
                }
	   }
	   if (strcmp(argv[i], "-o") == 0){                       // Output
	      DI_filename=argv[i+1];
	      DI_File = fopen(DI_filename, "w");
              printf("#INPUT -o:   \t%s\n",DI_filename);
              if (DI_File == NULL){
		printf("--Problems opening the OUTPUT file %s\n",DI_filename);
		exit(2);
                }
	   }
	   if (strcmp(argv[i], "-plot") == 0){                       // Output
              root=argv[i+1];
              printf("#INPUT -plot:\t%s\n",root);
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
#          Bernat Anton & Baldo Oliva\n\
#          Structural Bioinformatics lab (GRIB)\n\
#          Universitat Pompeu Fabra\n\
#          Catalonia\n\
#          Contact: baldo.oliva@upf.edu  (+34 933160509)\n\
#          Last Update: February 2018\n\
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
           \n\t -o   \t\t Output file with MI and DI results\
           \n\t -plot\t\t Root name of the outputs with MI, DI and Contact-Map data for GNU plot.\
           \n\t -swap\t\t Swap the AA (or raAA) without affecting gaps \
           \n\t -v   \t\t Verbose \
	   \n\t -h   \t\t This help (if MSA file is null help is also shown)\n\n");
	   exit(0);
	}

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

//Define plot names
        cmap_filename=(char*)calloc(linel,sizeof(char));
        MI_top_filename=(char*)calloc(linel,sizeof(char));
        DI_top_filename=(char*)calloc(linel,sizeof(char));
        gnu_filename=(char*)calloc(linel,sizeof(char));
        gnu_gif=(char*)calloc(linel,sizeof(char));
        sprintf(cmap_filename,"%s_ra%1d_cmap.dat",root,ra_cluster);
        sprintf(MI_top_filename,"%s_ra%1d_mi.dat",root,ra_cluster);
        sprintf(DI_top_filename,"%s_ra%1d_di.dat",root,ra_cluster);
        sprintf(gnu_filename,"%s_ra%1d_gnu.dat",root,ra_cluster);
        sprintf(gnu_gif,"%s_ra%1d_cmap.gif",root,ra_cluster);
	DI_top_File = fopen(DI_top_filename, "w");
	MI_top_File = fopen(MI_top_filename, "w");
	cmap_File   = fopen(cmap_filename, "w");
        gnu_File    = fopen(gnu_filename, "w");
        printf("#INPUT -plot GNU:    \t%s\n",gnu_filename);
        printf("#INPUT -plot CMAP:   \t%s\n",cmap_filename);
        printf("#INPUT -plot MI Top: \t%s\n",MI_top_filename);
        printf("#INPUT -plot DI Top: \t%s\n",DI_top_filename);
        if (DI_top_File == NULL){
		printf("--Problems opening the OUTPUT file %s\n",DI_top_filename);
		exit(2);
                }
        if (MI_top_File == NULL){
		printf("--Problems opening the OUTPUT file %s\n",MI_top_filename);
		exit(2);
                }
        if (cmap_File == NULL){
		printf("--Problems opening the OUTPUT file %s\n",cmap_filename);
		exit(2);
                }
        if (gnu_File == NULL){
		printf("--Problems opening the OUTPUT file %s\n",gnu_filename);
		exit(2);
                }

 
//Detect the size of the MSA matrix
        numseqs = 0;
        L = 0;
        line = (char*) calloc(linel,sizeof(char));
        while (fgets(line,MAXSEQ,MSA_File) != NULL){
              
              if (line[0] == '>'){ numseqs++; if (verbose) printf("#Read sequence %s",line);continue; }
              if (numseqs == 1){
			i = 0;
			while (line[i] != '\n'){ i++; }
			L = L+i;
              }
        }
        rewind(MSA_File);
        free(line);
        if (verbose) printf("#Sequence length %d\n",L);

//Construct the MSA matrix
        if (verbose) {printf("-- Reading MSA file\n");}
        MSA = ReadMSA(MSA_File, numseqs, ra_cluster, L,mini_aa,aa,q,swap,verbose);
        fclose(MSA_File);

//Define the column-position to use pos[j]=i, with 0<i<L and 0<j<l, being l the real size of columns to study
        if (verbose){printf("-- Select Columns\n");}
        pos = UseColumns(MSA,mini_aa,aa,q,numseqs, L, &l,verbose);
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


//Use the structure of the protein when this is known
        if (PDB_File != NULL){
          if (verbose){printf("-- Read the PDB file\n");}
          protein = readpdb(PDB_File,verbose);
          fclose(PDB_File);
          //Calculate the contact-maps for CB, CA and closest distance
          if (verbose){printf("-- Calculate contact maps\n");}
	  //Allocate memory
	  nr=protein.number_of_res;
          contact_ca =(double**) malloc((nr+1)*sizeof(double*));
          contact_cb =(double**) malloc((nr+1)*sizeof(double*));
          contact_min=(double**) malloc((nr+1)*sizeof(double*));
          for (i=0; i<=nr; i++){contact_ca[i]=(double*) calloc((nr+1),sizeof(double));}
          for (i=0; i<=nr; i++){contact_cb[i]=(double*) calloc((nr+1),sizeof(double));}
          for (i=0; i<=nr; i++){contact_min[i]=(double*) calloc((nr+1),sizeof(double));}
          for (i=0; i<=nr; i++){for (j=0; j<=nr; j++){contact_ca[i][j]=-1.;contact_cb[i][j]=-1.;contact_min[i][j]=-1.;}}
	  //Get CMAP
	  type="CA";
          cmap_ca  = contact_map(protein,(int)GAP,type,(double)THRESHOLD_CA,contact_ca,verbose);
	  type="CB";
          cmap_cb  = contact_map(protein,(int)GAP,type,(double)THRESHOLD_CB,contact_cb,verbose);
	  type="MIN";
          cmap_min = contact_map(protein,(int)GAP,type,(double)THRESHOLD_MIN,contact_min,verbose);
        }

//Define the aminoacid number in structure for each column-position IF there is no structure positions are pos+1
        if (verbose){printf("-- Assign terciary structure positions\n");}
        xpos = AssignNumberInStructure(XSA_File,MSA,L,verbose);
        if (XSA_File != NULL){fclose(XSA_File);}

//Define the secondary structure order for each column-position
        if (verbose){printf("-- Assign secondary structure\n");}
        SSA = AssignStructure(SSA_File,MSA,L,verbose);
        if (SSA_File != NULL){fclose(SSA_File);}

//Allocate amino-acid effective frequencies per position
        if (verbose){printf("-- Allocate memory for residue-frequencies\n");}
	f = (double**) malloc(l*sizeof(double*));
	ff = (double****) malloc(l*sizeof(double***)); 
	if (f == NULL || ff == NULL){ printf("It was not possible to allocate memory\n"); exit(1); }
	for (i=0; i<l; i++){
		f[i] = (double*) calloc(q,sizeof(double));
		ff[i] = (double***) malloc(l*sizeof(double**));
		if (f[i] == NULL || ff[i] == NULL){ printf("It was not possible to allocate memory\n"); exit(1); }
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

//Allocate memory for DI and MI
	DI  = (double**) malloc(l*sizeof(double*));
	MI  = (double**) malloc(l*sizeof(double*));
	if (DI == NULL || MI == NULL ){ printf("It was not possible to allocate memory\n"); exit(1); }
	for (i=0; i<l; i++){
		DI[i]  = (double*) calloc(l,sizeof(double));
		MI[i]  = (double*) calloc(l,sizeof(double));
		if (DI[i] == NULL || MI[i] == NULL  ){ printf("It was not possible to allocate memory\n"); exit(1); }
	}

//Compute MI
        t_start = clock();
        MI_Compute(l, q, f, ff, MI);
        t_end   = clock();
        t = t_end - t_start;
        if (verbose) {printf("MI execution time t= %f seconds\n", ((float)t)/CLOCKS_PER_SEC);}

//Compute Direct Information;
        t_start = clock();
        MaxEnt (l, q, f, ff, DI);
        t_end   = clock();
        t = t_end - t_start;
        if (verbose) {printf("DI execution time t= %f seconds\n", ((float)t)/CLOCKS_PER_SEC);}

//Reorder by MI and DI
//      Calculate size
        n_index = 1;
        for (i=0; i<l; i++){
        for (j=i+1; j<l; j++){
            ss_i  = SSA[pos[i]].order;
            ss_j  = SSA[pos[j]].order;
            if (ss_i != ss_j  && abs(pos[i]-pos[j])>minimum_gap ){  n_index++; }
            }}

//      Request memory
        mi_rank   = (int*) calloc(n_index+1,sizeof(int));
        di_rank   = (int*) calloc(n_index+1,sizeof(int));
	DI_order  = (double*) calloc(n_index+1,sizeof(double));
	MI_order  = (double*) calloc(n_index+1,sizeof(double));
        mi_index  = (double*) calloc(n_index+1,sizeof(double));
        di_index  = (double*) calloc(n_index+1,sizeof(double));
	j_index   = (int*) calloc(n_index+1,sizeof(int));
	i_index   = (int*) calloc(n_index+1,sizeof(int));

//      Rank Order MI (small@1 -> largest@(index-1))
        max_rank_mi=RankOrder(minimum_gap,fragment_size,l,n_index,pos,MI,SSA,mi_rank,i_index,j_index,mi_index,MI_order);

//      Rank Order DI (small@1 -> largest@(index-1))
        max_rank_di=RankOrder(minimum_gap,fragment_size,l,n_index,pos,DI,SSA,di_rank,i_index,j_index,di_index,DI_order);

        if (PDB_File!=NULL){
//Analyses
// TP vs FP (%PPV)
         tp_mi_ca =0;
         tp_mi_cb =0;
         tp_mi_min=0;
         tp_mi_sum=0;
         fp_mi_ca =0;
         fp_mi_cb =0;
         fp_mi_min=0;
         fp_mi_sum=0;
         tp_di_ca =0;
         tp_di_cb =0;
         tp_di_min=0;
         tp_di_sum=0;
         fp_di_ca =0;
         fp_di_cb =0;
         fp_di_min=0;
         fp_di_sum=0;
         for (rank=1;rank<=TOPRANK;rank++){
            index_di = (int)di_index[di_rank[rank]];
            index_mi = (int)mi_index[mi_rank[rank]];
            i_mi     = i_index[index_mi];
            j_mi     = j_index[index_mi];
            i_di     = i_index[index_di];
            j_di     = j_index[index_di];
            if (NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_ca) < (double)THRESHOLD_CA){tp_mi_ca++;}else{fp_mi_ca++;}
            if (NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_ca) < (double)THRESHOLD_CA){tp_di_ca++;}else{fp_di_ca++;}
            if (NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_cb) < (double)THRESHOLD_CB){tp_mi_cb++;}else{fp_mi_cb++;}
            if (NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_cb) < (double)THRESHOLD_CB){tp_di_cb++;}else{fp_di_cb++;}
            if (NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_min) < (double)THRESHOLD_MIN){tp_mi_min++;}else{fp_mi_min++;}
            if (NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_min) < (double)THRESHOLD_MIN){tp_di_min++;}else{fp_di_min++;}
	    if (   NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_ca) < (double)THRESHOLD_CA
                || NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_cb) < (double)THRESHOLD_CB
                || NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_min) < (double)THRESHOLD_MIN){tp_mi_sum++;}else{fp_mi_sum++;}
	    if (   NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_ca) < (double)THRESHOLD_CA
                || NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_cb) < (double)THRESHOLD_CB
                || NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_min) < (double)THRESHOLD_MIN){tp_di_sum++;}else{fp_di_sum++;}
         }

//Plot
         cmap_sum=(int**) malloc((nr+1)*sizeof(int*));
         for (i=0; i<=nr; i++){cmap_sum[i]=(int*) calloc((nr+1),sizeof(int));}
         for (i=0; i<=nr; i++){cmap_sum[i][j]=0;}
         for (i=0; i<=nr; i++){
         for (j=0; j<=nr; j++){
             if (cmap_ca[i][j] == 1 || cmap_cb[i][j] == 1  || cmap_min[i][j] == 1) cmap_sum[i][j]=1;
         }}
         fprintf(cmap_File,"#pos1\tpos2\tradius\n");
         fprintf(MI_top_File,"#pos1\tpos2\n");
         fprintf(DI_top_File,"#pos1\tpos2\n");
         for (i=1; i<=nr; i++){
         for (j=i+1; j<=nr; j++){
         if  (cmap_sum[i][j]) {
             fprintf(cmap_File,"%3d\t%3d\t%3d\n",i,j,2);
             fprintf(cmap_File,"%3d\t%3d\t%3d\n",j,i,2);
         }}}
         for (rank=1;rank<=TOPRANK;rank++){
            index_di = (int)di_index[di_rank[rank]];
            index_mi = (int)mi_index[mi_rank[rank]];
            i_mi     = i_index[index_mi];
            j_mi     = j_index[index_mi];
            i_di     = i_index[index_di];
            j_di     = j_index[index_di];
            if (xpos[pos[i_mi]]<xpos[pos[j_mi]]){fprintf(MI_top_File,"%3d\t%3d\n",xpos[pos[i_mi]],xpos[pos[j_mi]]);}else{fprintf(MI_top_File,"%3d\t%3d\n",xpos[pos[j_mi]],xpos[pos[i_mi]]);}
            if (xpos[pos[i_mi]]>xpos[pos[j_mi]]){fprintf(DI_top_File,"%3d\t%3d\n",xpos[pos[i_di]],xpos[pos[j_di]]);}else{fprintf(DI_top_File,"%3d\t%3d\n",xpos[pos[j_di]],xpos[pos[i_di]]);}
         }
         fprintf(gnu_File,"set term gif\n");
         fprintf(gnu_File,"set xrange [0:]\n");
         fprintf(gnu_File,"set yrange [0:]\n");
         fprintf(gnu_File,"set style fill transparent solid 0.5 noborder\n");
         fprintf(gnu_File,"set output '%s' \n",gnu_gif);
         fprintf(gnu_File,"plot '%s' using 1:2:(sqrt($3)) title 'Close contacts' with circles lc rgb 'gray', '%s' title 'MI' lt rgb 'red', '%s' title 'DI' lt rgb 'blue'\n", cmap_filename , MI_top_filename, DI_top_filename);
         fclose(MI_top_File);
         fclose(DI_top_File);
         fclose(cmap_File);
//End analyses when we have a PDB file
        }
 

//Write Output
        t =  clock() - t_init;
        fprintf(DI_File, "#Length:            \t%d\n",L);
        fprintf(DI_File, "#Effective Length:  \t%d\n",l);
        fprintf(DI_File, "#Num. Sequences:    \t%d\n",numseqs);
        fprintf(DI_File, "#Num. effective:    \t%.2f\n",Meff);
        fprintf(DI_File, "#Execution time:    \t%f\n",((float)t)/CLOCKS_PER_SEC);
        fprintf(DI_File, "#True Posititves MI:\t%d (CA)\t%5d (CB)\t%5d (MIN)\t%5d (Global)\n",tp_mi_ca,tp_mi_cb,tp_mi_min,tp_mi_sum);
        fprintf(DI_File, "#False Posititves MI:\t%d (CA)\t%5d (CB)\t%5d (MIN)\t%5d (Global)\n",fp_mi_ca,fp_mi_cb,fp_mi_min,fp_mi_sum);
        fprintf(DI_File, "#True Posititves DI:\t%d (CA)\t%5d (CB)\t%5d (MIN)\t%5d (Global)\n",tp_di_ca,tp_di_cb,tp_di_min,tp_di_sum);
        fprintf(DI_File, "#False Posititves DI:\t%d (CA)\t%5d (CB)\t%5d (MIN)\t%5d (Global)\n",fp_di_ca,fp_di_cb,fp_di_min,fp_di_sum);
        fprintf(DI_File, "\n");
	fprintf(DI_File, "#Rank     \tPosition\tPosition\tMI      \tPosition\tPosition\tDI      ");
	if (PDB_File!=NULL){
          fprintf(DI_File, "\tCA-Distance MI\t    Nearest MI\tCA-Distance DI\t    Nearest DI");
          fprintf(DI_File, "\tCB-Distance MI\t    Nearest MI\tCB-Distance DI\t    Nearest DI");
          fprintf(DI_File, "\tMINDistance MI\t    Nearest MI\tMINDistance DI\t    Nearest DI");
          fprintf(DI_File, "\tSuccess MI\tSuccess DI\n");
        }else{
          fprintf(DI_File, "\n");
        }
        for (rank=1;rank<n_index;rank++){
            index_di = (int)di_index[di_rank[rank]];
            index_mi = (int)mi_index[mi_rank[rank]];
            i_mi     = i_index[index_mi];
            j_mi     = j_index[index_mi];
            i_di     = i_index[index_di];
            j_di     = j_index[index_di];
            if (rank < max_rank_di && rank < max_rank_mi){
               fprintf(DI_File, "%10d\t%8d\t%8d\t%8.6f\t%8d\t%8d\t%8.6f",rank,xpos[pos[i_mi]],xpos[pos[j_mi]],MI_order[mi_rank[rank]],xpos[pos[i_di]],xpos[pos[j_di]],DI_order[di_rank[rank]]);
            }else if (rank < max_rank_di){
               fprintf(DI_File, "%10d\t%8s\t%8s\t%8s\t%8d\t%8d\t%8.6f",rank,"-","-","-",xpos[pos[i_di]],xpos[pos[j_di]],DI_order[di_rank[rank]]);
            }else if (rank < max_rank_mi){
               fprintf(DI_File, "%10d\t%8d\t%8d\t%8.6f\t%8s\t%8s\t%8s",rank,xpos[pos[i_mi]],xpos[pos[j_mi]],MI_order[mi_rank[rank]],"-","-","-");
            }else{
               break;
            }
	    if (PDB_File!=NULL){
             fprintf(DI_File, "\t%14.3f\t%14.3f\t%14.3f\t%14.3f",contact_ca[xpos[pos[i_mi]]][xpos[pos[j_mi]]],NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_ca),contact_ca[xpos[pos[i_di]]][xpos[pos[j_di]]],NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_ca));
             fprintf(DI_File, "\t%14.3f\t%14.3f\t%14.3f\t%14.3f",contact_cb[xpos[pos[i_mi]]][xpos[pos[j_mi]]],NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_cb),contact_cb[xpos[pos[i_di]]][xpos[pos[j_di]]],NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_cb));
             fprintf(DI_File, "\t%14.3f\t%14.3f\t%14.3f\t%14.3f",contact_min[xpos[pos[i_mi]]][xpos[pos[j_mi]]],NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_min),contact_min[xpos[pos[i_di]]][xpos[pos[j_di]]],NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_min));
	     if (  NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_ca) < (double)THRESHOLD_CA
                || NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_cb) < (double)THRESHOLD_CB
                || NearestContact(xpos[pos[i_mi]],xpos[pos[j_mi]],fragment_size,L,contact_min) < (double)THRESHOLD_MIN){
	        fprintf(DI_File, "\t%10s","yes");
	     }else{
	        fprintf(DI_File, "\t%10s","no");
	     }
	     if (  NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_ca) < (double)THRESHOLD_CA
                || NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_cb) < (double)THRESHOLD_CB
                || NearestContact(xpos[pos[i_di]],xpos[pos[j_di]],fragment_size,L,contact_min) < (double)THRESHOLD_MIN){
	        fprintf(DI_File, "\t%10s\n","yes");
	     }else{
	        fprintf(DI_File, "\t%10s\n","no");
	     }
            }else{
             fprintf(DI_File, "\n");
            }
        }
        fclose(DI_File);

//Free memory
	if (PDB_File!=NULL){
          for (i=0; i<=protein.number_of_res; i++){
	    free(cmap_ca[i]);
	    free(cmap_cb[i]);
	    free(cmap_min[i]);
	    free(cmap_sum[i]);
	  }
	  for (i=0; i<=nr;i++) free(contact_ca[i]);
	  for (i=0; i<=nr;i++) free(contact_cb[i]);
	  for (i=0; i<=nr;i++) free(contact_min[i]);
  	  free(cmap_ca);
	  free(cmap_cb);
	  free(cmap_min);
	  free(cmap_sum);
	  free(contact_ca);
	  free(contact_cb);
	  free(contact_min);
	}
	for (i=0; i<l; i++){
          for (j=0; j<l; j++){
            for (A=0; A<q; A++){ free(ff[i][j][A]);}
            free(ff[i][j]);
            }
          free(DI[i]);
          free(MI[i]);
          free(f[i]);
          free(ff[i]);
          }
        free(DI);
        free(MI);
        free(f);
        free(ff);
        free(MI_order);
        free(DI_order);
        free(mi_rank);
        free(di_rank);
        free(mi_index);
        free(di_index);
	free(i_index);
        free(j_index);
        free(pos);
        free(xpos);
        free(cmap_filename);
        free(MI_top_filename);
        free(DI_top_filename);
        free(gnu_filename);
        free(gnu_gif);

//Exit
        printf("Done\n");
        exit(0);

 }


