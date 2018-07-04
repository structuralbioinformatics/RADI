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

void Results_MI_DI(PDB_File, DI_filename, cmap_filename, cmap_off_filename, MI_top_filename, DI_top_filename,
                   gnu_filename, gnu_gif, cmap_ca, contact_ca, cmap_cb, contact_cb, cmap_min, contact_min, nr,
                   l, L, minimum_gap, numseqs, fragment_size, Meff, pos, xpos, SSA, MI, DI, time, verbose)
  FILE *PDB_File;
  secondary_structure *SSA;
  int nr, l, L, numseqs, verbose, minimum_gap, fragment_size, *xpos, *pos;
  int **cmap_ca, **cmap_cb, **cmap_min;
  double **DI, **MI, **contact_ca, **contact_cb, **contact_min, *time, Meff;
  char *gnu_gif, *cmap_filename, *cmap_off_filename, *MI_top_filename, *DI_top_filename, *DI_filename, *gnu_filename;
{

  FILE *DI_File, *DI_top_File, *MI_top_File, *cmap_File, *cmap_offFile, *gnu_File;
  int i, j, k, n_index, ss_i, ss_j, rank, index_mi, index_di, swap, print_return;
  int max_rank_mi, max_rank_di, i_mi, j_mi, i_di, j_di, nt_core, ct_core, toprank;
  int *di_rank, *mi_rank, *j_index, *i_index;
  int tp_mi_ca, tp_mi_cb, tp_mi_min, tp_mi_sum, tp_di_ca, tp_di_cb, tp_di_min, tp_di_sum;
  int fp_mi_ca, fp_mi_cb, fp_mi_min, fp_mi_sum, fp_di_ca, fp_di_cb, fp_di_min, fp_di_sum;
  double *mi_index, *di_index, *DI_order, *MI_order;
  int RankOrder();
  double NearestContact();
  void sort2();

  //Open plot/output names
  DI_File     = NULL;
  DI_File     = fopen(DI_filename, "w");
  if (DI_File == NULL){
    printf("--Problems opening the OUTPUT file %s\n",DI_filename);
    exit(2);
  }
  MI_top_File = NULL;
  DI_top_File = NULL;
  cmap_File   = NULL;
  cmap_offFile= NULL;
  gnu_File    = NULL;
  DI_top_File = fopen(DI_top_filename, "w");
  MI_top_File = fopen(MI_top_filename, "w");
  cmap_File   = fopen(cmap_filename, "w");
  cmap_offFile= fopen(cmap_off_filename, "w");
  gnu_File    = fopen(gnu_filename, "w");
  if (DI_top_File == NULL){
    printf("--Problems opening the OUTPUT file %s\n", DI_top_filename);
    exit(2);
  }
  if (MI_top_File == NULL){
    printf("--Problems opening the OUTPUT file %s\n", MI_top_filename);
    exit(2);
  }
  if (cmap_File == NULL){
    printf("--Problems opening the OUTPUT file %s\n", cmap_filename);
    exit(2);
  }
  if (cmap_offFile == NULL){
    printf("--Problems opening the OUTPUT file %s\n", cmap_off_filename);
    exit(2);
  }
  if (gnu_File == NULL){
    printf("--Problems opening the OUTPUT file %s\n", gnu_filename);
    exit(2);
  }

  //Reorder by MI and DI
  //      Calculate size
  if (verbose) printf("\t-- Parameters: minimum distance in sequence %d\n", minimum_gap);
  n_index = 1;
  for (i=0; i<l; i++){
    for (j=i+1; j<l; j++){
      ss_i  = SSA[pos[i]].order;
      ss_j  = SSA[pos[j]].order;
      if (ss_i != ss_j  && abs(pos[i]-pos[j])>minimum_gap ){ n_index++; }
    }
  }

  //      Request memory
  mi_rank   = (int*) calloc(n_index+1, sizeof(int));
  di_rank   = (int*) calloc(n_index+1, sizeof(int));
  DI_order  = (double*) calloc(n_index+1, sizeof(double));
  MI_order  = (double*) calloc(n_index+1, sizeof(double));
  mi_index  = (double*) calloc(n_index+1, sizeof(double));
  di_index  = (double*) calloc(n_index+1, sizeof(double));
  j_index   = (int*) calloc(n_index+1, sizeof(int));
  i_index   = (int*) calloc(n_index+1, sizeof(int));

  //      Rank Order MI (small@1 -> largest@(index-1))
  if (verbose) printf("\t-- Rank top MI scores\n");
  max_rank_mi = RankOrder(minimum_gap, fragment_size, l, n_index, pos, MI, SSA, mi_rank, i_index, j_index, mi_index, MI_order);
  if (verbose) printf("\t-- Maximum number of MI %d\n", max_rank_mi);

  //      Rank Order DI (small@1 -> largest@(index-1))
  if (verbose) printf("\t-- Rank top DI scores\n");
  max_rank_di=RankOrder(minimum_gap, fragment_size, l, n_index, pos, DI, SSA, di_rank, i_index, j_index, di_index, DI_order);
  if (verbose) printf("\t-- Maximum number of DI %d\n", max_rank_di);

  //Analyses
  if (fmin(max_rank_di, max_rank_mi)<TOPRANK){ toprank=fmin(max_rank_di, max_rank_mi)/2; }else{ toprank=TOPRANK; }
  if (verbose) printf("\t-- TOPRANK  %d\n",toprank);
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
  for (rank=1; rank<=toprank; rank++){
    index_di = (int)di_index[di_rank[rank]];
    index_mi = (int)mi_index[mi_rank[rank]];
    i_mi     = i_index[index_mi];
    j_mi     = j_index[index_mi];
    i_di     = i_index[index_di];
    j_di     = j_index[index_di];
    if (NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_ca) < (double)THRESHOLD_CA){ tp_mi_ca++; }else{ fp_mi_ca++; }
    if (NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_ca) < (double)THRESHOLD_CA){ tp_di_ca++; }else{ fp_di_ca++; }
    if (NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_cb) < (double)THRESHOLD_CB){ tp_mi_cb++; }else{ fp_mi_cb++; }
    if (NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_cb) < (double)THRESHOLD_CB){ tp_di_cb++; }else{ fp_di_cb++; }
    if (NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_min) < (double)THRESHOLD_MIN){ tp_mi_min++; }else{ fp_mi_min++; }
    if (NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_min) < (double)THRESHOLD_MIN){ tp_di_min++; }else{ fp_di_min++; }
    if (NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_ca) < (double)THRESHOLD_CA
      || NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size,L ,contact_cb) < (double)THRESHOLD_CB
      || NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size,L ,contact_min) < (double)THRESHOLD_MIN){ tp_mi_sum++; }else{ fp_mi_sum++; }
    if (NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size,L ,contact_ca) < (double)THRESHOLD_CA
      || NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size,L ,contact_cb) < (double)THRESHOLD_CB
      || NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size,L ,contact_min) < (double)THRESHOLD_MIN){ tp_di_sum++; }else{ fp_di_sum++; }
  }
  if (verbose) printf("\t-- Successes on MI: %d\n", tp_mi_sum);
  if (verbose) printf("\t-- Successes on DI: %d\n", tp_di_sum);

  //Plot
  fprintf(cmap_File, "#pos1\tpos2\tradius\n");
  fprintf(cmap_offFile, "#pos1\tpos2\tradius\n");
  fprintf(MI_top_File, "#pos1\tpos2\n");
  fprintf(DI_top_File, "#pos1\tpos2\n");
  nt_core=L;
  ct_core=1;
  for (k=0; k<l; k++){ if (xpos[pos[k]] > ct_core) ct_core=xpos[pos[k]];}
  for (k=0; k<l; k++){ if (xpos[pos[k]] < nt_core) nt_core=xpos[pos[k]];}
  for (i=1; i<=nr; i++){
    for (j=i+1; j<=nr; j++){
      if  (cmap_ca[i][j] == 1 || cmap_cb[i][j] == 1  || cmap_min[i][j] == 1) {
        fprintf(cmap_File,"%3d\t%3d\t%3d\n",i,j,2);
        fprintf(cmap_File,"%3d\t%3d\t%3d\n",j,i,2);
        if (i<nt_core || j<nt_core){
          fprintf(cmap_offFile,"%3d\t%3d\t%3d\n",i,j,2);
          fprintf(cmap_offFile,"%3d\t%3d\t%3d\n",j,i,2);
        }
        if (i>ct_core || j>ct_core){
          fprintf(cmap_offFile,"%3d\t%3d\t%3d\n",i,j,2);
          fprintf(cmap_offFile,"%3d\t%3d\t%3d\n",j,i,2);
        }
      }
    }
  }
  fflush(cmap_File);
  fflush(cmap_offFile);
  for (rank=1; rank<=toprank; rank++){
    index_mi = (int)mi_index[mi_rank[rank]];
    i_mi     = i_index[index_mi];
    j_mi     = j_index[index_mi];
    if (xpos[pos[i_mi]] < xpos[pos[j_mi]]){
      fprintf(MI_top_File, "%3d\t%3d\n", xpos[pos[i_mi]], xpos[pos[j_mi]]);
    }else{
      fprintf(MI_top_File, "%3d\t%3d\n", xpos[pos[j_mi]], xpos[pos[i_mi]]);
    }
    index_di = (int)di_index[di_rank[rank]];
    i_di     = i_index[index_di];
    j_di     = j_index[index_di];
    if (xpos[pos[i_mi]] > xpos[pos[j_mi]]){
      fprintf(DI_top_File, "%3d\t%3d\n", xpos[pos[i_di]], xpos[pos[j_di]]);
    }else{
      fprintf(DI_top_File, "%3d\t%3d\n", xpos[pos[j_di]], xpos[pos[i_di]]);
    }
  }
  fflush(MI_top_File);
  fflush(DI_top_File);
  fprintf(gnu_File,"set term gif\n");
  fprintf(gnu_File,"set xrange [0:]\n");
  fprintf(gnu_File,"set yrange [0:]\n");
  fprintf(gnu_File,"set style fill transparent solid 0.5 noborder\n");
  fprintf(gnu_File,"set output '%s' \n", gnu_gif);
  fprintf(gnu_File,"plot '%s' using 1:2:(sqrt($3)) title 'Close contacts' with circles lc rgb 'gray', '%s' title 'MI' lt rgb 'red', '%s' title 'DI' lt rgb 'blue', '%s' using 1:2:(sqrt($3)) title 'Off-side contacts' with circles lc rgb '#FFC0CB'\n", cmap_filename , MI_top_filename, DI_top_filename, cmap_off_filename);


  //Write Output
  fprintf(DI_File, "#Length:            \t%d\n", L);
  fprintf(DI_File, "#Effective Length:  \t%d\n", l);
  fprintf(DI_File, "#Coverage:          \t%.2f\n", (double)l/(double)L);
  fprintf(DI_File, "#Num. Sequences:    \t%d\n", numseqs);
  fprintf(DI_File, "#Num. effective:    \t%.2f\n", Meff);
  fprintf(DI_File, "#Execution time MI: \t%f\n", time[0]);
  fprintf(DI_File, "#Execution time DI: \t%f\n", time[1]);
  fprintf(DI_File, "#True Posititves MI:\t%d (CA)\t%5d (CB)\t%5d (MIN)\t%5d (Global)\n",tp_mi_ca, tp_mi_cb, tp_mi_min, tp_mi_sum);
  fprintf(DI_File, "#False Posititves MI:\t%d (CA)\t%5d (CB)\t%5d (MIN)\t%5d (Global)\n",fp_mi_ca, fp_mi_cb, fp_mi_min, fp_mi_sum);
  fprintf(DI_File, "#True Posititves DI:\t%d (CA)\t%5d (CB)\t%5d (MIN)\t%5d (Global)\n",tp_di_ca, tp_di_cb, tp_di_min, tp_di_sum);
  fprintf(DI_File, "#False Posititves DI:\t%d (CA)\t%5d (CB)\t%5d (MIN)\t%5d (Global)\n",fp_di_ca, fp_di_cb, fp_di_min, fp_di_sum);
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
      fprintf(DI_File, "%10d\t%8d\t%8d\t%8.6f\t%8d\t%8d\t%8.6f",rank, xpos[pos[i_mi]], xpos[pos[j_mi]], MI_order[mi_rank[rank]], xpos[pos[i_di]], xpos[pos[j_di]], DI_order[di_rank[rank]]);
    }else if (rank < max_rank_di){
      fprintf(DI_File, "%10d\t%8s\t%8s\t%8s\t%8d\t%8d\t%8.6f",rank, "-", "-", "-", xpos[pos[i_di]], xpos[pos[j_di]], DI_order[di_rank[rank]]);
    }else if (rank < max_rank_mi){
      fprintf(DI_File, "%10d\t%8d\t%8d\t%8.6f\t%8s\t%8s\t%8s", rank, xpos[pos[i_mi]], xpos[pos[j_mi]], MI_order[mi_rank[rank]], "-", "-", "-");
    }else{
       break;
    }
    if (PDB_File!=NULL) {
     if (rank < max_rank_di && rank < max_rank_mi){
      fprintf(DI_File, "\t%14.3f\t%14.3f\t%14.3f\t%14.3f", contact_ca[xpos[pos[i_mi]]][xpos[pos[j_mi]]],
              NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_ca), contact_ca[xpos[pos[i_di]]][xpos[pos[j_di]]],
              NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_ca));
     }else if (rank < max_rank_di){
      fprintf(DI_File, "\t%14s\t%14s\t%14.3f\t%14.3f", "-",
              "-", contact_ca[xpos[pos[i_di]]][xpos[pos[j_di]]],
              NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_ca));
     }else if (rank < max_rank_mi){
      fprintf(DI_File, "\t%14.3f\t%14.3f\t%14s\t%14s", contact_ca[xpos[pos[i_mi]]][xpos[pos[j_mi]]],
              NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_ca), "-",
              "-");
     }else{
      fprintf(DI_File, "\t%14s\t%14s\t%14s\t%14s", "-","-","-","-");
     }
     if (rank < max_rank_di && rank < max_rank_mi){
      fprintf(DI_File, "\t%14.3f\t%14.3f\t%14.3f\t%14.3f", contact_cb[xpos[pos[i_mi]]][xpos[pos[j_mi]]],
              NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_cb), contact_cb[xpos[pos[i_di]]][xpos[pos[j_di]]],
              NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_cb));
     }else if (rank < max_rank_di){
      fprintf(DI_File, "\t%14s\t%14s\t%14.3f\t%14.3f", "-",
              "-", contact_cb[xpos[pos[i_di]]][xpos[pos[j_di]]],
              NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_cb));
     }else if (rank < max_rank_mi){
      fprintf(DI_File, "\t%14.3f\t%14.3f\t%14s\t%14s", contact_cb[xpos[pos[i_mi]]][xpos[pos[j_mi]]],
              NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_cb), "-",
              "-");
     }else{
      fprintf(DI_File, "\t%14s\t%14s\t%14s\t%14s", "-","-","-","-");
     }
     if (rank < max_rank_di && rank < max_rank_mi){
      fprintf(DI_File, "\t%14.3f\t%14.3f\t%14.3f\t%14.3f", contact_min[xpos[pos[i_mi]]][xpos[pos[j_mi]]],
              NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_min), contact_min[xpos[pos[i_di]]][xpos[pos[j_di]]],
              NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_min));
     }else if (rank < max_rank_di){
      fprintf(DI_File, "\t%14s\t%14s\t%14.3f\t%14.3f", "-",
              "-", contact_min[xpos[pos[i_di]]][xpos[pos[j_di]]],
              NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_min));
     }else if (rank < max_rank_mi){
      fprintf(DI_File, "\t%14.3f\t%14.3f\t%14s\t%14s", contact_min[xpos[pos[i_mi]]][xpos[pos[j_mi]]],
              NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_min), "-",
              "-");
     }else{
      fprintf(DI_File, "\t%14s\t%14s\t%14s\t%14s", "-","-","-","-");
     }
     if (NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_ca) < (double)THRESHOLD_CA
        || NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_cb) < (double)THRESHOLD_CB
        || NearestContact(xpos[pos[i_mi]], xpos[pos[j_mi]], fragment_size, L, contact_min) < (double)THRESHOLD_MIN){
        if (rank < max_rank_mi){fprintf(DI_File, "\t%10s","yes");}else{fprintf(DI_File, "\t%10s","-");}
     }else{
        if (rank < max_rank_mi){fprintf(DI_File, "\t%10s","no");}else{fprintf(DI_File, "\t%10s","-");}
     }
     if (NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_ca) < (double)THRESHOLD_CA
        || NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_cb) < (double)THRESHOLD_CB
        || NearestContact(xpos[pos[i_di]], xpos[pos[j_di]], fragment_size, L, contact_min) < (double)THRESHOLD_MIN){
         if (rank < max_rank_di){fprintf(DI_File, "\t%10s\n", "yes");}else{fprintf(DI_File, "\t%10s","-");}
     }else{
         if (rank < max_rank_di){fprintf(DI_File, "\t%10s\n", "no");}else{fprintf(DI_File, "\t%10s","-");}
     }
    }else{
      fprintf(DI_File, "\n");
    }
  }

  //Close files
  fclose(DI_File);
  fclose(cmap_File);
  fclose(cmap_offFile);
  fclose(gnu_File);
  fclose(MI_top_File);
  fclose(DI_top_File);

  //Free memory
  free(MI_order);
  free(DI_order);
  free(mi_rank);
  free(di_rank);
  free(mi_index);
  free(di_index);
  free(i_index);
  free(j_index);

}
