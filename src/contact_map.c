#include "raDI.h"


int **contact_map(protein, gap, type, cutoff, d, verbose)
  PROT protein;
  int verbose,gap;
  char *type;
  double cutoff,**d;
{
  int   i, j, ic, n, m, skp, nr;
  int   ia[MAXATOM], ja[MAXATOM], nc[MAXCHAIN];
  int   **cmap;
  float x, y, z, dx, dy, dz;

  nr = protein.number_of_res;
  if (verbose) printf("\t-- Check contacts %d X %d between %s-atom-distance < %f A\n", nr, nr, type, cutoff);
  if (verbose) printf("\t-- For a max. possible total of Atoms=%d\n", protein.number_of_atoms);
  cmap = (int**)malloc((nr + 1) * sizeof(int*));
  for (i=0; i<=nr; i++){ cmap[i]=(int*) calloc((nr+1), sizeof(int)); }
  for (i=0; i<=nr; i++){ cmap[i][j]=0; }

  x  = 0.0;
  y  = 0.0;
  z  = 0.0;
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;
  ic = 0;

  for (i=0; i<protein.number_of_atoms; i++){
    if (verbose) printf("%d ", protein.atoms_of_prot[i].res_number);
    if(!strcmp(protein.atoms_of_prot[i].name, "CA  ") && (!strcmp(type, "ca") || !strcmp(type, "CA")) ){
      if (verbose) printf("\t-- It's CA\n");
      x=protein.atoms_of_prot[i].x;
      y=protein.atoms_of_prot[i].y;
      z=protein.atoms_of_prot[i].z;
    }
    if ((!strcmp(protein.atoms_of_prot[i].name, "CB  ") || (!strcmp(protein.atoms_of_prot[i].name, "CA  ") && !strcmp(protein.atoms_of_prot[i].res, "GLY")) )
     && (!strcmp(type,"cb") || !strcmp(type,"CB")) ){
      if (verbose) printf("\t-- It's CB\n");
      x=protein.atoms_of_prot[i].x;
      y=protein.atoms_of_prot[i].y;
      z=protein.atoms_of_prot[i].z;
    }
    if (((strcmp(protein.atoms_of_prot[i].name, "C   ") && strcmp(protein.atoms_of_prot[i].name, "N   ") && strcmp(protein.atoms_of_prot[i].name, "O   ")
       && strcmp(protein.atoms_of_prot[i].name,"CA  ")) || (!strcmp(protein.atoms_of_prot[i].res,"GLY")) ) && (!strcmp(type,"min")||!strcmp(type,"MIN")) ){
      if (verbose) printf("\t-- It's BB\n");
      x=protein.atoms_of_prot[i].x;
      y=protein.atoms_of_prot[i].y;
      z=protein.atoms_of_prot[i].z;
    }
    if (i+1<protein.number_of_atoms) {
      if (strcmp(protein.atoms_of_prot[i].chain,protein.atoms_of_prot[i+1].chain)) {
        nc[ic]=protein.atoms_of_prot[i].res_number;
        ic++;
        if (ic==MAXCHAIN){
          printf("Too many chains: %s vs %s in residue %d\n", protein.atoms_of_prot[i].chain, protein.atoms_of_prot[i+1].chain, protein.atoms_of_prot[i].res_number);
          exit(0);
        }
      }
    }
    for (j=i+1; j<protein.number_of_atoms; j++){
      n=protein.atoms_of_prot[i].res_number;
      m=protein.atoms_of_prot[j].res_number;
      if (abs(protein.atoms_of_prot[i].res_number - protein.atoms_of_prot[j].res_number) > gap || strcmp(protein.atoms_of_prot[i].chain, protein.atoms_of_prot[j].chain)){
        if (!strcmp(protein.atoms_of_prot[j].name,"CA  ") && (!strcmp(type,"ca") || !strcmp(type,"CA")) && !strcmp(protein.atoms_of_prot[i].name,"CA  ")){
          dx=protein.atoms_of_prot[j].x - x ;
          dy=protein.atoms_of_prot[j].y - y;
          dz=protein.atoms_of_prot[j].z - z;
          d[protein.atoms_of_prot[i].res_number][protein.atoms_of_prot[j].res_number]=(double)sqrt(dx*dx + dy*dy + dz*dz);
          d[protein.atoms_of_prot[j].res_number][protein.atoms_of_prot[i].res_number]=(double)sqrt(dx*dx + dy*dy + dz*dz);
          ia[protein.atoms_of_prot[i].res_number]=i;
          ja[protein.atoms_of_prot[j].res_number]=j;
        }
        if (   (!strcmp(protein.atoms_of_prot[j].name,"CB  ") || (!strcmp(protein.atoms_of_prot[j].name,"CA  ") && !strcmp(protein.atoms_of_prot[j].res,"GLY")) ) &&
               (!strcmp(protein.atoms_of_prot[i].name,"CB  ") || (!strcmp(protein.atoms_of_prot[i].name,"CA  ") && !strcmp(protein.atoms_of_prot[i].res,"GLY")) )
            && (!strcmp(type,"cb")||!strcmp(type,"CB")) ){
          dx=protein.atoms_of_prot[j].x - x ;
          dy=protein.atoms_of_prot[j].y - y;
          dz=protein.atoms_of_prot[j].z - z;
          d[protein.atoms_of_prot[i].res_number][protein.atoms_of_prot[j].res_number]=(double)sqrt(dx*dx + dy*dy + dz*dz);
          d[protein.atoms_of_prot[j].res_number][protein.atoms_of_prot[i].res_number]=(double)sqrt(dx*dx + dy*dy + dz*dz);
          ia[protein.atoms_of_prot[i].res_number]=i;
          ja[protein.atoms_of_prot[j].res_number]=j;
        }
        if ((!strcmp(type,"MIN")||!strcmp(type,"min")) ){
          dx=protein.atoms_of_prot[j].x - x ;
          dy=protein.atoms_of_prot[j].y - y;
          dz=protein.atoms_of_prot[j].z - z;
          if ( d[protein.atoms_of_prot[i].res_number][protein.atoms_of_prot[j].res_number] <= 0) {
            d[protein.atoms_of_prot[i].res_number][protein.atoms_of_prot[j].res_number] = (double)sqrt(dx*dx + dy*dy + dz*dz);
            d[protein.atoms_of_prot[j].res_number][protein.atoms_of_prot[i].res_number] = (double)sqrt(dx*dx + dy*dy + dz*dz);
            ia[protein.atoms_of_prot[i].res_number]=i;
            ja[protein.atoms_of_prot[j].res_number]=j;
          }
          if ((double)sqrt(dx*dx + dy*dy +  dz*dz) < d[protein.atoms_of_prot[i].res_number][protein.atoms_of_prot[j].res_number]){
            d[protein.atoms_of_prot[i].res_number][protein.atoms_of_prot[j].res_number]=(double)sqrt(dx*dx + dy*dy + dz*dz);
            d[protein.atoms_of_prot[j].res_number][protein.atoms_of_prot[i].res_number]=(double)sqrt(dx*dx + dy*dy + dz*dz);
            ia[protein.atoms_of_prot[i].res_number]=i;
            ja[protein.atoms_of_prot[j].res_number]=j;
          }
        }
      }
    }
  }

  for (n=1; n<=protein.number_of_res; n++){
    skp=0;
    for (i=0; i<ic; i++){
      if (n==nc[i]){
        skp=1;
        break;
      }
    }
    for (m=n+1; m<=protein.number_of_res; m++){
      if (m>=n+gap || skp==1){
        if (d[n][m]<=cutoff && d[n][m]>0.0){
          if (verbose) printf("\t-- CONTACT (type %s) %5d\t%5d\t%10.5e\t%c\t%c\t%s\t%s\n", type, n, m, d[n][m], protein.atoms_of_prot[ia[n]].res_code[0], protein.atoms_of_prot[ja[m]].res_code[0], protein.atoms_of_prot[ia[n]].name, protein.atoms_of_prot[ja[m]].name);
          cmap[n][m]=1;
          cmap[m][n]=1;
        }
      }
    }
  }

  //report
  return cmap;

}
