#include "raDI.h"

int RankOrder(minimum_gap,fragment_size,l,n_index,pos,MI,SSA,mi_rank,i_index,j_index,mi_index,MI_order)
int minimum_gap,fragment_size,n_index,l;
int *mi_rank,*pos,*i_index,*j_index;
double *MI_order,*mi_index;
double **MI;
secondary_structure *SSA;
{

	int i, j, k,index,i_mi,j_mi,i_di,j_di,index_mi,rank,ss_i,ss_j;
        int max_rank,keep,pos_i,pos_j;
	int *mi_pairs;
        void sort2();


//      Request memory
        mi_pairs  = (int*) calloc(n_index+1,sizeof(int));

//      Rank Order MI (small@1 -> largest@(index-1))
        index = 1;
        for (i=0; i<l; i++){
        for (j=i+1; j<l; j++){
            ss_i  = SSA[pos[i]].order;
            ss_j  = SSA[pos[j]].order;
            if (ss_i != ss_j && abs(pos[i]-pos[j])>minimum_gap){ 
               mi_index[index]  = (double) index;
               i_index[index]   = i;
               j_index[index]   = j;
               MI_order[index]  = MI[i][j];
               index++;
               }
        }} 
        sort2(index-1,MI_order,mi_index);
        for (i=0; i<n_index;i++){mi_pairs[i]=0;}
        mi_pairs[n_index-1]=1;
        for (i=n_index-2;i>0;i--){
          keep=1;
          index_mi   = (int)mi_index[i];
          i_mi       = i_index[index_mi];
          j_mi       = j_index[index_mi];
          pos_i      = pos[i_mi];
          pos_j      = pos[j_mi];
          for (j=n_index-1;j>i;j--){
            index_mi   = (int)mi_index[j];
            i_mi       = i_index[index_mi];
            j_mi       = j_index[index_mi];
            if (mi_pairs[j] && 
                abs(pos_i -pos[i_mi])< fragment_size/2 && 
                abs(pos_j -pos[j_mi])< fragment_size/2)  { keep=0; break;}
            }
          if(keep){mi_pairs[i]=1;}
        }
        for (i=n_index-1;i>0;i--){mi_rank[n_index-i]=i;}
        rank=1;
        for (i=n_index-1;i>0;i--){
          if (mi_pairs[i]){
             mi_rank[rank]=i;
             rank++;
             }
        }
        max_rank=rank;

//Free memory
	free(mi_pairs);

  return max_rank;

 }


