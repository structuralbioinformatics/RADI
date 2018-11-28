#include "raDI.h"


double NearestContact(n,m,length,L,d)
int n,m,length,L;
double **d;
{
  int i,j,start,n_minim,n_maxim,m_minim,m_maxim;
  double x;

  n_minim = (int)fmax(1,n-length/2);
  n_maxim = (int)fmin(L,n+length/2);
  m_minim = (int)fmax(1,m-length/2);
  m_maxim = (int)fmin(L,m+length/2);

  start=1;
  x=100.0;
  if (n==0 || m==0) {return x;}
  for (i=n_minim;i<=n_maxim;i++){
  for (j=m_minim;j<=m_maxim;j++){
   if (start && d[i][j]>0) {x=d[i][j];start=0;}
   if (d[i][j]>0) {x=fmin(x,d[i][j]);}
  }}

  return x;
}


