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
  for (i=n_minim;i<=n_maxim;i++){
  for (j=m_minim;j<=m_maxim;j++){
   if (start) {x=d[i][j];start=0;}
   x=fmin(x,d[i][j]);
  }}

  return x;
}
