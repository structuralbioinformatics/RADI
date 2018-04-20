#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"svdcmp.h"

#define NR_END 1
#define FREE_ARG char*

/* This routine computes the SVD of a square matrix */
/* A is the matrix, of dimension m x n, A = UWV^{t}, where U substitutes A in the output, V (not its transpose) is given as output, same as the diagonal matrix W, given as a vector */

static double at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static double maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void svdcmp(a,m,n,w,v)
double **a,*w,**v;
int m,n;
{
	int flag,i,its,j,jj,k,l,nm,Maxit;
	double c,f,h,s,x,y,z;
	double anorm=0.0,g=0.0,scale=0.0,num_zero=1e-35;
	double *rv1,*vector();
	void nrerror(),free_vector();
	Maxit = 1000;
	if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
	rv1=vector(1,n);
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k-1][i-1]);
			if (scale>num_zero) { /* Check in your email for NaN */
				for (k=i;k<=m;k++) {
					a[k-1][i-1] /= scale;
					s += a[k-1][i-1]*a[k-1][i-1];
				}
				f=a[i-1][i-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i-1][i-1]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k-1][i-1]*a[k-1][j-1];
						f=s/h;
						for (k=i;k<=m;k++) a[k-1][j-1] += f*a[k-1][i-1];
					}
				}
				for (k=i;k<=m;k++) a[k-1][i-1] *= scale;
			}
		}
		w[i-1]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i-1][k-1]);
			if (scale>num_zero) { /* Check in your email for NaN */
				for (k=l;k<=n;k++) {
					a[i-1][k-1] /= scale;
					s += a[i-1][k-1]*a[i-1][k-1];
				}
				f=a[i-1][l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i-1][l-1]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i-1][k-1]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j-1][k-1]*a[i-1][k-1];
						for (k=l;k<=n;k++) a[j-1][k-1] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i-1][k-1] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i-1])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j-1][i-1]=(a[i-1][j-1]/a[i-1][l-1])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i-1][k-1]*v[k-1][j-1];
					for (k=l;k<=n;k++) v[k-1][j-1] += s*v[k-1][i-1];
				}
			}
			for (j=l;j<=n;j++) v[i-1][j-1]=v[j-1][i-1]=0.0;
		}
		v[i-1][i-1]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i-1];
		if (i < n)
			for (j=l;j<=n;j++) a[i-1][j-1]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k-1][i-1]*a[k-1][j-1];
					f=(s/a[i-1][i-1])*g;
					for (k=i;k<=m;k++) a[k-1][j-1] += f*a[k-1][i-1];
				}
			}
			for (j=i;j<=m;j++) a[j-1][i-1] *= g;
		} else {
			for (j=i;j<=m;j++) a[j-1][i-1]=0.0;
		}
		++a[i-1][i-1];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=Maxit;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm-1])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i-1];
					h=PYTHAG(f,g);
					w[i-1]=h;
					h=1.0/h;
					c=g*h;
					s=(-f*h);
					for (j=1;j<=m;j++) {
						y=a[j-1][nm-1];
						z=a[j-1][i-1];
						a[j-1][nm-1]=y*c+z*s;
						a[j-1][i-1]=z*c-y*s;
					}
				}
			}
			z=w[k-1];
			if (l == k) {
				if (z < 0.0) {
					w[k-1] = -z;
					for (j=1;j<=n;j++) v[j-1][k-1]=(-v[j-1][k-1]);
				}
				break;
			}
			if (its == Maxit) nrerror("No convergence in %d SVDCMP iterations", Maxit);
			x=w[l-1];
			nm=k-1;
			y=w[nm-1];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i-1];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj-1][j-1];
					z=v[jj-1][i-1];
					v[jj-1][j-1]=x*c+z*s;
					v[jj-1][i-1]=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j-1]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj-1][j-1];
					z=a[jj-1][i-1];
					a[jj-1][j-1]=y*c+z*s;
					a[jj-1][i-1]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k-1]=x;
		}
	}
	free_vector(rv1,1,n);
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	printf("Numerical Recipes run-time error...\n");
	printf("%s\n",error_text);
	printf("...now exiting to system...\n");
	exit(1);
}

double *vector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


#undef SIGN
#undef MAX
#undef PYTHAG
