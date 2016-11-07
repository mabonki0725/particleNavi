#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : - fabs(a))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1 = (a),maxarg2 = (b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1 = (a),iminarg2 = (b),(iminarg1 < (iminarg2) ? (iminarg1) : iminarg2))

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

int svdcomp(double **a, int nRows, int nCols, double *w, double **v);
/**********************
  SVDŒŸØ
  by m.nakai
**********************/
#if 0
int main(argc,argv)
int argc;
char *argv[];
{
    int i;
#if 1
    double **a,**aw;
    double s[2];
    double **v,**vT;
    double **ss;
    a = (double **)comMalloc(sizeof(double *)*4);
    aw= (double **)comMalloc(sizeof(double *)*4);

    v = (double **)comMalloc(sizeof(double *)*2);
    ss =(double **)comMalloc(sizeof(double *)*2);
    vT= (double **)comMalloc(sizeof(double *)*2);

    for(i=0;i<4;i++) {
      a[i]=(double *)comMalloc(sizeof(double)*2);
     aw[i]=(double *)comMalloc(sizeof(double)*2);
    }
    for(i=0;i<2;i++) {
      v[i]=(double *)comMalloc(sizeof(double)*2);
     vT[i]=(double *)comMalloc(sizeof(double)*2);

     ss[i]=(double *)comMalloc(sizeof(double)*2);
    }
    for(i=0;i<4;i++) {
      a[i][0] = 1.0 + i*2.0;
      a[i][1] = 2.0 + i*2.0;
    }

    /* “ÁˆÙ’l•ª‰ð  a->U s->diag(ƒ°) v->V*/
    svdcomp(a,4,2,s,v);
    ss[0][0]=s[0];
    ss[1][1]=s[1];

    //A=U*ss*V' -> Œ³‚Ì’l‚É‚È‚é
    mxTrns(v,2,2,vT);
    mxMult(a,4,2,2,ss,aw);
    mxMult(aw,4,2,2,vT,a);
#else
    double **a;
    double s[3];
    double **v;
    a = (double **)comMalloc(sizeof(double *)*5);
    v = (double **)comMalloc(sizeof(double *)*3);
    for(i=0;i<5;i++) {
      a[i]=(double *)comMalloc(sizeof(double)*3);
    }
    for(i=0;i<3;i++) {
      v[i]=(double *)comMalloc(sizeof(double)*3);
    }

    for(i=1;i<5;i++) {
      a[i][1] = 1.0 + (i-1)*2.0;
      a[i][2] = 2.0 + (i-1)*2.0;
    }

    svdcmp(a,4,2,s,v);
#endif
    return(0);

}
#endif
// prints an arbitrary size matrix to the standard output
void printMatrix(double **a, int rows, int cols);
void printMatrix(double **a, int rows, int cols) {
  int i,j;

  for(i=0;i<rows;i++) {
    for(j=0;j<cols;j++) {
      printf("%.4f ",a[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

// prints an arbitrary size vector to the standard output
void printVector(double *v, int size);
void printVector(double *v, int size) {
  int i;

  for(i=0;i<size;i++) {
    printf("%.4f ",v[i]);
  }
  printf("\n\n");
}

// calculates sqrt( a^2 + b^2 ) with decent precision
double pythag(double a, double b);
double pythag(double a, double b) {
  double absa,absb;

  absa = fabs(a);
  absb = fabs(b);

  if(absa > absb)
    return(absa * sqrt(1.0 + SQR(absb/absa)));
  else
    return(absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

/*
  Modified from Numerical Recipes in C
  Given a matrix a[nRows][nCols], svdcmp() computes its singular value 
  decomposition, A = U * W * Vt.  A is replaced by U when svdcmp 
  returns.  The diagonal matrix W is output as a vector w[nCols].
  V (not V transpose) is output as the matrix V[nCols][nCols].
*/
int svdcomp(double **a, int nRows, int nCols, double *w, double **v);
int svdcomp(double **a, int nRows, int nCols, double *w, double **v) {
  int flag,i,its,j,jj,k,l=0,nm=0;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

  rv1 = malloc(sizeof(double)*nCols);
  if(rv1 == NULL) {
  	printf("svdcmp(): Unable to allocate vector\n");
  	return(-1);
  }

  g = scale = anorm = 0.0;
  for(i=0;i<nCols;i++) {
    l = i+1;
    rv1[i] = scale*g;
    g = s = scale = 0.0;
    if(i < nRows) {
      for(k=i;k<nRows;k++) scale += fabs(a[k][i]);
      if(scale) {
	for(k=i;k<nRows;k++) {
	  a[k][i] /= scale;
	  s += a[k][i] * a[k][i];
	}
	f = a[i][i];
	g = -SIGN(sqrt(s),f);
	h = f * g - s;
	a[i][i] = f - g;
	for(j=l;j<nCols;j++) {
	  for(s=0.0,k=i;k<nRows;k++) s += a[k][i] * a[k][j];
	  f = s / h;
	  for(k=i;k<nRows;k++) a[k][j] += f * a[k][i];
	}
	for(k=i;k<nRows;k++) a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if(i < nRows && i != nCols-1) {
      for(k=l;k<nCols;k++) scale += fabs(a[i][k]);
      if(scale)  {
	for(k=l;k<nCols;k++) {
	  a[i][k] /= scale;
	  s += a[i][k] * a[i][k];
	}
	f = a[i][l];
	g = - SIGN(sqrt(s),f);
	h = f * g - s;
	a[i][l] = f - g;
	for(k=l;k<nCols;k++) rv1[k] = a[i][k] / h;
	for(j=l;j<nRows;j++) {
	  for(s=0.0,k=l;k<nCols;k++) s += a[j][k] * a[i][k];
	  for(k=l;k<nCols;k++) a[j][k] += s * rv1[k];
	}
	for(k=l;k<nCols;k++) a[i][k] *= scale;
      }
    }
    anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));

    fflush(stdout);
  }

  for(i=nCols-1;i>=0;i--) {
    if(i < nCols-1) {
      if(g) {
	for(j=l;j<nCols;j++)
	  v[j][i] = (a[i][j] / a[i][l]) / g;
	for(j=l;j<nCols;j++) {
	  for(s=0.0,k=l;k<nCols;k++) s += a[i][k] * v[k][j];
	  for(k=l;k<nCols;k++) v[k][j] += s * v[k][i];
	}
      }
      for(j=l;j<nCols;j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
    fflush(stdout);
  }

  for(i=IMIN(nRows,nCols) - 1;i >= 0;i--) {
    l = i + 1;
    g = w[i];
    for(j=l;j<nCols;j++) a[i][j] = 0.0;
    if(g) {
      g = 1.0 / g;
      for(j=l;j<nCols;j++) {
	for(s=0.0,k=l;k<nRows;k++) s += a[k][i] * a[k][j];
	f = (s / a[i][i]) * g;
	for(k=i;k<nRows;k++) a[k][j] += f * a[k][i];
      }
      for(j=i;j<nRows;j++) a[j][i] *= g;
    }
    else
      for(j=i;j<nRows;j++) a[j][i] = 0.0;
    ++a[i][i];
    printf(".");
    fflush(stdout);
  }

  for(k=nCols-1;k>=0;k--) {
    for(its=0;its<30;its++) {
      flag = 1;
      for(l=k;l>=0;l--) {
	nm = l-1;
	if((fabs(rv1[l]) + anorm) == anorm) {
	  flag =  0;
	  break;
	}
	if((fabs(w[nm]) + anorm) == anorm) break;
      }
      if(flag) {
	c = 0.0;
	s = 1.0;
	for(i=l;i<=k;i++) {
	  f = s * rv1[i];
	  rv1[i] = c * rv1[i];
	  if((fabs(f) + anorm) == anorm) break;
	  g = w[i];
	  h = pythag(f,g);
	  w[i] = h;
	  h = 1.0 / h;
	  c = g * h;
	  s = -f * h;
	  for(j=0;j<nRows;j++) {
	    y = a[j][nm];
	    z = a[j][i];
	    a[j][nm] = y * c + z * s;
	    a[j][i] = z * c - y * s;
	  }
	}
      }
      z = w[k];
      if(l == k) {
	if(z < 0.0) {
	  w[k] = -z;
	  for(j=0;j<nCols;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if(its == 29) printf("no convergence in 30 svdcmp iterations\n");
      x = w[l];
      nm = k-1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = pythag(f,1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f))) - h)) / x;
      c = s = 1.0;
      for(j=l;j<=nm;j++) {
	i = j+1;
	g = rv1[i];
	y = w[i];
	h = s * g;
	g = c * g;
	z = pythag(f,h);
	rv1[j] = z;
	c = f/z;
	s = h/z;
	f = x * c + g * s;
	g = g * c - x * s;
	h = y * s;
	y *= c;
	for(jj=0;jj<nCols;jj++) {
	  x = v[jj][j];
	  z = v[jj][i];
	  v[jj][j] = x * c + z * s;
	  v[jj][i] = z * c - x * s;
	}
	z = pythag(f,h);
	w[j] = z;
	if(z) {
	  z = 1.0 / z;
	  c = f * z;
	  s = h * z;
	}
	f = c * g + s * y;
	x = c * y - s * g;
	for(jj=0;jj < nRows;jj++) {
	  y = a[jj][j];
	  z = a[jj][i];
	  a[jj][j] = y * c + z * s;
	  a[jj][i] = z * c - y * s;
	}
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
    printf(".");
    fflush(stdout);
  }
  
  free(rv1);
  
  return(0);
}


void svdbksb(double **u, double *w, double **v, int nRows, int nCols, double *b, double *x);
void svdbksb(double **u, double *w, double **v, int nRows, int nCols, double *b, double *x)
{
  int jj,j,i;
  double s;
  double *tmp;

  tmp = malloc(sizeof(double) * nCols);

  for(j=0;j<nCols;j++) { /* multiply b by U transpose */
    s = 0.0;
    if(w[j]) {
      for(i=0;i<nRows;i++)
	s += u[i][j] * b[i];
      s /= w[j];
    }
    tmp[j] = s;
  }

  for(j=0;j<nCols;j++) {
    s = 0.0;
    for(jj=0;jj<nCols;jj++)
      s += v[j][jj] * tmp[jj];
    x[j] = s;
  }

  free(tmp);
}
#undef IMIN
#undef FMAX
#undef SQRT
#undef SIGN
