#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "usrlib.h"

/**********************
  共通メモリー割り当て
***********************/
char *comMalloc(size)
int size;
{
   char *p;
   p=(char *)malloc(size);
   if(p) memset(p,'\0',size); 
   else  printf("malloc0 領域は確保出来ません!!\n");
   return(p);
}
/********************
  共通文字割り当て
*********************/
char *comAssign(adr,str)
char **adr;
char *str;
{
    int len;

    len=strlen(str);
    if(!len) return(NULL);

    *adr=(char *)malloc(len+1);
    memset(*adr,'\0',len+1);
    strcpy(*adr,str);
    return(*adr);
}
/*********************
  文字列のリスト
**********************/
int comCmember(memb,num,pc)
char **memb;
int num;
char *pc;
{
   int i;
   for(i=0;i<num;i++) {
     if(!strcmp(memb[i],pc)) break;
   }
   if(i == num) return(-1);
   return(i);
}
/*********************
  文字列のリスト
**********************/
int comNmember(memb,num,id)
int *memb;
int num;
int id;
{
   int i;
   for(i=0;i<num;i++) {
     if(memb[i] == id) break;
   }
   if(i == num) return(-1);
   return(i);
}
/**************
数値のソート
***************/
int comDsort(v,order,n)
double v[];
int  order;
int  n;
{
   int gap,i,j;
   double dwrk;

   for(gap = n/2;gap > 0;gap /= 2) {
     for(i = gap;i < n; i++) { 
       for(j = i-gap; j >= 0; j -= gap) {
    	 if((order == FLOW && v[j] <= v[j+gap]) ||
	        (order != FLOW && v[j] >= v[j+gap])) break;
	       dwrk=v[j];
	       v[j]=v[j+gap];
	       v[j+gap]=dwrk;
       } 
     }
   }
   return(n);
}
/**************
数値のソート順番
***************/
int comDsortJun(org,order,n,jun)
double org[];
int  order;
int  n;
int  jun[];
{
   int gap,i,j;
   double dwrk;
   double *v;
   int jwrk;

   v=(double *)comMalloc(sizeof(double)*n);
   for(i=0;i<n;i++) v[i]=org[i];

   for(gap = n/2;gap > 0;gap /= 2) {
     for(i = gap;i < n; i++) { 
       for(j = i-gap; j >= 0; j -= gap) {
    	 if((order == FLOW && v[j] <= v[j+gap]) ||
	        (order != FLOW && v[j] >= v[j+gap])) break;
	       dwrk=v[j];
	       v[j]=v[j+gap];
	       v[j+gap]=dwrk;

           jwrk=jun[j];
           jun[j]=jun[j+gap];
           jun[j+gap]=jwrk;
       } 
     }
   }
   free(v);
   return(n);
}
/************************
  ２次元配列の確保
*************************/
void **comMxAlloc(int n,int m,int size)
{
    int i;
    void **data;

    data=(void **)comMalloc(sizeof(void *)*n);
    for(i=0;i<n;i++) {
      data[i]=(void *)comMalloc(size*m);
    }
    return(data);
}
void comMxFree(void **data,int n,int m)
{
    int i;
    for(i=0;i<n;i++) {
      free(data[i]);
    }
    free(data);
}
/*************
  逆正規分布
**************/
double comQnorm(double u)
{
	static double a[9] = {	 1.24818987e-4, -1.075204047e-3, 5.198775019e-3,
							-0.019198292004, 0.059054035642,-0.151968751364,
							 0.319152932694,-0.5319230073,   0.797884560593};
	static double b[15] = {	-4.5255659e-5,   1.5252929e-4,  -1.9538132e-5,
							-6.76904986e-4,  1.390604284e-3,-7.9462082e-4,
							-2.034254874e-3, 6.549791214e-3,-0.010557625006,
							 0.011630447319,-9.279453341e-3, 5.353579108e-3,
							-2.141268741e-3, 5.35310549e-4,  0.999936657524};
	double w, y, z;
	int i;

	if(u == 0.)	return 0.5;
	y = u / 2.;
	if(y < -6.)	return 0.;
	if(y > 6.)		return 1.;
	if(y < 0.)		y = - y;
	if(y < 1.)
	{
		w = y * y;
		z = a[0];
		for(i = 1; i < 9; i++)		z = z * w + a[i];
		z *= (y * 2.);
	}
	else
	{
		y -= 2.;
		z = b[0];
		for(i = 1; i < 15; i++)	z = z * y + b[i];
	}

	if(u < 0.)	return (1. - z) / 2.;
	return (1. + z) / 2.;
}
/****************
   正規分布
*****************/
double comPnorm(double qn)
{
	static double b[11] = {	 1.570796288,     0.03706987906,  -0.8364353589e-3,
							-0.2250947176e-3, 0.6841218299e-5, 0.5824238515e-5,
							-0.104527497e-5,  0.8360937017e-7,-0.3231081277e-8,
							 0.3657763036e-10,0.6936233982e-12};
	double w1, w3;
	int i;

	if(qn < 0. || 1. < qn)
	{
		fprintf(stderr, "Error : qn <= 0 or qn >= 1  in pnorm()!\n");
		return 0.;
	}
	if(qn == 0.5)	return 0.;

	w1 = qn;
	if(qn > 0.5)	w1 = 1. - w1;
	w3 = -log(4. * w1 * (1. - w1));
	w1 = b[0];
	for(i = 1; i < 11; i++)	w1 += (b[i] * pow(w3, (double)i));
	if(qn > 0.5)	return sqrt(w1 * w3);
	return -sqrt(w1 * w3);
}


