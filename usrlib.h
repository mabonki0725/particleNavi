#define FLOW 1
#define REVS 0
#define PAI 3.14159265358

char *comMalloc(int);
char *comAssign(char **,char *);
int  comCmember(char **,int,char *);
int  comNmember(int *,int,int);
int  comDsort(double *,int,int);
int  comDsortJun(double *,int,int,int *);
void **comMxAlloc(int,int,int);
void comMxFree(void **,int,int);
double comQnorm(double);
double comPnorm(double);