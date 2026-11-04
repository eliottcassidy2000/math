// coverage_stream.c: stream "h e B strong" rows (allpairs output) and mark attained values in bitsets:
// two-step h+2B, near-join h+2^k B (k<=8), tower heights 1..10 (h A_k - e C_k + B D_k); values above CAP ignored.
// Report, for each family, the largest X with every odd in [69, X] attained, and misses below 69.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define CAP (1ULL<<31)
static unsigned char *two,*nj,*tw;
static inline void mark(unsigned char*b,u64 v){ if(v<CAP) b[v>>3]|=1<<(v&7); }
static inline int get(unsigned char*b,u64 v){ return v<CAP && (b[v>>3]>>(v&7)&1); }
int main(void){ size_t sz=CAP/8; two=calloc(sz,1); nj=calloc(sz,1); tw=calloc(sz,1);
  long long A[11]={1,1,1},C[11]={0,1,0},D[11]={0,1,2}; for(int k=3;k<=10;k++){A[k]=A[k-1]+A[k-2]+A[k-3];C[k]=C[k-1]+C[k-2]+C[k-3];D[k]=D[k-1]+D[k-2]+D[k-3];}
  u64 h,e,B; int st; u64 rows=0;
  while(scanf("%llu %llu %llu %d",&h,&e,&B,&st)==4){ rows++; mark(two,h+2*B); for(int k=1;k<=8;k++) mark(nj,h+(B<<k)); for(int k=1;k<=10;k++){ long long v=(long long)h*A[k]-(long long)e*C[k]+(long long)B*D[k]; if(v>0) mark(tw,(u64)v);} }
  fprintf(stderr,"rows %llu\n",rows);
  const char*names[4]={"two-step h+2B","near-join h+2^k B (k<=8)","tower heights 1..10","union"};
  for(int f=0;f<4;f++){ unsigned char*b; u64 X=69; 
    #define IN(v) (f==0?get(two,v):f==1?get(nj,v):f==2?get(tw,v):(get(two,v)||get(nj,v)||get(tw,v)))
    while(X+2<CAP && (IN(X+2)|| X+2==7 || X+2==21)) X+=2;
    printf("%-28s: every odd in [69, %llu] attained; missed below 69:",names[f],X);
    for(u64 v=1;v<69;v+=2) if(!IN(v) && v!=7 && v!=21) printf(" %llu",v);
    printf("; next gaps:"); int c=0; for(u64 v=X+2; v<X+4000 && c<6; v+=2) if(!IN(v)){printf(" %llu",v);c++;} printf("\n"); }
  return 0; }
