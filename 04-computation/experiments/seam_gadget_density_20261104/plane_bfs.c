/* Pure plane dynamics on P = {S = 0} (HYP-9471 P6): BFS of (H, R) under the 2x2 restrictions [[a,b],[c,d]]
   (a = M[0][0], b = M[0][2], c = M[2][0], d = M[2][2]) of every S-row-zero matrix in the input list.
   Seeds: the twelve realisable states with S = 0 of order <= 5. States with H <= X and R <= X are kept.
   Prints the realised fraction of lattice points (H odd, R even) in two bands of slope [2, 2.75]. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
int main(int argc,char**argv){ ll X=atoll(argv[1]); FILE*f=fopen(argv[2],"r");
  static ll A[600000][4]; int n=0; ll m[9];
  while(fscanf(f,"%lld %lld %lld %lld %lld %lld %lld %lld %lld",&m[0],&m[1],&m[2],&m[3],&m[4],&m[5],&m[6],&m[7],&m[8])==9){
    if(m[3]||m[4]||m[5]) continue; if(m[0]*m[8]-m[2]*m[6]==0) continue;
    int dup=0; /* cheap dedupe of consecutive duplicates only; exact dedupe below */
    A[n][0]=m[0];A[n][1]=m[2];A[n][2]=m[6];A[n][3]=m[8]; n++; if(n>=600000) break; }
  fclose(f);
  /* exact dedupe */
  int cmp(const void*p,const void*q){ const ll*x=p,*y=q; for(int i=0;i<4;i++){ if(x[i]<y[i]) return -1; if(x[i]>y[i]) return 1;} return 0; }
  qsort(A,n,sizeof(A[0]),cmp); int k=0; for(int i=0;i<n;i++){ if(k==0||cmp(A[i],A[k-1])) { memcpy(A[k],A[i],sizeof(A[0])); k++; } } n=k;
  size_t W=(size_t)(X+1)*(size_t)(X/2+1); unsigned char*seen=calloc(W/8+1,1);
  #define IDX(h,r) ((size_t)(h)*(size_t)(X/2+1)+(size_t)((r)/2))
  ll *qh=malloc(sizeof(ll)*2*20000000); size_t qs=0,qe=0;
  ll seeds[12][2]={{1,2},{1,4},{1,8},{1,16},{3,6},{3,12},{3,24},{5,10},{5,12},{5,14},{5,18},{5,32}};
  for(int i=0;i<12;i++){ ll h=seeds[i][0],r=seeds[i][1]; size_t id=IDX(h,r); if(!(seen[id>>3]>>(id&7)&1)){ seen[id>>3]|=1<<(id&7); qh[2*qe]=h; qh[2*qe+1]=r; qe++; } }
  while(qs<qe){ ll h=qh[2*qs], r=qh[2*qs+1]; qs++;
    for(int i=0;i<n;i++){ ll h2=A[i][0]*h+A[i][1]*r; if(h2>X) continue; ll r2=A[i][2]*h+A[i][3]*r; if(r2>X||r2<0) continue;
      if(r2&1){ fprintf(stderr,"odd R'\n"); return 1; }
      size_t id=IDX(h2,r2); if(seen[id>>3]>>(id&7)&1) continue; seen[id>>3]|=1<<(id&7);
      if(qe>=20000000){ fprintf(stderr,"queue full\n"); return 1;} qh[2*qe]=h2; qh[2*qe+1]=r2; qe++; } }
  double fill[2]; ll lo[2]={801,1601}, hi[2]={999,1799};
  for(int b=0;b<2;b++){ ll tot=0,have=0; for(ll h=lo[b];h<=hi[b];h+=2){ ll r0=2*h; ll r1=(ll)(2.75*h); for(ll r=r0;r<=r1;r+=2){ tot++; size_t id=IDX(h,r); have+=seen[id>>3]>>(id&7)&1; } } fill[b]=(double)have/tot; }
  printf("X=%lld maps=%d states=%zu fill[801,999]=%.4f fill[1601,1799]=%.4f ratio=%.4f\n",X,n,qe,fill[0],fill[1],fill[1]/fill[0]);
  return 0; }
