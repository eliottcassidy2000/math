/* closure4.c - four-move closure (cone, E_2-join, cone over G u K_1, R) from exhaustive bases, for PROOF.md remark 7.1.
   usage: closure4 NM base9.txt NB  (base9.txt = n d U pendant pendantNU from mkbase.py columns 1,2,3,4,5) */
/* Value-level closure of even-graph realisations by order n, sets:
   E[n] all, W[n] no universal vertex, P[n] has a pendant vertex, S[n] pendant whose neighbour is non-universal.
   Bases: exhaustive for n <= NB (file). Rules for n > NB:
     E_n >= E_{n-1}, W_{n-1}+2^{n-1} (cone), E_{n-2}+3*2^{n-2}-2 (E_2-join)
     W_n >= E_{n-1}, W_{n-2}+3*2^{n-2}-2, S_n
     P_n >= P_{n-1}, S_n, E_{n-2}+2^{n-1} (cone over G u K_1)
     S_n >= S_{n-1}, P_{n-1}, S_{n-2}+3*2^{n-3} (R: hub misses the pendant, one leaf)
   Report misses of E_n against odd[1,M(n)] minus {3,7,13} minus X(n), and of W_n against odd < 2^(n-1). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define NMAXX 30
static uint8_t *E[NMAXX+1],*W[NMAXX+1],*P[NMAXX+1],*S[NMAXX+1];
static int NB,NM;
static inline void setb(uint8_t*a,uint64_t k){a[k>>4] |= 1<<((k>>1)&7);} /* odd k only */
static inline int getb(uint8_t*a,uint64_t k){return (a[k>>4]>>((k>>1)&7))&1;}
static uint64_t sz(int n){return (((uint64_t)1<<n)>>4)+1;}
static void orset(uint8_t*dst,uint8_t*src,int nsrc){ uint64_t m=sz(nsrc); for(uint64_t i=0;i<m;i++) dst[i]|=src[i]; }
static void shiftadd(uint8_t*dst,int ndst,uint8_t*src,int nsrc,int64_t off){
  uint64_t lim=(uint64_t)1<<nsrc, dl=(uint64_t)1<<ndst;
  for(uint64_t k=1;k<lim;k+=2) if(getb(src,k)){ int64_t v=(int64_t)k+off; if(v>0 && (uint64_t)v<dl) setb(dst,(uint64_t)v);} }
int main(int argc,char**argv){
  NM=atoi(argv[1]); NB=atoi(argv[3]);
  for(int n=1;n<=NM;n++){E[n]=calloc(sz(n),1);W[n]=calloc(sz(n),1);P[n]=calloc(sz(n),1);S[n]=calloc(sz(n),1);}
  FILE*f=fopen(argv[2],"r"); int n,U,pe,pnu; unsigned long long d;
  while(fscanf(f,"%d %llu %d %d %d",&n,&d,&U,&pe,&pnu)==5){ if(n>NB||n>NM) continue;
    setb(E[n],d); if(!U) setb(W[n],d); if(pe) setb(P[n],d); if(pnu) setb(S[n],d); }
  fclose(f);
  /* padding monotonicity within the base range */
  for(int n=2;n<=NB && n<=NM;n++){ orset(E[n],E[n-1],n-1); orset(W[n],E[n-1],n-1); orset(P[n],P[n-1],n-1); orset(S[n],S[n-1],n-1); orset(S[n],P[n-1],n-1);}
  for(int n=NB+1;n<=NM;n++){
    orset(S[n],S[n-1],n-1); orset(S[n],P[n-1],n-1); shiftadd(S[n],n,S[n-2],n-2,(int64_t)3<<(n-3));
    orset(P[n],P[n-1],n-1); orset(P[n],S[n],n); shiftadd(P[n],n,E[n-2],n-2,(int64_t)1<<(n-1));
    orset(W[n],E[n-1],n-1); shiftadd(W[n],n,W[n-2],n-2,((int64_t)3<<(n-2))-2); orset(W[n],S[n],n);
    orset(E[n],E[n-1],n-1); orset(E[n],W[n],n); shiftadd(E[n],n,W[n-1],n-1,(int64_t)1<<(n-1)); shiftadd(E[n],n,E[n-2],n-2,((int64_t)3<<(n-2))-2);
  }
  for(int n=3;n<=NM;n++){
    uint64_t M = (n&1)? ((uint64_t)1<<n)-n : ((uint64_t)1<<n)-n-1;
    long missE=0, extraE=0, missW=0; uint64_t firstmiss=0;
    for(uint64_t k=1;k<((uint64_t)1<<n);k+=2){
      int exc = (k==3||k==7||k==13) || ((n&1) && n>=5 && k==((uint64_t)1<<n)-n-2);
      int want = (k<=M) && !exc;
      int has = getb(E[n],k);
      if(want && !has){ missE++; if(!firstmiss) firstmiss=k; }
      if(!want && has) extraE++;
      if(k < ((uint64_t)1<<(n-1)) && !(k==3||k==7||k==13) && !getb(W[n],k)) missW++;
    }
    printf("n=%d M=%llu missE=%ld extraE=%ld firstmissE=%llu missW(<2^(n-1))=%ld\n",n,(unsigned long long)M,missE,extraE,(unsigned long long)firstmiss,missW);
    fflush(stdout);
  }
  return 0;
}
