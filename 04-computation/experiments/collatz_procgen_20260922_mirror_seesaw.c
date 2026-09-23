// collatz_procgen_20260922_mirror_seesaw.c
// Q1 mirror of loops_seesaw.c: level-r forward certificates of the E-game on classes mod 2^r.
//   f_b(y) = least number of x3-moves over E-paths from y making exactly b halvings (decided by y mod 2^b),
//   recursion (dimension lane, sec. 1.1): f_0 = 0; odd y: f_p(y) = 1 + f_p(3y+1);
//   even y: f_p(y) = min_j [2j + f_{p-1}(y_j/2)], y_0 = y, y_{j+1} = 9 y_j + 4 (mod 2^p).
//   Best certified factor of a class x mod 2^r: F(x) = min_{b<=r} 3^{f_b(x)}/2^b (a certificate if < 1).
//   Exceptional at level r: f_b(x) > floor(b log_3 2) for every b <= r.
// Output: exceptional count; worst certified factor; histogram of factors; worst factor by distance d to the
//   projected exceptional set; for the -1 thread exit at precision m (price P(m) = 3^{A*}/2^{m+1}, A* given
//   on the command line as a list), the number of landing classes y mod 2^r whose best factor is >= 1/P(m).
// Usage: mirror_seesaw r "m:A* m:A* ..."
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
static uint8_t **T;   // T[p][y/2] = f_p(y) for even y mod 2^p
static inline int fP(int p, uint64_t y){
  if(p==0) return 0; uint64_t m=(1ULL<<p)-1; y&=m;
  if(y&1) return 1 + T[p][((3*y+1)&m)>>1]; return T[p][y>>1];
}
int main(int argc,char**argv){
  int r=atoi(argv[1]); const char *mlist=argc>2?argv[2]:"";
  T=calloc(r+1,sizeof(uint8_t*));
  T[1]=malloc(1); T[1][0]=0;
  for(int p=2;p<=r;p++){
    uint64_t n=1ULL<<(p-1); T[p]=malloc(n); uint64_t mp=(1ULL<<p)-1;
    for(uint64_t i=0;i<n;i++){ uint64_t y=2*i; int best=255;
      for(int j=0;2*j<best;j++){ int v=fP(p-1,y>>1)+2*j; if(v<best) best=v; y=(9*y+4)&mp; }
      T[p][i]=best; }
  }
  double L3=log(3.0), L2=log(2.0), TH=L2/L3;
  uint64_t NC=1ULL<<r;
  // best factor exponent pair per class: store log-factor as float and (a,b) of the best certificate
  float *lf=malloc(sizeof(float)*NC); uint8_t *exc=calloc(NC,1);
  uint64_t nexc=0;
  for(uint64_t x=0;x<NC;x++){
    double best=1e9; int isexc=1;
    for(int b=1;b<=r;b++){ int a=fP(b,x); double l=a*L3-b*L2; if(l<best) best=l; if(a<=(int)floor(b*TH)) isexc=0; }
    lf[x]=(float)best; exc[x]=isexc; nexc+=isexc;
  }
  printf("level r=%d: classes %llu, exceptional %llu\n",r,(unsigned long long)NC,(unsigned long long)nexc);
  // worst certified factor
  double worst=-1e9; uint64_t wx=0;
  for(uint64_t x=0;x<NC;x++) if(!exc[x] && lf[x]>worst){ worst=lf[x]; wx=x; }
  printf("worst certified factor: %.6f at class %llu mod 2^%d (= %lld as a signed residue)\n",exp(worst),(unsigned long long)wx,r,
         (long long)(wx>=NC/2? (long long)wx-(long long)NC : (long long)wx));
  // histogram
  double cuts[]={0.25,0.5,0.75,0.9,0.95,0.99,1.0};
  printf("certified classes with factor >= c:");
  for(int i=0;i<7;i++){ uint64_t c=0; for(uint64_t x=0;x<NC;x++) if(!exc[x] && lf[x]>=log(cuts[i])-1e-12) c++; printf("  %.2f:%llu",cuts[i],(unsigned long long)c); }
  printf("\n");
  // distance to projected exceptional set: d(x) = largest j with x mod 2^j = e mod 2^j for some exceptional e
  // compute via marking prefixes
  uint8_t *dist=calloc(NC,1);
  for(int j=1;j<=r;j++){
    uint64_t M=(1ULL<<j)-1; uint8_t *mark=calloc(1ULL<<j,1);
    for(uint64_t x=0;x<NC;x++) if(exc[x]) mark[x&M]=1;
    for(uint64_t x=0;x<NC;x++) if(mark[x&M]) dist[x]=j;
    free(mark);
  }
  printf("worst certified factor by distance d to the exceptional set (d = shared low bits):\n");
  for(int j=0;j<=r;j++){
    double w=-1e9; uint64_t c=0; for(uint64_t x=0;x<NC;x++) if(!exc[x] && dist[x]==j){ c++; if(lf[x]>w) w=lf[x]; }
    if(c) printf("  d=%2d: %9llu classes, worst factor %.6f\n",j,(unsigned long long)c,exp(w));
  }
  // landing-class failures for the -1 exit at precision m
  printf("-1 thread exit at precision m: price P(m)=3^A*/2^(m+1); landing y is uniform mod 2^r;\n");
  printf("  fail = #classes y mod 2^%d with best factor >= 1/P(m) (incl. exceptional)\n",r);
  const char *s=mlist;
  while(*s){
    int m,A; int nch=0; if(sscanf(s,"%d:%d%n",&m,&A,&nch)!=2) break; s+=nch; while(*s==' ') s++;
    double lp=A*L3-(m+1)*L2;   // log P
    uint64_t fail=0; for(uint64_t x=0;x<NC;x++) if(exc[x] || lf[x] >= -lp - 1e-12) fail++;
    printf("  m=%3d A*=%3d P=%.6f 1/P=%.6f fail=%9llu (%.3e)\n",m,A,exp(lp),exp(-lp),(unsigned long long)fail,(double)fail/NC);
  }
  return 0;
}
