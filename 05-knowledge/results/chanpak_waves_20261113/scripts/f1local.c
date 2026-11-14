// Wave 23: is the excess of F1 misses on a square class t*m^2 explained by local factors?
// Local model: for primes 3 <= p <= PMAX, class means of r(N) over the window; predicted log-mass
// lm(N) = sum_p log(mean_p[N mod p] / mean); bin all odd N of the window by lm into NB equal-count bins;
// expected misses of a subset = sum over its members of the miss rate of their bin.  Compare with observed.
// Usage: f1local K lo PMAX
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef long long i64; typedef unsigned long long u64;
static int cmpf(const void*a,const void*b){ float x=*(const float*)a,y=*(const float*)b; return (x>y)-(x<y); }
static int isprime(i64 n){ if(n<2) return 0; for(i64 d=2;d*d<=n;d++) if(n%d==0) return 0; return 1; }
int main(int argc,char**argv){
  int K=argc>1?atoi(argv[1]):27, lo=argc>2?atoi(argv[2]):22, PMAX=argc>3?atoi(argv[3]):97;
  i64 X=1LL<<K, L=1LL<<lo, Hmax=X/2+2;
  uint32_t *psi=malloc(Hmax*sizeof(uint32_t)); psi[0]=0;
  for(i64 h=1;h<Hmax;h++){ i64 q=h>>1; psi[h]=(h&1)? psi[q]+(uint32_t)(q>>1)+1 : psi[q]+(uint32_t)q; }
  uint8_t *r=calloc(X+1,1);
  for(int s=0;s<K;s++) for(i64 h1=0;h1<Hmax;h1++){
    i64 y1=(2*h1+1)<<s; if(y1>X) break; i64 g1=2*(i64)psi[h1]+1; i64 lim=X/y1; i64 h2max=(lim-1)/2; if(h2max<0) break;
    for(i64 h2=0;h2<=h2max;h2++){ i64 N=g1*(2*(i64)psi[h2]+1)+y1*(2*h2+1)-(s==0); if(N<=X && r[N]<255) r[N]++; } }
  free(psi);
  // primes
  int P[64],np=0; for(int p=3;p<=PMAX;p++) if(isprime(p)) P[np++]=p;
  double *cls[64]; u64 *cnt[64];
  double tot=0; u64 n=0;
  for(int i=0;i<np;i++){ cls[i]=calloc(P[i],sizeof(double)); cnt[i]=calloc(P[i],sizeof(u64)); }
  for(i64 N=L|1;N<X;N+=2){ tot+=r[N]; n++; for(int i=0;i<np;i++){ int c=N%P[i]; cls[i][c]+=r[N]; cnt[i][c]++; } }
  double mean=tot/n;
  float *lt[64]; for(int i=0;i<np;i++){ lt[i]=malloc(P[i]*sizeof(float)); for(int c=0;c<P[i];c++) lt[i][c]=(float)log((cls[i][c]/cnt[i][c])/mean); }
  // predicted log mass and equal-count bins
  u64 W=n; float *lm=malloc(W*sizeof(float)); u64 k=0;
  for(i64 N=L|1;N<X;N+=2){ float s=0; for(int i=0;i<np;i++) s+=lt[i][N%P[i]]; lm[k++]=s; }
  float *srt=malloc(W*sizeof(float)); memcpy(srt,lm,W*sizeof(float)); qsort(srt,W,sizeof(float),cmpf);
  const int NB=100; float edge[NB+1]; for(int b=0;b<NB;b++) edge[b]=srt[(u64)((double)b*W/NB)]; edge[NB]=1e30f; free(srt);
  u64 bn[NB]={0}, bm[NB]={0};
  #define BIN(v) ({ int lo_=0,hi_=NB-1; while(lo_<hi_){ int md=(lo_+hi_+1)/2; if((v)>=edge[md]) lo_=md; else hi_=md-1; } lo_; })
  k=0; for(i64 N=L|1;N<X;N+=2){ int b=BIN(lm[k]); bn[b]++; if(!r[N]) bm[b]++; k++; }
  printf("F1 window odd N in [2^%d,2^%d): n=%llu misses rate %.5f; local model over %d primes 3..%d, %d bins\n",lo,K,n,0.0,np,PMAX,NB);
  { u64 mm=0; for(int b=0;b<NB;b++) mm+=bm[b]; printf("  overall miss rate %.5f; bin miss rates: lowest-mass bin %.5f, median bin %.5f, highest-mass bin %.5f\n",
      (double)mm/n,(double)bm[0]/bn[0],(double)bm[NB/2]/bn[NB/2],(double)bm[NB-1]/bn[NB-1]); }
  // subsets
  int ts[]={1,3,5,7,11,13,15,21,35};
  for(int ti=0;ti<9;ti++){ i64 t=ts[ti]; u64 cnt_=0,obs=0; double expd=0, expflat=0;
    for(i64 m=1;;m+=2){ i64 N=t*m*m; if(N>=X) break; if(N<L) continue; i64 a=m,b=t; while(b){i64 c=a%b;a=b;b=c;} if(t>1&&a!=1) continue;
      u64 idx=(N-(L|1))/2; int bb=BIN(lm[idx]); cnt_++; expd+=(double)bm[bb]/bn[bb]; if(!r[N]) obs++; }
    printf("  t*m^2, t=%3lld: n=%6llu observed misses %4llu, expected under local model %.1f, ratio %.2f\n",t,cnt_,obs,expd,expd>0?obs/expd:-1.0); }
  // primes and Rabinowitsch-rough as sanity checks
  { u64 cnt_=0,obs=0; double expd=0; u64 c2=0,o2=0; double e2=0;
    for(i64 N=L|1;N<X;N+=2){ if(N%3==0||N%5==0||N%7==0) continue; if(!isprime(N)) continue; u64 idx=(N-(L|1))/2; int bb=BIN(lm[idx]);
      cnt_++; expd+=(double)bm[bb]/bn[bb]; if(!r[N]) obs++;
      if((N&3)==3 && isprime((N+1)/4)){ c2++; e2+=(double)bm[bb]/bn[bb]; if(!r[N]) o2++; }
      if(cnt_>=400000) break; }
    printf("  primes (first %llu of the window): observed %llu, expected %.1f, ratio %.2f\n",cnt_,obs,expd,obs/expd);
    printf("  Rabinowitsch-rough primes among them: n=%llu observed %llu expected %.1f ratio %.2f\n",c2,o2,e2,e2>0?o2/e2:-1.0); }
  return 0; }
