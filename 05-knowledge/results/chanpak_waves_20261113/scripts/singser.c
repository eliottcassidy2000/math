// Singular-series check for the strict 11-block orbit: over odd N in [X/2, X], the ratio
//   rho_p = mean r over N not divisible by p  /  mean r over all N
// for small primes p, the product of these ratios, and the observed mean over primes / overall mean.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
static i64 X; static unsigned *cnt; static int L0, NB, EMPTYOK; static unsigned PB[64]; static int PLEN[64];
static int isprefix(unsigned code,int len){ for(int b=0;b<NB;b++) if(PLEN[b]>=len && (PB[b]&((1u<<len)-1))==code) return 1; return 0; }
static int isterm(unsigned code,int len){ if(len==0) return EMPTYOK; for(int b=0;b<NB;b++) if(PLEN[b]==len && PB[b]==code) return 1; return 0; }
static void dfs(i64 H,i64 E,i64 L,i64 h,int run,unsigned code,int init){
  int term = init || isterm(code,run); if(term) cnt[H>>1]++;
  i64 D=L-E; if(D>X-H) return;
  if(init ? run<L0 : isprefix(code,run+1)) { if(2*L<=X-H) dfs(H,0,2*L,h,run+1,code,init); }
  { unsigned c2=code|(1u<<run); if((init ? run<L0 : isprefix(c2,run+1)) && !(D==0 && L==h && E==h)) dfs(H+D,h,2*L-h,h,run+1,c2,init); }
  if(term) dfs(H+D,H-E,H+L,H,0,0,0);
}
int main(int argc,char**argv){
  X=atoll(argv[1]); L0=atoi(argv[2]); NB=0; EMPTYOK=0;
  for(int i=3;i<argc;i++){ const char*w=argv[i]; if(!strcmp(w,"e")){ EMPTYOK=1; continue; } unsigned c=0; int l=(int)strlen(w); for(int t=0;t<l;t++) if(w[t]=='Z') c|=1u<<t; PB[NB]=c; PLEN[NB]=l; NB++; }
  cnt=calloc(X/2+2,sizeof(unsigned)); unsigned char*comp=calloc(X/2+2,1);
  for(i64 p=3;p*p<=X;p+=2) if(!comp[p>>1]) for(i64 q=p*p;q<=X;q+=2*p) comp[q>>1]=1;
  dfs(1,1,1,1,0,0,1);
  int ps[]={3,5,7,11,13,17,19,23,29,31,37,41,43}; int np=13;
  double all=0, nall=0, pr=0, npr=0, rough=0, nrough=0;
  double s[13]={0}, ns[13]={0};
  for(i64 n=(X/2)|1;n<=X;n+=2){ double r=cnt[n>>1]; all+=r; nall++;
    if(!comp[n>>1]){ pr+=r; npr++; }
    int isrough=1; for(int k=0;k<np;k++){ if(n%ps[k]){ s[k]+=r; ns[k]++; } else isrough=0; }
    if(isrough){ rough+=r; nrough++; } }
  double mean=all/nall, prod=1;
  printf("X=%lld: overall mean %.1f; primes mean %.1f (ratio %.3f); N coprime to all p<=43: mean %.1f (ratio %.3f)\n",X,mean,pr/npr,pr/npr/mean,rough/nrough,rough/nrough/mean);
  for(int k=0;k<np;k++){ double rr=(s[k]/ns[k])/mean; prod*=rr; printf("  p=%2d: coprime-class ratio %.4f   cumulative product %.4f\n",ps[k],rr,prod); }
  return 0; }
