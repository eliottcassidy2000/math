/* Auditor's own test of HYP-9669 on one tournament: all 2^n switchings S of V (no vertex fixed), all nonempty U.
   usage: audit_2adic n kind [param]
     kind: rand <seed> | trans | neartrans | circ <connection-bitmask over 1..n-1> | paley (n prime = 3 mod 4)
           | paleysrc (n-1 prime = 3 mod 4: Paley on 0..n-2 plus vertex n-1 as a source)
   Output: H, hc, violations of v2(c_U(H)) >= ceil(|U|/2) and v2(c_U(hc)) >= ceil(|U|/2)-1, min v2 by |U|,
           and the F2-degree of binary digits 1..4 of H on the cube F_2^V. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n, full; static int A[16][16]; /* A[i][j]=1 iff i->j */
static uint64_t g[1<<14][14];
static void counts(int out[], long long *H, long long *hc){
  memset(g,0,sizeof(uint64_t)*(full+1)*14);
  for(int v=0;v<n;v++) g[1<<v][v]=1;
  for(int W=1;W<=full;W++) for(int v=0;v<n;v++){ uint64_t c=g[W][v]; if(!c) continue; int nx=out[v]&~W; while(nx){int w=__builtin_ctz(nx); nx&=nx-1; g[W|(1<<w)][w]+=c;} }
  long long h=0; for(int v=0;v<n;v++) h+=g[full][v]; *H=h;
  memset(g,0,sizeof(uint64_t)*(full+1)*14); g[1][0]=1;
  for(int W=1;W<=full;W+=2) for(int v=0;v<n;v++){ uint64_t c=g[W][v]; if(!c) continue; int nx=out[v]&~W; while(nx){int w=__builtin_ctz(nx); nx&=nx-1; g[W|(1<<w)][w]+=c;} }
  long long c=0; for(int v=1;v<n;v++) if((out[v]>>0)&1) c+=g[full][v]; *hc=c;
}
static int v2(long long x){ if(!x) return 99; int k=0; while(!(x&1)){x>>=1;k++;} return k; }
static unsigned long long rs=88172645463325252ULL; static unsigned long long rnd(){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
int main(int argc,char**argv){
  n=atoi(argv[1]); full=(1<<n)-1; const char*k=argv[2];
  memset(A,0,sizeof(A));
  if(!strcmp(k,"rand")){ rs^=(unsigned long long)atoll(argv[3])*0x9E3779B97F4A7C15ULL; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(rnd()&1) A[i][j]=1; else A[j][i]=1; } }
  else if(!strcmp(k,"trans")||!strcmp(k,"neartrans")){ for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) A[i][j]=1; if(!strcmp(k,"neartrans")){ A[0][n-1]=0; A[n-1][0]=1; } }
  else if(!strcmp(k,"circ")){ int C=atoi(argv[3]); for(int i=0;i<n;i++) for(int d=1;d<n;d++) if((C>>d)&1) A[i][(i+d)%n]=1; }
  else if(!strcmp(k,"paley")||!strcmp(k,"paleysrc")){ int p = !strcmp(k,"paley")? n : n-1; int qr[32]={0}; for(int x=1;x<p;x++) qr[(x*x)%p]=1;
      for(int i=0;i<p;i++) for(int j=0;j<p;j++) if(i!=j && qr[((j-i)%p+p)%p]) A[i][j]=1;
      if(p==n-1) for(int j=0;j<p;j++) A[n-1][j]=1; }
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) if(A[i][j]+A[j][i]!=1){ fprintf(stderr,"not a tournament %d %d\n",i,j); return 1; }
  long long *cH=malloc(sizeof(long long)*(full+1)), *cC=malloc(sizeof(long long)*(full+1)), *Hv=malloc(sizeof(long long)*(full+1));
  for(int S=0;S<=full;S++){ int out[16]; for(int i=0;i<n;i++){ out[i]=0; for(int j=0;j<n;j++) if(j!=i){ int a=A[i][j]; if(((S>>i)&1)!=((S>>j)&1)) a=!a; if(a) out[i]|=1<<j; } }
    counts(out,&cH[S],&cC[S]); Hv[S]=cH[S]; }
  long long H0=cH[0], hc0=cC[0];
  for(int i=0;i<n;i++) for(int S=0;S<=full;S++) if((S>>i)&1){ cH[S]-=cH[S^(1<<i)]; cC[S]-=cC[S^(1<<i)]; }
  int mH[20],mC[20]; long long vH=0,vC=0; for(int i=0;i<20;i++){mH[i]=99;mC[i]=99;}
  for(int U=1;U<=full;U++){ int s=__builtin_popcount(U); int a=v2(cH[U]),b=v2(cC[U]); if(a<mH[s])mH[s]=a; if(b<mC[s])mC[s]=b; if(a<(s+1)/2) vH++; if(b<(s+1)/2-1) vC++; }
  printf("n=%d %s %s H=%lld hc=%lld | viol H %lld hc %lld | minv2 H:",n,k,argc>3?argv[3]:"",H0,hc0,vH,vC);
  for(int s=1;s<=n;s++) printf(" %d:%d",s,mH[s]); printf(" | hc:"); for(int s=1;s<=n;s++) printf(" %d:%d",s,mC[s]);
  /* digit degrees: F2 Moebius of digit k */
  printf(" | digit deg:");
  uint8_t *d=malloc(full+1);
  for(int kk=1;kk<=4;kk++){ for(int S=0;S<=full;S++) d[S]=(Hv[S]>>kk)&1; for(int i=0;i<n;i++) for(int S=0;S<=full;S++) if((S>>i)&1) d[S]^=d[S^(1<<i)];
    int deg=-1; for(int S=0;S<=full;S++) if(d[S] && __builtin_popcount(S)>deg) deg=__builtin_popcount(S); printf(" k%d:%d",kk,deg); }
  printf("\n"); return 0;
}
