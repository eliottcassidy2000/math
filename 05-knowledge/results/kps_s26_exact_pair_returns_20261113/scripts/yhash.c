/* yhash.c: enumerate strong order-n tournaments with <= K backward arcs (w.r.t. 0<1<..<n-1) and report those whose
   full path-cover profile lies in a (large) target set (file lines: c3 F_1..F_n), with their Hamiltonian-cycle count C.
   usage: yhash n K targets.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXN 12
static int n,K; static int pi[100],pj[100]; static int npairs;
typedef struct { u64 F[MAXN]; int used; } tent;
static tent *tt; static size_t ttsz;
static u64 *hset; static size_t hsz;
static u64 Hs[1<<MAXN]; static u64 endp[(1<<MAXN)*MAXN]; static u64 g[2][1<<MAXN]; static u64 st0[(1<<MAXN)*MAXN];
static long long tested=0,strongc=0,hpass=0,matches=0,matchesC1=0;
static u64 hF(const u64*F){ u64 h=1469598103934665603ULL; for(int j=0;j<n;j++){ h^=F[j]; h*=1099511628211ULL; } return h; }
static void tadd(const u64*F){ size_t k=hF(F)%ttsz; while(tt[k].used){ if(!memcmp(tt[k].F,F,sizeof(u64)*n)) return; k=(k+1)%ttsz; } tt[k].used=1; memcpy(tt[k].F,F,sizeof(u64)*n); }
static int tin(const u64*F){ size_t k=hF(F)%ttsz; while(tt[k].used){ if(!memcmp(tt[k].F,F,sizeof(u64)*n)) return 1; k=(k+1)%ttsz; } return 0; }
static int inH(u64 h){ size_t k=(h*0x9E3779B97F4A7C15ULL)%hsz; while(hset[k]){ if(hset[k]==h) return 1; k=(k+1)%hsz; } return 0; }
static void addH(u64 h){ size_t k=(h*0x9E3779B97F4A7C15ULL)%hsz; while(hset[k]){ if(hset[k]==h) return; k=(k+1)%hsz; } hset[k]=h; }
static int sel[64];
static void test(int k){
  unsigned out[MAXN]; for(int i=0;i<n;i++){ out[i]=0; for(int j=i+1;j<n;j++) out[i]|=1u<<j; }
  for(int a=0;a<k;a++){ int i=pi[sel[a]], j=pj[sel[a]]; out[i]&=~(1u<<j); out[j]|=1u<<i; }
  tested++;
  unsigned FULL=(1u<<n)-1; unsigned seen=1, fr=1;
  while(fr){ unsigned nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=out[v]; nf&=~seen; seen|=nf; fr=nf; }
  if(seen!=FULL) return;
  unsigned in[MAXN]={0}; for(int v=0;v<n;v++) for(int w=0;w<n;w++) if(out[v]>>w&1) in[w]|=1u<<v;
  seen=1; fr=1; while(fr){ unsigned nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=in[v]; nf&=~seen; seen|=nf; fr=nf; }
  if(seen!=FULL) return;
  strongc++;
  unsigned N=1u<<n;
  for(unsigned S=1;S<N;S++){ u64 tot=0;
    for(int v=0;v<n;v++) if(S>>v&1){ unsigned P=S&~(1u<<v); u64 a=0;
      if(!P) a=1; else for(int w=0;w<n;w++) if((P>>w&1)&&(out[w]>>v&1)) a+=endp[(size_t)P*n+w];
      endp[(size_t)S*n+v]=a; tot+=a; }
    Hs[S]=tot; }
  if(!inH(Hs[N-1])) return;
  hpass++;
  u64 F[MAXN]; F[0]=Hs[N-1];
  for(unsigned S=0;S<N;S++) g[0][S]=Hs[S];
  int cur=0;
  for(int j=2;j<=n;j++){ int nx=cur^1;
    for(unsigned S=0;S<N;S++){ u64 s=0; for(unsigned B=S;B;B=(B-1)&S){ unsigned A=S&~B; if(!A) continue; s+=g[cur][A]*Hs[B]; } g[nx][S]=s; }
    cur=nx; F[j-1]=g[cur][N-1]; }
  if(!tin(F)) return;
  /* Hamiltonian cycles */
  memset(st0,0,sizeof(u64)*N*n); st0[1u*n+0]=1;
  for(unsigned S=1;S<N;S++){ if(!(S&1)) continue;
    for(int v=0;v<n;v++){ u64 a=st0[(size_t)S*n+v]; if(!a) continue;
      for(int w=0;w<n;w++) if(!(S>>w&1)&&(out[v]>>w&1)) st0[(size_t)(S|1u<<w)*n+w]+=a; } }
  u64 C=0; for(int v=1;v<n;v++) if(out[v]&1u) C+=st0[(size_t)(N-1)*n+v];
  matches++; if(C==1) matchesC1++;
  printf("MATCH C %llu H %llu arcs",C,F[0]); for(int a=0;a<k;a++) printf(" %d>%d",pj[sel[a]],pi[sel[a]]);
  printf(" rows"); for(int v=0;v<n;v++) printf(" %u",out[v]); putchar(10); fflush(stdout);
}
static void rec(int start,int depth,int k){ if(depth==k){ test(k); return; } for(int p=start;p<npairs;p++){ sel[depth]=p; rec(p+1,depth+1,k); } }
int main(int argc,char**argv){
  n=atoi(argv[1]); K=atoi(argv[2]); FILE*f=fopen(argv[3],"r");
  ttsz=1<<18; tt=calloc(ttsz,sizeof(tent)); hsz=1<<18; hset=calloc(hsz,sizeof(u64));
  char line[8192]; long nt=0;
  while(fgets(line,sizeof line,f)){ char*p=line; strtol(p,&p,10); u64 F[MAXN]; for(int j=0;j<n;j++) F[j]=strtoull(p,&p,10); tadd(F); addH(F[0]); nt++; }
  fprintf(stderr,"targets %ld\n",nt);
  npairs=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ pi[npairs]=i; pj[npairs]=j; npairs++; }
  for(int k=1;k<=K;k++){ rec(0,0,k); fprintf(stderr,"k=%d done: tested %lld strong %lld Hpass %lld matches %lld (C=1: %lld)\n",k,tested,strongc,hpass,matches,matchesC1); fflush(stderr); }
  return 0;
}
