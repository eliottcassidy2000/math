/* ysearch.c: search order-n tournaments with <= K backward arcs (w.r.t. 0<1<..<n-1) whose full
   path-cover profile equals one of the target profiles (file: c3 F_1..F_n per line).
   usage: ysearch n K targets.txt   -> prints matches: K-set and the profile. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXN 13
static int n,K; static int pi[100],pj[100]; static int npairs;
static int ntg; static u64 tg[4000][MAXN]; static int tc3[4000];
static char c3ok[400]; static u64 *hset; static size_t hsz;
static u64 Hs[1<<MAXN]; static u64 endp[(1<<MAXN)*MAXN]; static u64 g[2][1<<MAXN];
static long long tested=0, strongc=0, c3pass=0, hpass=0, matches=0;
static int inH(u64 h){ size_t k=(h*0x9E3779B97F4A7C15ULL)%hsz; while(hset[k]){ if(hset[k]==h) return 1; k=(k+1)%hsz; } return 0; }
static void addH(u64 h){ size_t k=(h*0x9E3779B97F4A7C15ULL)%hsz; while(hset[k]){ if(hset[k]==h) return; k=(k+1)%hsz; } hset[k]=h; }
static int sel[20];
static void test(int k){
  unsigned out[MAXN]; for(int i=0;i<n;i++){ out[i]=0; for(int j=i+1;j<n;j++) out[i]|=1u<<j; }
  for(int a=0;a<k;a++){ int i=pi[sel[a]], j=pj[sel[a]]; out[i]&=~(1u<<j); out[j]|=1u<<i; }
  tested++;
  /* strong: vertex 0 reaches all and all reach 0 */
  unsigned FULL=(1u<<n)-1; unsigned seen=1, fr=1;
  while(fr){ unsigned nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=out[v]; nf&=~seen; seen|=nf; fr=nf; }
  if(seen!=FULL) return;
  unsigned in[MAXN]={0}; for(int v=0;v<n;v++) for(int w=0;w<n;w++) if(out[v]>>w&1) in[w]|=1u<<v;
  seen=1; fr=1; while(fr){ unsigned nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=in[v]; nf&=~seen; seen|=nf; fr=nf; }
  if(seen!=FULL) return;
  strongc++;
  int c3=n*(n-1)*(n-2)/6; for(int v=0;v<n;v++){ int s=__builtin_popcount(out[v]); c3-=s*(s-1)/2; }
  if(c3<0||c3>=400||!c3ok[c3]) return;
  c3pass++;
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
  for(int t=0;t<ntg;t++){ if(tc3[t]!=c3) continue; int ok=1; for(int j=0;j<n;j++) if(tg[t][j]!=F[j]){ ok=0; break; }
    if(ok){ matches++; printf("MATCH target %d c3 %d arcs",t,c3); for(int a=0;a<k;a++) printf(" %d>%d",pj[sel[a]],pi[sel[a]]); printf(" rows"); for(int v=0;v<n;v++) printf(" %u",out[v]); printf(" F"); for(int j=0;j<n;j++) printf(" %llu",F[j]); putchar(10); fflush(stdout); } }
}
static void rec(int start,int depth,int k){
  if(depth==k){ test(k); return; }
  for(int p=start;p<npairs;p++){ sel[depth]=p; rec(p+1,depth+1,k); }
}
int main(int argc,char**argv){
  n=atoi(argv[1]); K=atoi(argv[2]); FILE*f=fopen(argv[3],"r");
  hsz=1<<16; hset=calloc(hsz,sizeof(u64));
  char line[4096];
  while(fgets(line,sizeof line,f)){ char*p=line; tc3[ntg]=strtol(p,&p,10); for(int j=0;j<n;j++) tg[ntg][j]=strtoull(p,&p,10); c3ok[tc3[ntg]]=1; addH(tg[ntg][0]); ntg++; }
  npairs=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ pi[npairs]=i; pj[npairs]=j; npairs++; }
  for(int k=1;k<=K;k++){ rec(0,0,k); fprintf(stderr,"k=%d done: tested %lld strong %lld c3pass %lld Hpass %lld matches %lld\n",k,tested,strongc,c3pass,hpass,matches); }
  return 0;
}
