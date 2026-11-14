// Order-n census of strong tournaments with no two vertex-disjoint triangles (pairwise-meeting odd cycles).
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
typedef unsigned long long u64;
static int n; static unsigned adj[16];
static int strong(void){ unsigned full=(1u<<n)-1, seen=1, fr=1; while(fr){unsigned nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=adj[v]; fr=nf&~seen; seen|=nf;} if(seen!=full) return 0;
  unsigned radj[16]={0}; for(int i=0;i<n;i++) for(int j=0;j<n;j++) if(adj[i]>>j&1) radj[j]|=1u<<i;
  seen=1; fr=1; while(fr){unsigned nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=radj[v]; fr=nf&~seen; seen|=nf;} return seen==full; }
static u64 H(void){ static u64 dp[1<<16][16]; unsigned full=(1u<<n)-1; for(unsigned m=0;m<=full;m++) for(int v=0;v<n;v++) dp[m][v]=0;
  for(int v=0;v<n;v++) dp[1u<<v][v]=1; for(unsigned m=1;m<=full;m++) for(int v=0;v<n;v++){ u64 c=dp[m][v]; if(!c) continue; unsigned nb=adj[v]&~m; while(nb){int w=__builtin_ctz(nb); nb&=nb-1; dp[m|1u<<w][w]+=c;} }
  u64 h=0; for(int v=0;v<n;v++) h+=dp[full][v]; return h; }
int main(int argc,char**argv){
  char line[256]; long long total=0, strongc=0, omega=0, hub=0, nohub=0; u64 hv[64]; int nhv=0; int taus[8]={0};
  while(fgets(line,sizeof line,stdin)){
    if(line[0]=='>') continue; int len=strlen(line); while(len&&(line[len-1]=='\n')) line[--len]=0;
    n=0; while(n*(n-1)/2<len) n++; memset(adj,0,sizeof adj); int k=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[k]=='1') adj[i]|=1u<<j; else adj[j]|=1u<<i; k++; }
    total++; if(!strong()) continue; strongc++;
    unsigned tri[200]; int nt=0, disjoint=0;
    for(int a=0;a<n&&!disjoint;a++) for(int b=a+1;b<n&&!disjoint;b++) for(int c=b+1;c<n&&!disjoint;c++){
      int cyc=((adj[a]>>b&1)&&(adj[b]>>c&1)&&(adj[c]>>a&1))||((adj[a]>>c&1)&&(adj[c]>>b&1)&&(adj[b]>>a&1));
      if(!cyc) continue; unsigned t=(1u<<a)|(1u<<b)|(1u<<c);
      for(int q=0;q<nt;q++) if(!(tri[q]&t)){disjoint=1;break;} tri[nt++]=t; }
    if(disjoint) continue; omega++;
    int ishub=0; for(int v=0;v<n&&!ishub;v++){ int all=1; for(int q=0;q<nt;q++) if(!(tri[q]>>v&1)){all=0;break;} if(all) ishub=1; }
    if(ishub){hub++; continue;}
    nohub++; int tau=9; for(int u=0;u<n&&tau>2;u++) for(int v=u+1;v<n;v++){ unsigned S=(1u<<u)|(1u<<v); int all=1; for(int q=0;q<nt;q++) if(!(tri[q]&S)){all=0;break;} if(all){tau=2;break;} }
    taus[tau<8?tau:7]++;
    u64 h=H(); int seen=0; for(int q=0;q<nhv;q++) if(hv[q]==h) seen=1; if(!seen&&nhv<64) hv[nhv++]=h;
  }
  printf("order %d: tournaments %lld, strong %lld, no two disjoint triangles %lld, with hub %lld, hub-free %lld (tau=2: %d, tau>2: %d)\nhub-free H values:",n,total,strongc,omega,hub,nohub,taus[2],nohub-taus[2]);
  for(int q=0;q<nhv;q++) printf(" %llu",hv[q]); printf("\n"); return 0; }
