/* random search: even graphs on 10 vertices, no universal vertex, d in a target set (cover.py E0 misses at n=10) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int n=10;
static int targets[64]={0};
static unsigned rng=12345u;
static unsigned rnd(void){ rng^=rng<<13; rng^=rng>>17; rng^=rng<<5; return rng; }
static int dcount(int *adj){ int full=(1<<n)-1,c=0; int cl[16]; for(int v=0;v<n;v++) cl[v]=adj[v]|(1<<v);
  for(int S=0;S<(1<<n);S++){ int dom=0,s=S; while(s){int v=__builtin_ctz(s); s&=s-1; dom|=cl[v];} if(dom==full) c++; } return c; }
/* automorphism enumeration with sign check: returns 1 if even */
static int adjg[16], perm[16], used[16], order_[16], deg[16], oddfound;
static void bt(int i){ if(oddfound) return; if(i==n){ int rev=0; for(int u=0;u<n;u++) for(int v=u+1;v<n;v++) if((adjg[u]>>v)&1 && perm[u]>perm[v]) rev++; if(rev&1) oddfound=1; return; }
  int v=order_[i]; for(int w=0;w<n;w++){ if(used[w]||deg[w]!=deg[v]) continue; int ok=1; for(int j=0;j<i;j++){int u=order_[j]; if(((adjg[v]>>u)&1)!=((adjg[w]>>perm[u])&1)){ok=0;break;}} if(ok){perm[v]=w;used[w]=1;bt(i+1);used[w]=0;perm[v]=-1;} } }
static int is_even(int *adj){ memcpy(adjg,adj,sizeof(int)*n); for(int v=0;v<n;v++){deg[v]=__builtin_popcount(adj[v]); order_[v]=v; used[v]=0; perm[v]=-1;}
  for(int a=0;a<n;a++) for(int b=a+1;b<n;b++) if(deg[order_[b]]>deg[order_[a]]){int t=order_[a];order_[a]=order_[b];order_[b]=t;}
  oddfound=0; bt(0); return !oddfound; }
int main(int argc,char**argv){ long iters=atol(argv[1]); n=atoi(argv[2]); rng=(unsigned)atoi(argv[3]); for(int i=4;i<argc;i++) targets[i-4]=atoi(argv[i]); int found[64]={0}; /* was [16]: overflowed with 18/20 targets (audit2 Finding 3) */ int adj[16];
  for(long it=0;it<iters;it++){ memset(adj,0,sizeof(adj)); int p=30+rnd()%45; /* edge prob percent */
    for(int u=0;u<n;u++) for(int v=u+1;v<n;v++) if((int)(rnd()%100)<p){adj[u]|=1<<v;adj[v]|=1<<u;}
    int univ=0; for(int v=0;v<n;v++) if(__builtin_popcount(adj[v])==n-1) univ=1; if(univ) continue;
    int d=dcount(adj); for(int t=0;targets[t];t++) if(d==targets[t] && !found[t]){ if(is_even(adj)){ found[t]=1; printf("FOUND d=%d edges:",d); for(int u=0;u<n;u++) for(int v=u+1;v<n;v++) if((adj[u]>>v)&1) printf(" %d-%d",u,v); printf("\n"); fflush(stdout);} } }
  for(int t=0;targets[t];t++) printf("target %d %s\n",targets[t],found[t]?"realised":"not found"); return 0; }
