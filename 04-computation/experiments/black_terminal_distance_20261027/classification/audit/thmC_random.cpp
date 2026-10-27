// Random test of Theorem C (PROOF.md Sec.5) at n = 10..14: random self-converse strong nonregular T with
// lambda = delta must have a BLACK Ham path starting at a min-score vertex or ending at a max-score vertex.
// Independent code (no tour.hpp).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <random>
#include <algorithm>
using namespace std;
int n; unsigned out[20]; int sc[20];
bool arc(int u,int v){return out[u]>>v&1;}
int fcut(unsigned X){int f=0; for(int x=0;x<n;x++) if(X>>x&1) f+=__builtin_popcount(out[x]&~X&((1u<<n)-1)); return f;}
bool isBlue(const int*Q){ int r[20]; for(int i=0;i<n;i++) r[Q[i]]=Q[n-1-i];
  for(int u=0;u<n;u++) for(int v=0;v<n;v++) if(u!=v && arc(u,v) && !arc(r[v],r[u])) return false; return true; }
// DFS over Ham paths with start in startSet and end in endSet; return true if a black one is found
long long nodes; 
bool dfs(int*Q,int d,unsigned used,unsigned endSet){
  if(++nodes>20000000) return false;
  if(d==n){ if((endSet>>Q[n-1]&1) && !isBlue(Q)) return true; return false; }
  unsigned c=out[Q[d-1]]&~used;
  while(c){int w=__builtin_ctz(c); c&=c-1; Q[d]=w; if(dfs(Q,d+1,used|1u<<w,endSet)) return true;}
  return false;
}
bool strongT(){ unsigned full=(1u<<n)-1; for(unsigned X=1;X<full;X++) if(fcut(X)==0) return false; return true; }
int main(int argc,char**argv){
  int nmin=atoi(argv[1]), nmax=atoi(argv[2]); long long trials=atoll(argv[3]); mt19937_64 rng(12345);
  for(n=nmin;n<=nmax;n++){
    long long tested=0, fails=0, c2=0, skipped=0; int hist[20]={0};
    for(long long tr=0;tr<trials;tr++){
      memset(out,0,sizeof out);
      int sig[20]; for(int i=0;i<n;i++) sig[i]=n-1-i;
      // random bias profile to vary score sequences
      double bias[20]; uniform_real_distribution<double> U(0,1); double spread=U(rng)*0.9;
      for(int i=0;i<n;i++) bias[i]=0.5+spread*(U(rng)-0.5);
      bool done[20][20]; memset(done,0,sizeof done);
      for(int u=0;u<n;u++) for(int v=u+1;v<n;v++){ if(done[u][v]) continue;
        // orientation u->v with prob bias[u]/(bias[u]+bias[v])
        bool uv = U(rng) < bias[u]/(bias[u]+bias[v]+1e-12)*1.0;
        int a=u,b=v; if(!uv) swap(a,b); out[a]|=1u<<b; done[u][v]=done[v][u]=true;
        int sa=sig[a], sb=sig[b]; // need sig(b)->sig(a)
        if(!((sa==u&&sb==v)||(sa==v&&sb==u))){ out[sb]|=1u<<sa; done[min(sa,sb)][max(sa,sb)]=done[max(sa,sb)][min(sa,sb)]=true; }
      }
      // sanity: tournament & sigma anti-aut
      bool okT=true; for(int u=0;u<n;u++) for(int v=u+1;v<n;v++) if(arc(u,v)==arc(v,u)) okT=false;
      for(int u=0;u<n&&okT;u++) for(int v=0;v<n;v++) if(u!=v && arc(u,v)!=arc(sig[v],sig[u])) okT=false;
      if(!okT){printf("gen error\n"); return 1;}
      for(int v=0;v<n;v++) sc[v]=__builtin_popcount(out[v]);
      int dmin=99,dmax=-1; for(int v=0;v<n;v++){dmin=min(dmin,sc[v]); dmax=max(dmax,sc[v]);}
      if(dmin==dmax) continue; if(!strongT()) continue;
      int lam=1<<20; unsigned full=(1u<<n)-1; for(unsigned X=1;X<full;X++) lam=min(lam,fcut(X));
      if(lam!=dmin) continue;
      unsigned mset=0,Mset=0; for(int v=0;v<n;v++){ if(sc[v]==dmin) mset|=1u<<v; if(sc[v]==dmax) Mset|=1u<<v; }
      // is (c2)? all min deletions strong
      bool allStrong=true; for(int v=0;v<n;v++) if(mset>>v&1){ // strong T-v?
        unsigned F=full&~(1u<<v); bool st=true; for(unsigned X=F; X; X=(X-1)&F){ if(X==F) continue; int f=0; for(int x=0;x<n;x++) if(X>>x&1) f+=__builtin_popcount(out[x]&F&~X); if(f==0){st=false;break;} } if(!st){allStrong=false;break;} }
      if(allStrong) c2++;
      hist[__builtin_popcount(mset)]++;
      tested++; bool found=false; int Q[20];
      for(int v=0;v<n&&!found;v++) if(mset>>v&1){ Q[0]=v; nodes=0; if(dfs(Q,1,1u<<v,full)) found=true; }
      if(!found){ // paths ending in M: search all starts
        for(int v=0;v<n&&!found;v++){ Q[0]=v; nodes=0; if(dfs(Q,1,1u<<v,Mset)) found=true; } }
      if(!found){ fails++; printf("n=%d NO singleton-type black descent found (maybe node cap) scores:",n); for(int v=0;v<n;v++) printf(" %d",sc[v]); printf("\n"); }
    }
    printf("n=%d tested(lambda=delta, nonregular SC strong)=%lld  (c2 subcase=%lld)  fails=%lld  mu-hist:",n,tested,c2,fails);
    for(int i=1;i<=n/2;i++) printf(" %d:%d",i,hist[i]); printf("\n");
  }
}
