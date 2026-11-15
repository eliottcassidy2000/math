// For unit cuts X of tournaments of order n: compute Z, E, R with delta = E + R - Z, by enumerating Hamiltonian paths
// and 2-path covers explicitly (small n). Print distribution of (sign, Z, E, R) and verify the identity.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
static int n, out[16];
static int path[16]; static i64 Hc;
static i64 Zc, Ec;
static void dfs(int len, int used, int X){
  if(len==n){ Hc++; int steps=0; for(int i=0;i+1<n;i++) if(((X>>path[i])&1) && !((X>>path[i+1])&1)) steps++;
    if(steps==0) Zc++; else Ec+=steps-1; return; }
  int last=path[len-1]; int nb=out[last]&~used;
  while(nb){ int u=__builtin_ctz(nb); nb&=nb-1; path[len]=u; dfs(len+1,used|1<<u,X); }
}
// R: pairs (P1,P2) both nonempty covering V, end(P1) notin X, start(P2) in X, start(P2) -> end(P1)
// count via DP: E1[S][v] = #paths covering S ending at v ; S2[S][w] = #paths covering S starting at w
static i64 dpE[1<<15][15], dpS[1<<15][15];
int main(int argc,char**argv){ n=atoi(argv[1]); char line[4096]; int N=1<<n;
  long long tab[2][4][4][4]={{{{0}}}}; long long bad=0, cnt=0, ecut=0; i64 mind=1LL<<60, mindE=1LL<<60, mind2=1LL<<60, minRZ=1LL<<60;
  long long zbig[2]={0}, rbig[2]={0}, ebig[2]={0};
  long long zhist[2][8]={{0}}, rhist[2][8]={{0}}, ehist[2][8]={{0}};
  while(fgets(line,sizeof line,stdin)){ if(line[0]=='>') continue;
    int k=0; for(int i=0;i<n;i++) out[i]=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[k]=='1') out[i]|=1<<j; else out[j]|=1<<i; k++; }
    // DP
    memset(dpE,0,sizeof(i64)*N*15); memset(dpS,0,sizeof(i64)*N*15);
    int in[16]; for(int v=0;v<n;v++) in[v]=0; for(int v=0;v<n;v++) for(int u=0;u<n;u++) if((out[v]>>u)&1) in[u]|=1<<v;
    for(int v=0;v<n;v++){ dpE[1<<v][v]=1; dpS[1<<v][v]=1; }
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ i64 x=dpE[S][v]; if(x){ int nb=out[v]&~S; while(nb){int u=__builtin_ctz(nb); nb&=nb-1; dpE[S|1<<u][u]+=x;} }
      i64 y=dpS[S][v]; if(y){ int nb=in[v]&~S; while(nb){int u=__builtin_ctz(nb); nb&=nb-1; dpS[S|1<<u][u]+=y;} } }
    i64 H=0; for(int v=0;v<n;v++) H+=dpE[N-1][v];
    for(int X=1;X<N-1;X++){
      // H(T+x) via 2-covers
      i64 tot=0; for(int v=0;v<n;v++){ if((X>>v)&1) tot+=dpS[N-1][v]; else tot+=dpE[N-1][v]; }
      i64 Rr=0;
      for(int S=1;S<N-1;S++){ int C=(N-1)^S; for(int u=0;u<n;u++) if(((S>>u)&1) && !((X>>u)&1) && dpE[S][u]) for(int w=0;w<n;w++) if(((C>>w)&1) && ((X>>w)&1) && dpS[C][w]){
            tot+=dpE[S][u]*dpS[C][w]; if((out[w]>>u)&1) Rr+=dpE[S][u]*dpS[C][w]; } }
      i64 delta=tot-2*H;
      Hc=0; Zc=0; Ec=0; for(int s=0;s<n;s++){ path[0]=s; dfs(1,1<<s,X); }
      cnt++; int sg=delta>0?0:1; if(Ec>0){ ecut++; if(delta<mind) mind=delta; if(delta-Ec<mindE) mindE=delta-Ec; if(delta-2*Ec+1<mind2) mind2=delta-2*Ec+1; if(Rr-Zc<minRZ) minRZ=Rr-Zc; }
      if(Ec+Rr-Zc!=delta) bad++;
      zhist[sg][Zc<7?Zc:7]++; rhist[sg][Rr<7?Rr:7]++; ehist[sg][Ec<7?Ec:7]++;
    }
  }
  printf("order %d: all cuts %lld; cuts with E>=1: %lld; min delta over them %lld; min (delta-E) %lld; min (delta-2E+1) %lld; min (R-Z) %lld\n",n,cnt,ecut,mind,mindE,mind2,minRZ);
  for(int s=0;s<2;s++){ printf(" sign %s: Z hist",s?"-1":"+1"); for(int i=0;i<8;i++) printf(" %lld",zhist[s][i]); printf(" | E hist"); for(int i=0;i<8;i++) printf(" %lld",ehist[s][i]); printf(" | R hist"); for(int i=0;i<8;i++) printf(" %lld",rhist[s][i]); printf("\n"); }
  return 0;}
