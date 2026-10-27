// Independent check of the last-block lemma (PROOF.md Sec. 9), no shared code with tour.hpp.
// For every strong R on m vertices and every A (B = V\A nonempty):
//  H1: G = Aut(R;A) transitive on Ham(R -> B)  <=>  |Ham(R->B)| == |G|  (free action)
//  H2: every nonempty W: e(W -> V\W) + |W cap B| >= 2
//  H2w: |B| >= 2 only (the only instance of H2 used in the written proof)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
using namespace std;
int m; int adj[16][16];
vector<vector<int>> auts;
int perm[16]; bool usedv[16];
void autrec(int k){
  if(k==m){auts.push_back(vector<int>(perm,perm+m));return;}
  for(int c=0;c<m;c++){ if(usedv[c]) continue; bool ok=true;
    for(int j=0;j<k&&ok;j++) if(adj[k][j]!=adj[c][perm[j]]) ok=false;
    if(!ok) continue; perm[k]=c; usedv[c]=true; autrec(k+1); usedv[c]=false; }
}
int main(int argc,char**argv){
  m=atoi(argv[1]); string line; long long nt=0, cntH1H2=0, cntH1w=0, cntH1H2_Anonempty=0, cntH1w_An=0;
  map<string,long long> examplesW;
  while(getline(cin,line)){
    if(line.empty()||(line[0]!='0'&&line[0]!='1')) continue;
    nt++; int k=0; memset(adj,0,sizeof adj);
    for(int i=0;i<m;i++)for(int j=i+1;j<m;j++){ if(line[k]=='1') adj[i][j]=1; else adj[j][i]=1; k++; }
    auts.clear(); memset(usedv,0,sizeof usedv); autrec(0);
    int full=(1<<m)-1;
    // dp count of paths
    vector<vector<long long>> dp(1<<m, vector<long long>(m,0));
    for(int v=0;v<m;v++) dp[1<<v][v]=1;
    for(int mask=1;mask<=full;mask++) for(int v=0;v<m;v++){ long long c=dp[mask][v]; if(!c) continue;
      for(int w=0;w<m;w++) if(!(mask>>w&1) && adj[v][w]) dp[mask|1<<w][w]+=c; }
    vector<long long> hend(m); for(int e=0;e<m;e++) hend[e]=dp[full][e];
    vector<int> outW(1<<m,0);
    for(int W=1;W<=full;W++){ int f=0; for(int a=0;a<m;a++) if(W>>a&1) for(int b=0;b<m;b++) if(!(W>>b&1)&&adj[a][b]) f++; outW[W]=f; }
    for(int A=0;A<full;A++){
      int B=full&~A; long long hB=0; for(int e=0;e<m;e++) if(B>>e&1) hB+=hend[e];
      long long g=0; for(auto&p:auts){ bool ok=true; for(int v=0;v<m&&ok;v++) if(((A>>v)&1)!=((A>>p[v])&1)) ok=false; if(ok) g++; }
      bool H1=(hB==g);
      bool H2=true; for(int W=1;W<=full&&H2;W++) if(outW[W]+__builtin_popcount(W&B)<2) H2=false;
      bool H2w=__builtin_popcount(B)>=2;
      if(H1&&H2){ cntH1H2++; if(A) {cntH1H2_Anonempty++; printf("H1+H2 with A nonempty: %s A=%x\n",line.c_str(),A);} }
      if(H1&&H2w){ cntH1w++; if(A){cntH1w_An++; if(cntH1w_An<=5) printf("H1+|B|>=2 with A nonempty: %s A=%x |G|=%lld hB=%lld\n",line.c_str(),A,g,hB);} }
    }
  }
  printf("m=%d strong=%lld  H1+H2: %lld (A nonempty: %lld)   H1+|B|>=2: %lld (A nonempty: %lld)\n",m,nt,cntH1H2,cntH1H2_Anonempty,cntH1w,cntH1w_An);
}
