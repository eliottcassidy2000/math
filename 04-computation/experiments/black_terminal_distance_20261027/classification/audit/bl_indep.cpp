// Independent check of Lemma BL (PROOF.md Sec. 9).  For every strong S on a vertices and z in S:
// if Stab(z) is transitive on Ham(S -> z), the last strong component R_q of S - z has |R_q| >= 2 and (H2),
// then the conclusion must hold.  Counts (S,z) up to Aut(S)-orbits of z.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;
int a; int adj[16][16];
vector<vector<int>> auts; int perm[16]; bool usedv[16];
void autrec(int k){ if(k==a){auts.push_back(vector<int>(perm,perm+a));return;}
  for(int c=0;c<a;c++){ if(usedv[c]) continue; bool ok=true; for(int j=0;j<k&&ok;j++) if(adj[k][j]!=adj[c][perm[j]]) ok=false;
    if(!ok) continue; perm[k]=c; usedv[c]=true; autrec(k+1); usedv[c]=false; } }
int main(int argc,char**argv){
  a=atoi(argv[1]); string line; long long nt=0, hits=0, bad=0, transOnly=0, transH2fail=0, transRq1=0;
  while(getline(cin,line)){
    if(line.empty()||(line[0]!='0'&&line[0]!='1')) continue; nt++;
    int k=0; memset(adj,0,sizeof adj);
    for(int i=0;i<a;i++)for(int j=i+1;j<a;j++){ if(line[k]=='1') adj[i][j]=1; else adj[j][i]=1; k++; }
    auts.clear(); memset(usedv,0,sizeof usedv); autrec(0);
    int full=(1<<a)-1;
    vector<vector<long long>> dp(1<<a, vector<long long>(a,0));
    for(int v=0;v<a;v++) dp[1<<v][v]=1;
    for(int mask=1;mask<=full;mask++) for(int v=0;v<a;v++){ long long c=dp[mask][v]; if(!c) continue;
      for(int w=0;w<a;w++) if(!(mask>>w&1) && adj[v][w]) dp[mask|1<<w][w]+=c; }
    vector<int> seen(a,0);
    for(int z=0;z<a;z++){
      if(seen[z]) continue;
      for(auto&p:auts) seen[p[z]]=1; // orbit representative
      long long stab=0; for(auto&p:auts) if(p[z]==z) stab++;
      long long h=dp[full][z];
      if(h!=stab) continue;
      transOnly++;
      // condensation of S - z: reachability
      int rest=full&~(1<<z);
      vector<int> reach(a,0);
      for(int v=0;v<a;v++){ if(v==z) continue; int r=1<<v, fr=r; while(fr){int x=__builtin_ctz(fr); fr&=fr-1; for(int y=0;y<a;y++) if((rest>>y&1)&&!(r>>y&1)&&adj[x][y]){r|=1<<y; fr|=1<<y;}} reach[v]=r; }
      // components: comp(v) = {w: v reaches w and w reaches v}
      vector<int> comp(a,0);
      for(int v=0;v<a;v++){ if(v==z) continue; int c=0; for(int w=0;w<a;w++) if(w!=z && (reach[v]>>w&1) && (reach[w]>>v&1)) c|=1<<w; comp[v]=c; }
      // order components: R_1 = component reaching everything; last = reached by everything
      vector<int> comps; for(int v=0;v<a;v++) if(v!=z && find(comps.begin(),comps.end(),comp[v])==comps.end()) comps.push_back(comp[v]);
      // sort by size of reach (larger reach = earlier)
      sort(comps.begin(),comps.end(),[&](int c1,int c2){ return __builtin_popcount(reach[__builtin_ctz(c1)])>__builtin_popcount(reach[__builtin_ctz(c2)]); });
      int q=comps.size(); int Rq=comps[q-1];
      int Aq=0; for(int v=0;v<a;v++) if((Rq>>v&1)&&adj[z][v]) Aq|=1<<v; int Bq=Rq&~Aq;
      if(__builtin_popcount(Rq)<2){ transRq1++; continue; }
      bool H2=true;
      for(int W=Rq; W>0 && H2; W=(W-1)&Rq){ int f=0; for(int x=0;x<a;x++) if(W>>x&1) for(int y=0;y<a;y++) if((Rq>>y&1)&&!(W>>y&1)&&adj[x][y]) f++; if(f+__builtin_popcount(W&Bq)<2) H2=false; }
      if(!H2){ transH2fail++; continue; }
      hits++;
      // conclusion
      bool ok=(q>=2);
      auto isC3=[&](int c){ if(__builtin_popcount(c)!=3) return false; return true; }; // strong 3-set is C3
      for(int i=0;i<q && ok;i++){ int c=comps[i]; int sz=__builtin_popcount(c); if(!(sz==1||isC3(c))) ok=false;
        int Ai=0; for(int v=0;v<a;v++) if((c>>v&1)&&adj[z][v]) Ai|=1<<v;
        if(i<q-1){ if(!(Ai==0||Ai==c)) ok=false; if(i==0&&Ai!=c) ok=false; }
        else { if(Ai!=0) ok=false; if(sz!=3) ok=false; } }
      if(!ok){ bad++; printf("BL FAILS: %s z=%d\n",line.c_str(),z); }
    }
  }
  printf("a=%d strong=%lld  transitive (S,z) orbits=%lld  [|R_q|=1: %lld, H2 fails: %lld]  BL-hyp=%lld  violations=%lld\n",a,nt,transOnly,transRq1,transH2fail,hits,bad);
}
