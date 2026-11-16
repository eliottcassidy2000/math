#include<bits/stdc++.h>
using namespace std;
struct Counts {
 int n,F; vector<int>P,h;
 Counts(const vector<int>&a):n(a.size()),F((1<<n)-1),P((F+1)*n),h(F+1){
  h[0]=1; for(int S=1;S<=F;S++)for(int b=S;b;b&=b-1){
   int i=__builtin_ctz(b),R=S^(1<<i),v=!R;
   for(int x=R&a[i];x;x&=x-1)v+=P[R*n+__builtin_ctz(x)];
   P[S*n+i]=v;h[S]+=v;
  }
 }
};
struct W{vector<int>a;int q,H,X,E,O,g;};
vector<int> core(vector<int> a){a.resize(9);for(auto&v:a)v&=511;return a;}
int main(int argc,char**argv){
 vector<vector<int>> seeds={{228,29,304,501,801,770,822,374,3,399},{242,156,465,821,640,406,186,8,723,231},{314,32,387,710,526,412,439,19,666,231},{344,849,579,918,548,463,24,855,596,97}};
 for(auto a:seeds){Counts c(a);long long E=0,O=0;for(int S=512;S<=1023;S++)((__builtin_popcount(S)%2)?O:E)+=c.P[S*10+9]*c.h[1023^S];cout<<"seed "<<c.h[1023]<<" "<<E<<" "<<O<<" "<<c.h[511]<<"\n";}
 mt19937 rng(6409);int n=9,F=511;vector<pair<int,int>> arcs;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)arcs.emplace_back(i,j);
 map<pair<int,int>,array<W,3>> hits;set<uint64_t> seen;int T=argc>1?atoi(argv[1]):100000;
 for(int it=0;it<T;it++){
  vector<int>a(n);
  if(it<4)a=core(seeds[it]);
  else if(it%8<6){a=core(seeds[rng()%4]);int edits=1+rng()%8;while(edits--){auto [i,j]=arcs[rng()%36];a[i]^=1<<j;a[j]^=1<<i;}}
  else for(auto [i,j]:arcs){if(rng()&1)a[i]|=1<<j;else a[j]|=1<<i;}
  uint64_t code=0;for(auto[i,j]:arcs)code=2*code+((a[i]>>j)&1);
  if(!seen.insert(code).second)continue;
  Counts c(a);int g=c.h[F];if(g>6019)continue;
  int e[9]={},o[9]={};
  for(int S=1;S<=F;S++)for(int b=S;b;b&=b-1){
   int i=__builtin_ctz(b),v=c.P[S*n+i]*c.h[F^S];
   ((__builtin_popcount(S)%2)?e:o)[i]+=v;
  }
  int E[512]={},O[512]={};O[0]=g;
  for(int q=1;q<=F;q++){
   int i=__builtin_ctz(q),p=q&(q-1);E[q]=E[p]+e[i];O[q]=O[p]+o[i];
  }
  // Tournaments have H = g + 2*sum(even paths from contact to noncontact).
  int ev[9][9]={};
  vector<int> p((F+1)*n);
  for(int st=0;st<n;st++){
   fill(p.begin(),p.end(),0);p[(1<<st)*n+st]=1;
   for(int S=1<<st;S<=F;S++)if(S>>st&1)for(int b=S;b;b&=b-1){
    int j=__builtin_ctz(b),v=p[S*n+j];if(!v)continue;
    if(!(__builtin_popcount(S)%2))ev[st][j]+=v*c.h[F^S];
    for(int A=(F^S)&a[j];A;A&=A-1){int u=__builtin_ctz(A);p[(S|(1<<u))*n+u]+=v;}
   }
  }
  int cut[512]={};
  for(int q=1;q<F;q++){
   int i=__builtin_ctz(q),p=q&(q-1);cut[q]=cut[p]+o[i];for(int b=p;b;b&=b-1){int j=__builtin_ctz(b);cut[q]-=ev[i][j]+ev[j][i];}
   int H=g+2*cut[q];if(H<6015||H>6019)continue;
   int x=E[q]+O[q],odd=x>>__builtin_ctz(x);auto &row=hits[{g,odd}]; int k=(H-6015)/2;
   if(row[k].a.empty()){
    row[k]={a,q,H,x,E[q],O[q],g};
    if(!row[0].a.empty()&&!row[1].a.empty()&&!row[2].a.empty()){
     cout<<"TRIPLE iteration "<<it<<" g "<<g<<" odd "<<odd<<"\n";
     for(auto w:row){cout<<w.H<<" "<<w.E<<" "<<w.O<<" "<<w.g<<" rows";for(int j=0;j<n;j++)cout<<" "<<(w.a[j]|(((w.q>>j&1)^1)<<9));cout<<" "<<w.q<<"\n";}return 0;
    }
   }
  }
 }
 int pairs=0; for(auto&[key,r]:hits){int z=0;for(auto&w:r)z+=!w.a.empty();pairs+=z>=2;}
 cout<<"DONE "<<seen.size()<<" keys "<<hits.size()<<" pairs "<<pairs<<"\n";
}
