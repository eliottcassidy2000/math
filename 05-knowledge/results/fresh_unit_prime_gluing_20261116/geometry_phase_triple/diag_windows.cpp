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
 mt19937 rng(argc>2?atoi(argv[2]):370028);int n=9,F=511;vector<pair<int,int>> arcs;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)arcs.emplace_back(i,j);
 map<pair<int,int>,array<W,3>> hits;set<uint64_t> seen;int T=argc>1?atoi(argv[1]):100000;
 n=atoi(argv[1]);if(n<9||n>11)return 2;F=(1<<n)-1;
 arcs.clear();for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)arcs.emplace_back(i,j);
 vector<vector<pair<int,int>>> orbits;set<pair<int,int>> marked;
 auto sh=[](int i){return i<9?3*(i/3)+(i+1)%3:i;};
 for(auto [i,j]:arcs)if(!marked.count({i,j})){
  vector<pair<int,int>> o;int a=i,b=j;
  do {o.push_back({a,b});marked.emplace(min(a,b),max(a,b));a=sh(a);b=sh(b);}while(a!=i||b!=j);
  orbits.push_back(o);}
 int best=1000000000,events=0;set<int> mids; for(int it=0;it<(1<<orbits.size());it++){
  vector<int>a(n);for(int j=0;j<int(orbits.size());j++)for(auto [u,v]:orbits[j]){if(it>>j&1)swap(u,v);a[u]|=1<<v;}
  uint64_t code=n;for(auto[i,j]:arcs)code=2*code+((a[i]>>j)&1);
  if(!seen.insert(code).second)continue;
  Counts c(a);int g=c.h[F];if(g>6019)continue;
  int e[11]={},o[11]={};
  for(int S=1;S<=F;S++)for(int b=S;b;b&=b-1){
   int i=__builtin_ctz(b),v=c.P[S*n+i]*c.h[F^S];
   ((__builtin_popcount(S)%2)?e:o)[i]+=v;
  }
  int E[2048]={},O[2048]={};O[0]=g;
  for(int q=1;q<=F;q++){
   int i=__builtin_ctz(q),p=q&(q-1);E[q]=E[p]+e[i];O[q]=O[p]+o[i];
  }
  vector<uint64_t> sigs(F-1); for(int q=1;q<F;q++)sigs[q-1]=(uint64_t(E[q])<<32)|O[q];
  sort(sigs.begin(),sigs.end());bool viable=false;for(int j=2;j<int(sigs.size());j++)if(sigs[j]==sigs[j-2]){viable=true;break;}
  if(!viable)continue;
  // Tournaments have H = g + 2*sum(even paths from contact to noncontact).
  int ev[11][11]={};
  vector<int> p((F+1)*n);
  for(int st=0;st<n;st++){
   fill(p.begin(),p.end(),0);p[(1<<st)*n+st]=1;
   for(int S=1<<st;S<=F;S++)if(S>>st&1)for(int b=S;b;b&=b-1){
    int j=__builtin_ctz(b),v=p[S*n+j];if(!v)continue;
    if(!(__builtin_popcount(S)%2))ev[st][j]+=v*c.h[F^S];
    for(int A=(F^S)&a[j];A;A&=A-1){int u=__builtin_ctz(A);p[(S|(1<<u))*n+u]+=v;}
   }
  }
  int cut[2048]={};

  map<pair<int,int>,map<int,int>> same;
  for(int q=1;q<F;q++){
   int i=__builtin_ctz(q),p=q&(q-1);cut[q]=cut[p]+o[i];for(int b=p;b;b&=b-1){int j=__builtin_ctz(b);cut[q]-=ev[i][j]+ev[j][i];}
   int H=g+2*cut[q];same[{E[q],O[q]}][H]=q;
  }
  for(auto&[key,r]:same)for(auto [H,q]:r)if(r.count(H+2)&&r.count(H+4)){
   events++;mids.insert(H+2);
   if(abs(H+2-6017)<best){
    best=abs(H+2-6017);cout<<"BEST "<<H<<" "<<H+2<<" "<<H+4<<" E,O,g "<<key.first<<" "<<key.second<<" "<<g<<" masks "<<q<<" "<<r[H+2]<<" "<<r[H+4]<<" D";for(int x:a)cout<<" "<<x;cout<<"\n";
   }
   if(best==0)return 0;
  }
 }
 if(const char* path=getenv("DUMP")){ofstream out(path);for(auto&[key,r]:hits)for(auto&w:r)if(!w.a.empty()){
  out<<key.first<<" "<<key.second<<" "<<w.H<<" "<<w.E<<" "<<w.O;
  for(int j=0;j<int(w.a.size());j++)out<<" "<<(w.a[j]|(((w.q>>j&1)^1)<<w.a.size()));
  out<<" "<<w.q<<"\n";
 }}
 cout<<"END distinct "<<mids.size()<<" events "<<events<<"\n"; for(auto h:mids)cout<<h<<" ";cout<<"\n"; int pairs=0; for(auto&[key,r]:hits){int z=0;for(auto&w:r)z+=!w.a.empty();pairs+=z>=2;}
 cout<<"DONE "<<seen.size()<<" keys "<<hits.size()<<" pairs "<<pairs<<"\n";
}
