#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using I = int64_t;
void need(bool b,const char* s){if(!b)throw runtime_error(s);}
template<class T> void arr(const vector<T>& v){cout<<'[';for(size_t i=0;i<v.size();i++){if(i)cout<<',';cout<<v[i];}cout<<']';}
int main(int argc,char**argv){
 need(argc==3,"usage census <inherited output directory> <max n>");
 string dir=argv[1];int maxn=stoi(argv[2]);need(maxn>=1&&maxn<=8,"order 1..8");
 for(int n=1;n<=maxn;n++){
  ifstream f(dir+"/nodes_"+to_string(n)+".tsv");need(bool(f),"missing nodes input");
  string line;getline(f,line);int rep,tc,merged,sc,classes=0;
  while(f>>rep>>tc>>merged>>sc){
   classes++;int N=1<<n,full=N-1;vector<int> rows(n),deg(n);int bit=0;
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
    bool forward=j==i+1 || ((rep>>bit)&1);if(j>i+1)bit++;
    if(forward)rows[i]|=1<<j;else rows[j]|=1<<i;
   }
   for(int i=0;i<n;i++)deg[i]=__builtin_popcount((unsigned)rows[i]);
   vector<array<I,8>> ends(N),starts(N);
   for(int v=0;v<n;v++)ends[1<<v][v]=starts[1<<v][v]=1;
   for(int S=1;S<N;S++)if((S&(S-1))!=0)for(int v=0;v<n;v++)if(S>>v&1){
    int R=S^(1<<v);for(int u=0;u<n;u++)if(R>>u&1){
     if(rows[u]>>v&1)ends[S][v]+=ends[R][u];
     if(rows[v]>>u&1)starts[S][v]+=starts[R][u];
    }
   }
   I H=accumulate(ends[full].begin(),ends[full].begin()+n,I(0));
   need(H==accumulate(starts[full].begin(),starts[full].begin()+n,I(0)),"start/end checksum");
   array<array<I,8>,8> cover{},ep{};
   for(int S=1;S<full;S++)for(int i=0;i<n;i++)if(S>>i&1)for(int j=0;j<n;j++)if(!(S>>j&1))cover[i][j]+=ends[S][i]*starts[full^S][j];
   function<void(int,int,int)> paths=[&](int first,int last,int visited){
    if(visited==full){ep[first][last]++;return;}
    int avail=rows[last]&~visited;while(avail){int v=__builtin_ctz((unsigned)avail);avail&=avail-1;paths(first,v,visited|(1<<v));}
   };
   for(int v=0;v<n;v++)paths(v,v,1<<v);
   I directH=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)directH+=ep[i][j];need(directH==H,"independent literal paths");
   vector<array<int,8>> autos;array<int,8> p{};
   function<void(int,int)> autom=[&](int depth,int used){
    if(depth==n){autos.push_back(p);return;}
    for(int v=0;v<n;v++)if(!(used>>v&1)&&deg[v]==deg[depth]){
     bool ok=true;for(int u=0;u<depth;u++)if(((rows[u]>>depth)&1)!=((rows[p[u]]>>v)&1)){ok=false;break;}
     if(ok){p[depth]=v;autom(depth+1,used|(1<<v));}
    }
   };autom(0,0);I aut=autos.size();need(aut>0&&aut%2==1&&H==tc*aut,"H/Aut inherited fibre");
   vector<int> subsetOrbit(N);iota(subsetOrbit.begin(),subsetOrbit.end(),0);
   auto imageMask=[&](int S,const array<int,8>& g){int R=0;for(int v=0;v<n;v++)if(S>>v&1)R|=1<<g[v];return R;};
   for(int S=0;S<N;S++)for(auto &g:autos)subsetOrbit[S]=min(subsetOrbit[S],imageMask(S,g));
   vector<I> ext(N),U(N),W(N);vector<array<I,3>> seeds;vector<I> deletion;
   for(int v=0;v<n;v++){int R=full^(1<<v);I h=n==1?1:accumulate(ends[R].begin(),ends[R].begin()+n,I(0));deletion.push_back(h);}
   for(int S=0;S<N;S++){
    I endSelected=0;
    for(int j=0;j<n;j++)if(S>>j&1){
     U[S]+=starts[full][j];endSelected+=ends[full][j];
     for(int i=0;i<n;i++){U[S]+=cover[i][j];if(S>>i&1)W[S]+=cover[i][j];}
    }
    ext[S]=H+U[S]-endSelected-W[S];need(ext[S]>0&&ext[S]%2==1,"odd child count");
    if(W[S]==H-1&&U[S]>0){
     I c=U[S];int s=1;while(c%2==0){c/=2;s++;}seeds.push_back({S,c,s});
     int stab=0;for(auto&g:autos)stab+=imageMask(S,g)==S;need(stab==1,"unit seed stabilizer");
    }
   }
   need(seeds.size()%aut==0,"unit seed orbits");
   int c3=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int k=j+1;k<n;k++){
    int a=(rows[i]>>j)&1,b=(rows[j]>>k)&1,c=(rows[k]>>i)&1;c3+=a==b&&b==c;
   }
   int reach=1;for(int z=0;z<n;z++){int old=reach;for(int v=0;v<n;v++)if(old>>v&1)reach|=rows[v];}
   int back=1;for(int z=0;z<n;z++){int old=back;for(int v=0;v<n;v++)if(rows[v]&old)back|=1<<v;}
   bool strong=reach==full&&back==full;
   cout<<"{\"n\":"<<n<<",\"rep\":"<<rep<<",\"merged\":"<<merged<<",\"sc\":"<<sc<<",\"rows\":";arr(rows);
   cout<<",\"H\":"<<H<<",\"aut\":"<<aut<<",\"tc\":"<<tc<<",\"c3\":"<<c3<<",\"strong\":"<<int(strong)<<",\"delete_H\":";arr(deletion);
   cout<<",\"endpoint\":[";for(int i=0;i<n;i++){if(i)cout<<',';arr(vector<I>(ep[i].begin(),ep[i].begin()+n));}cout<<']';
   cout<<",\"extension\":";arr(ext);cout<<",\"attachment_orbit\":";arr(subsetOrbit);
   cout<<",\"seeds\":[";for(size_t i=0;i<seeds.size();i++){if(i)cout<<',';arr(vector<I>(seeds[i].begin(),seeds[i].end()));}cout<<"]}\n";
  }
  cerr<<"order "<<n<<": "<<classes<<" ordinary classes checked\n";
 }
}
