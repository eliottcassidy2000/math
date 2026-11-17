#include <array>
#include <vector>
#include <unordered_map>
#include <queue>
#include <iostream>
using namespace std;
using U=unsigned long long;
vector<pair<int,int>> pairs;
pair<int,int> read(U v) {
 int a[9]={},d[512][9]={},s[512][9]={},h[512]={};
 for(int k=0;k<36;k++){auto [i,j]=pairs[k];if(v>>k&1)swap(i,j);a[i]|=1<<j;}
 int reached=256;for(int k=0;k<9;k++)for(int i=0;i<9;i++)if(reached>>i&1)reached|=a[i];if(reached!=511)return {-1,-1};
 for(int i=0;i<9;i++)d[1<<i][i]=1;
 s[256][8]=1;h[0]=1;
 for(int m=1;m<512;m++)for(int i=0;i<9;i++)if((m>>i&1)&&(d[m][i]||s[m][i])){
   h[m]+=d[m][i];
   for(int z=a[i]&~m;z;z&=z-1){int j=__builtin_ctz(z),u=m|(1<<j);d[u][j]+=d[m][i];s[u][j]+=s[m][i];}
 }
 int W=0;
 for(int m=257;m<512;m++){int z=0;for(int j=0;j<9;j++)z+=s[m][j];W+=z*h[511^m];}
 return {W,h[511]-d[511][8]+h[255]-W};
}
int main(){
 for(int i=0;i<9;i++)for(int j=i+1;j<9;j++)pairs.emplace_back(i,j);
 int order[]={4,5,7,6,0,1,2,3,8},pos[9];
 for(int i=0;i<9;i++)pos[order[i]]=i;
 U start=0,end=0,locked=0;
 int goal[]={270,140,264,272,167,463,31,333,82};
 vector<int> edit,extra;
 for(int k=0;k<36;k++){
  auto [i,j]=pairs[k];
  if((pos[i]>pos[j]) ^ (i==4&&j==8))start|=1ULL<<k;
  if(!(goal[i]>>j&1))end|=1ULL<<k;
  if(abs(pos[i]-pos[j])==1 || (i==4&&j==8))locked|=1ULL<<k;
  if((start^end)>>k&1)edit.push_back(k);else if(!(locked>>k&1))extra.push_back(k);
 }
 cerr<<edit.size()<<" "<<extra.size()<<"\n";
 unordered_map<U,pair<int,int>> cache;
 auto rd=[&](U x){auto it=cache.find(x);if(it!=cache.end())return it->second;return cache[x]=read(x);};

 unordered_map<U,U> pre[2];
 vector<U> f[2]={{start},{end}};
 pre[0][start]=start;pre[1][end]=end;
 U hit=0;bool found=false;
 for(int t=0;t<60&&!found;t++){
  int d=f[0].size()>f[1].size();
  if(f[d].empty())break;
  vector<U> g;
  for(U x:f[d]){
   auto [a,b]=rd(x);
   for(int k=0;k<36;k++){
    U y=x^1ULL<<k;if(pre[d].count(y))continue;
    auto [c,e]=rd(y);if(47680LL*c+383LL*e<47680LL*128+383||47680LL*c+389LL*e>47680LL*896+389LL*7||e<=0)continue;
    int z=d?-1:1;
    if(z*(47680LL*(c-a)+383LL*(e-b))<=0 || z*(47680LL*(c-a)+389LL*(e-b))<=0)continue;
    pre[d][y]=x;g.push_back(y);
    if(pre[1-d].count(y)){hit=y;found=true;break;}
   }
   if(found)break;
  }
  cerr<<t<<" "<<d<<" "<<g.size()<<" "<<cache.size()<<"\n";f[d]=move(g);
 }
 if(!found){if(!f[0].empty()&&!f[1].empty()){cerr<<"depth cap; not exhaustive\n";return 1;}cout<<"exhausted: no route; cached reads "<<cache.size()<<"; forward reachable "<<pre[0].size()<<"; backward reachable "<<pre[1].size()<<"\n";return 0;}
 vector<U> p;
 for(U x=hit;;x=pre[0][x]){p.push_back(x);if(x==start)break;}
 for(int i=p.size()-1;i>=0;i--){auto [a,b]=rd(p[i]);if(i+1<(int)p.size()){auto [x,y]=pairs[__builtin_ctzll(p[i]^p[i+1])];cout<<x<<","<<y<<" ";}cout<<a<<" "<<b<<" "<<p[i]<<"\n";}
 for(U x=hit;x!=end;){U y=pre[1][x];auto [u,v]=pairs[__builtin_ctzll(x^y)];auto [a,b]=rd(y);cout<<u<<","<<v<<" "<<a<<" "<<b<<" "<<y<<"\n";x=y;}
}
