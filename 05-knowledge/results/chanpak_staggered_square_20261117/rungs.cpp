#include <array>
#include <algorithm>
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

struct Arm{vector<U> seq;vector<long long> h;};
unordered_map<U,pair<int,int>> cache;
auto rd=[&](U x){x&=(1ULL<<36)-1;auto it=cache.find(x);if(it!=cache.end())return it->second;return cache[x]=read(x);};
auto H=[&](U x){auto [a,b]=rd(x);return 47680LL*a+(x>>36&1?389:383)*b;};
vector<Arm> bank[2];vector<int> order0={0,1,2,3,4};
do{
 for(int side=0;side<2;side++){
  U x=start^(U(side)<<36);
  Arm a={{x},{H(x)}};
  bool inc=true;
  for(int k:order0){x^=1ULL<<edit[k];a.seq.push_back(x);a.h.push_back(H(x));if(a.h.back()<=a.h[a.h.size()-2])inc=false;}
  if(inc)bank[side].push_back(a);
 }
}while(next_permutation(order0.begin(),order0.end()));
unordered_map<U,vector<U>> rungs;
auto key=[&](U a,U b){int x=0,y=0;for(int j=0;j<5;j++){x|=int((a^start)>>edit[j]&1)<<j;y|=int((b^start)>>edit[j]&1)<<j;}return x+32*y;};
auto try_rung=[&](U a,U b)->vector<U>{
 if((a^b)==(1ULL<<36))return {a,b};
 auto lo=H(a),hi=H(b);
 unordered_map<U,U> pre[2];vector<U> f[2]={{a},{b}};
 pre[0][a]=a;pre[1][b]=b;U m=0; bool found=false;
 for(int i=0;i<24&&!found;i++){
  int d=f[0].size()>f[1].size();if(f[d].empty()||pre[0].size()+pre[1].size()>25000)break;
  vector<U> g;for(U x:f[d]){
   long long h=H(x);
   for(int k=0;k<37;k++){
    U y=x^(1ULL<<k);if(pre[d].count(y))continue;
    long long z=H(y);if(z<lo||z>hi||(!d&&z<=h)||(d&&z>=h))continue;
    pre[d][y]=x;g.push_back(y);
    if(pre[1-d].count(y)){m=y;found=true;break;}
   }
   if(found)break;
  }f[d]=move(g);
 }
 if(!found){if(!f[0].empty()&&!f[1].empty())throw 1;return {};}
 vector<U> p;for(U x=m;;x=pre[0][x]){p.push_back(x);if(x==a)break;}reverse(p.begin(),p.end());
 for(U x=m;x!=b;){x=pre[1][x];p.push_back(x);}return p;
};
int aligned=0;
for(auto &a:bank[0])for(auto &b:bank[1]){
 bool good=true;
 for(int i=0;i<6;i++)if(a.h[i]>=b.h[i])good=false;
 if(!good)continue;
 aligned++;
 for(int i=1;i<5&&good;i++){
  U x=a.seq[i],y=b.seq[i];int k=key(x,y);
  if(!rungs.count(k))rungs[k]=try_rung(x,y);
  if(rungs[k].empty())good=false;
 }
 if(good){
  cout<<"solved "<<aligned<<"\n";
  for(int i=0;i<6;i++){cout<<i<<": ";for(U x:try_rung(a.seq[i],b.seq[i]))cout<<x<<"@"<<H(x)<<" ";cout<<"\n";}
  return 1;
 }
}
cout<<"aligned "<<aligned<<" tested rungs "<<rungs.size()<<" cache "<<cache.size()<<"\n";
if(aligned!=322 || rungs.size()!=5)return 1;
for(auto &r:rungs)if(!r.second.empty()){cout<<"rung "<<r.first;for(U x:r.second)cout<<" "<<x<<"@"<<H(x);cout<<"\n";}
}
