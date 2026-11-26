#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <map>
using namespace std;
int n,k,D,U;
array<int,12> H{},deg{};
vector<int> nearsets;
map<uint64_t,uint64_t> shapes;
uint64_t scanned=0,admiss=0,hard=0,crit=0;
array<signed char,4096> memo;
bool connected(){
 int m=1,seen=0;
 while(m){int b=m&-m;m^=b;seen|=b;m|=H[countr_zero(unsigned(b))]&~(seen|m);}
 return seen==U;
}
bool room(){
 for(int s:nearsets){
  unsigned z=s;int e=0;
  while(z){int i=countr_zero(z);z&=z-1;e+=popcount(z&H[i]);if(e>=2)break;}
  if(e<2)return false;
 }
 return true;
}
int cover(int x){
 if(memo[x]>=0)return memo[x];
 const int v=countr_zero(unsigned(x)),b=1<<v,t=x^b;int best=1+cover(t);
 int h=t&H[v];
 for(int c=h;c;c=(c-1)&h){
  bool cl=true;
  for(unsigned z=c;z;z&=z-1){int w=countr_zero(z);if((c&~(1<<w)&H[w])!=(c&~(1<<w))){cl=false;break;}}
  if(cl)best=min(best,1+cover(t^c));
 }
 return memo[x]=best;
}
void test(){
 scanned++;
 if(!connected()||!room())return;
 admiss++;
 memo.fill(-1);memo[0]=0;
 if(cover(U)!=k)return;
 hard++;
 for(int v=0;v<n;v++)if(cover(U^(1<<v))>=k)return;
 crit++;
 vector<int> v(n),d(n);
 for(int i=0;i<n;i++)d[i]=deg[i];
 sort(d.begin(),d.end());
 uint64_t best=~uint64_t(0);
 auto pm=[&](auto&& pm,int i,int used)->void{
  if(i==n){
   uint64_t s=0;
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)s=(s<<1)|(H[v[i]]>>v[j]&1);
   best=min(best,s);return;
  }
  for(int j=0;j<n;j++)if(deg[j]==d[i]&&!(used>>j&1)){
   v[i]=j;pm(pm,i+1,used|(1<<j));
  }
 };
 pm(pm,0,0);
 shapes[best]++;

}
int lim;
void rows(int v);
void picks(int v,int j,int remain){
 if(j==n||!remain){rows(v+1);return;}
 picks(v,j+1,remain);
 if(deg[j]<lim){
  H[v]|=1<<j;H[j]|=1<<v;deg[j]++;deg[v]++;
  picks(v,j+1,remain-1);
  H[v]^=1<<j;H[j]^=1<<v;deg[j]--;deg[v]--;
 }
}
void rows(int v){if(v==n-1)test();else picks(v,v+1,lim-deg[v]);}
int main(int argc,char**argv){
 n=stoi(argv[1]);k=stoi(argv[2]);D=n-k;U=(1<<n)-1;
 for(int s=0;s<=U;s++)if(popcount(unsigned(s))==k)nearsets.push_back(s);
 for(lim=0;lim<=D;lim++){
  for(int j=1;j<=lim;j++){H[0]|=1<<j;H[j]=1;deg[j]=1;}
  deg[0]=lim;
  rows(1);
  H.fill(0);deg.fill(0);
 }
 for(auto [a,b]:shapes)cout<<a<<" "<<b<<"\n";
 cout<<"n,k "<<n<<" "<<k<<" candidates "<<scanned<<" connected and no near clique "<<admiss<<" exact chromatic "<<hard<<" vertex-critical "<<crit<<"\n";
}
