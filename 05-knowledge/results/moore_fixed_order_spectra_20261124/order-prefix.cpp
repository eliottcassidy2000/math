#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>
using namespace std;
using U=uint64_t;
struct Row {int h,e,R,root;U mask;};
using B=map<array<int,3>,Row>;
B bank[11];
void add(int n,U mask){
 int o[10]{},end[1024][10],start[1024][10],hs[1024];int q=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,q++)
  if(mask>>q&1)o[j]|=1<<i;else o[i]|=1<<j;
 int full=(1<<n)-1;hs[0]=1;
 for(int S=1;S<=full;S++){
  hs[S]=0;
  for(int b=S;b;b&=b-1){
   int i=countr_zero((unsigned)b),T=S^(1<<i),a=!T,z=!T;
   for(int g=T&o[i];g;g&=g-1)z+=start[T][countr_zero((unsigned)g)];
   for(int g=T&~o[i];g;g&=g-1)a+=end[T][countr_zero((unsigned)g)];
   hs[S]+=end[S][i]=a;start[S][i]=z;
  }
 }
 for(int i=0;i<n;i++){
  int r=0,T=full^(1<<i);
  for(int g=T;;g=(g-1)&T){
   r+=hs[g]*start[full^g][i];
   if(!g)break;
  }
  array<int,3> a={hs[full],end[full][i],r};
  bank[n].try_emplace(a,Row{a[0],a[1],a[2],i,mask});
 }
}
int main(int argc,char**argv){
 const int mx=argc>1?atoi(argv[1]):20, lim=1<<23;
 const int rounds=argc>2?atoi(argv[2]):1;
 assert(mx>=16&&mx<=21&&rounds>=1&&rounds<=16);
 add(1,0);
 for(int n=2;n<=8;n++){
  ifstream in("bank-"+to_string(n));string s;int nclass=0;
  while(in>>s){U m=0;for(size_t i=0;i<s.size();i++)m|=U(s[i]=='0')<<i;add(n,m);nclass++;}
  assert((nclass==array<int,9>{0,0,1,2,4,12,56,456,6880}[n]));
 }
 for(int n:{9,10}){
  uint32_t st=1325623;auto next=[&](){st=1664525u*st+1013904223u;return st;};
  int bits=n*(n-1)/2;
  for(int i=0;i<(n==9?8192:4096)*rounds;i++){
   U low=next(),high=next();add(n,(low|(high<<32))&((U(1)<<bits)-1));
  }
 }
 for(int n=1;n<=10;n++)cerr<<"root bank "<<n<<" "<<bank[n].size()<<"\n";
 vector<int> g(1<<mx,1);
 for(unsigned x=1;x<g.size();x++)g[x]=g[x/2]+(x%2?0:x/2+((x/2)&1));
 vector<uint8_t> order(lim/2,255);order[0]=1;
 auto put=[&](int n,long long H){
  assert(H>0&&(H&1));
  if(H<lim&&order[H/2]>n)order[H/2]=n;
 };
 // F1: two one-apex words.
 for(unsigned u=0;u<g.size();u++){
  int remain=mx-2-bit_width(u);
  if(remain<0||(1u<<remain)<=u)continue;
  for(unsigned v=u;v<(1u<<remain);v++)
   put(bit_width(u)+bit_width(v)+2,1LL*(u+1)*(v+1)+1LL*g[u]*g[v]-(!(u&1)&&!(v&1)));
 }
 // F2: with the second side the converse of a bank tournament.
 for(int n=1;n<=10;n++)for(auto &[p,x]:bank[n]){
  put(n,x.h);
  for(unsigned w=0;w<(1u<<(mx-n-1));w++){
   const long long H=1LL*(w+1)*x.R+1LL*g[w]*x.h-(!(w&1))*x.e;
   put(n+1+bit_width(w),H);
  }
 }
 for(int n=12;n<=mx;n++){int H=1;while(H<lim&&(H==7||H==21||order[H>>1]<=n))H+=2;cout<<"F1+F2 "<<n<<" first hole "<<H<<"\n";}
 // Ordinary joins with two finite banks.  Exclude the largest pairs to keep replay small.
 vector<Row> decks[11];
 for(int n=1;n<=10;n++){
  for(auto &p:bank[n])decks[n].push_back(p.second);
  sort(decks[n].begin(),decks[n].end(),[](auto x,auto y){return x.R<y.R;});
 }
 for(int a=1;a<=8;a++)for(int b=a;b<=min(10,mx-a);b++){
  if(a+b>17||a==8&&b==9)continue;
  for(const auto& x:decks[a])for(const auto& y:decks[b]){
   if(1LL*x.R*y.R>=lim)break; // h_x h_y >= e_x e_y
   put(a+b,1LL*x.h*y.h+1LL*x.R*y.R-1LL*x.e*y.e);
  }
 }
 for(int n=12;n<=mx;n++){int H=1;while(H<lim&&(H==7||H==21||order[H>>1]<=n))H+=2;cout<<"union "<<n<<" first hole "<<H<<" ; targets 20479 and 458141 minimum among this bank "<<int(order[20479/2])<<","<<int(order[458141/2])<<"\n";}
}
