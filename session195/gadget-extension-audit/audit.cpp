#include <array>
#include <algorithm>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_set>
#include <vector>
using Key=std::array<long long,8>;
using State=std::array<int,4>;
struct Hash { size_t operator()(State const& x) const {size_t z=1469598103934665603ULL;for(auto v:x){z^=v;z*=1099511628211ULL;}return z;} };
std::vector<Key> profiles(std::vector<int>const&a){
 int n=a.size(),lim=1<<n;std::vector<std::array<long long,6>>end(lim),start(lim);std::vector<long long>H(lim);H[0]=1;
 for(int i=0;i<n;i++)end[1<<i][i]=start[1<<i][i]=1;
 for(int S=1;S<lim;S++){
  for(int i=0;i<n;i++)if(S>>i&1){
   H[S]+=end[S][i];
   for(int j=0;j<n;j++)if(!(S>>j&1)){
    if(a[i]>>j&1)end[S|(1<<j)][j]+=end[S][i];
    if(a[j]>>i&1)start[S|(1<<j)][j]+=start[S][i];
   }
  }
 }
 std::vector<Key>out;int full=lim-1;
 for(int r=0;r<n;r++){
  long long L=0,R=0,B=0;int del=full^(1<<r);
  for(int S=0;S<lim;S++){
   L+=end[S][r]*H[full^S];R+=start[S][r]*H[full^S];
   if((S&del)==S)B+=H[S]*H[del^S];
  }
  out.push_back({H[full],start[full][r],L,end[full][r],R,B,n==1,H[del]});
 }
 return out;
}
std::vector<Key> inventory(int n){
 std::set<Key>bank;int edges=n*(n-1)/2;
 for(int mask=0;mask<(1<<edges);mask++){
  int k=0;std::vector<int>a(n);
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(mask>>k&1)a[i]|=1<<j;else a[j]|=1<<i;
  for(auto y:profiles(a))bank.insert(y);
 }
 return {bank.begin(),bank.end()};
}
struct Result{size_t states;std::vector<int>missing;};
Result closure(std::vector<Key>const&bank,int X){
 std::deque<State>q;std::unordered_set<State,Hash>seen;std::vector<bool>have(X+1);State first={1,1,1,1};seen.insert(first);q.push_back(first);have[1]=true;
 while(!q.empty()){
  State x=q.front();q.pop_front();long long h=x[0],e=x[1],r=x[2],d=x[3];
  auto take=[&](long long H,long long E,long long R,long long D){
   if(H<h || E<0 || E>D || D<1 || D>H || R<H)throw std::runtime_error("physical inequality");
   if(H>X)return;
   have[H]=true;
   if(R>X)return;
   State next={(int)H,(int)E,(int)R,(int)D};if(seen.insert(next).second)q.push_back(next);
  };
  take(h+r-e,d,2*r-d,d);take(h,0,2*r,d);
  for(auto y:bank)take(y[0]*h-y[1]*e+y[2]*r,y[3]*h-y[6]*e,y[4]*h+y[5]*r,y[7]*h);
 }
 Result out{seen.size(),{}};for(int v=1;v<=X;v+=2)if(!have[v])out.missing.push_back(v);return out;
}
void print(std::ostream&f,Key const&k){for(int j=0;j<8;j++){if(j)f<<' ';f<<k[j];}f<<'\n';}
void report(char const*tag,int X,Result const&r){std::cout<<tag<<' '<<X<<' '<<r.states;for(auto v:r.missing)std::cout<<' '<<v;std::cout<<'\n';}
int main(int argc,char**argv){
 int X=argc>1?std::stoi(argv[1]):20000;std::set<Key>baseSet;
 for(int n=1;n<=5;n++)for(auto k:inventory(n))baseSet.insert(k);
 std::vector<Key>base(baseSet.begin(),baseSet.end()),six=inventory(6);if(base.size()!=65||six.size()!=280)throw std::runtime_error("inventory cardinality");
 std::ofstream f("independent-base.txt");for(auto k:base)print(f,k);f.close();f.open("independent-six.txt");for(auto k:six)print(f,k);f.close();
 auto old=closure(base,173);report("BASE",173,old);
 f.open("independent-single.txt");int c67=0,c173=0;
 for(size_t j=0;j<six.size();j++){
  auto b=base;b.push_back(six[j]);auto r=closure(b,173);
  bool a=std::find(r.missing.begin(),r.missing.end(),67)==r.missing.end(),c=std::find(r.missing.begin(),r.missing.end(),173)==r.missing.end();
  if(a&&c)throw std::runtime_error("single suffices");
  c67+=a;c173+=c;
  f<<j<<' '<<r.states<<' '<<a<<' '<<c<<'\n';
 }
 f.close();std::cout<<"SINGLE "<<six.size()<<' '<<c67<<' '<<c173<<'\n';
 Key ga={29,5,43,4,50,84,0,9},gb={23,1,75,8,24,84,0,9};
 if(std::find(six.begin(),six.end(),ga)==six.end()||std::find(six.begin(),six.end(),gb)==six.end())throw std::runtime_error("missing physical matrix");
 base.push_back(ga);base.push_back(gb);report("EXTENDED",1000,closure(base,1000));if(X!=1000)report("EXTENDED",X,closure(base,X));
}
