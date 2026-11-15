// Exact mixed-menu audit by orbit-complete vertex extensions.
// Build: g++ -O3 -std=c++17 movable_extensions.cpp -o /tmp/movable-extensions
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <utility>
#include <vector>
using namespace std;
using R = vector<unsigned>;

struct P {int H; vector<int> d; vector<bool> stable; bool strong;};
static P profile(const R& a){
 int n=a.size(), full=(1<<n)-1, dp[1024][10]={}, h[1024]={};
 h[0]=1;
 for(int s=1;s<=full;s++){
   unsigned cur=s;
   while(cur){
     int i=__builtin_ctz(cur); cur&=cur-1;
     unsigned t=s^(1u<<i), fol=t&a[i];
     if(!t)dp[s][i]=1;
     else while(fol){
       int j=__builtin_ctz(fol); fol&=fol-1;
       dp[s][i]+=dp[t][j];
     }
     h[s]+=dp[s][i];
   }
 }
 vector<int> d(n,-h[full]);
 for(int s=1;s<=full;s++) if(__builtin_parity(s)==0)
   for(int i=0;i<n;i++) d[i]+=2*h[full^s]*dp[s][i];
 vector<bool> stable(n,true);bool strong=true;
 for(int i=0;i<n;i++){
   strong &= dp[full][i]>0;
   for(int j=0;j<n;j++) if(i!=j&& !(a[i]>>j&1u) && dp[full^(1u<<i)][j]) stable[i]=false;
 }
 return {h[full],d,stable,strong};
}

static R single(const R& a,int q){
 int n=a.size();
 R b=a;
 for(int i=0;i<n;i++)if(i!=q)b[i]|=1u<<n;
 b.push_back(1u<<q); return b;
}
static R reverse_arcs(R a){unsigned f=(1u<<a.size())-1;for(int i=0;i<(int)a.size();i++)a[i]^=f^(1u<<i);return a;}
static R leaf(const R& a,int q,bool rev){return rev?reverse_arcs(single(reverse_arcs(a),q)):single(a,q);}
static int mask(const P& f,const P& r,bool st){int m=0;for(auto a:{&f,&r})for(int q=0;q<(int)a->d.size();q++)if(!st||a->stable[q]){
 if(a->d[q]==1)m|=1;
 if(a->d[q]==-1)m|=2;
}return m;}

int main(){
 vector<R> seeds;
 constexpr int n=7;
 // Each permutation is a bit permutation possibly with sign changes.
 vector<array<unsigned,21>> bits;
 vector<unsigned> flip;
 array<int,n> p;
 iota(p.begin(),p.end(),0);
 int idx[n][n]={},e=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)idx[i][j]=idx[j][i]=e++;
 do{
  array<unsigned,21> b{};unsigned f=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   unsigned c=1u<<idx[p[i]][p[j]];
   b[idx[i][j]]=c;if(p[i]>p[j])f|=c;
  }
  bits.push_back(b);flip.push_back(f);
 }while(next_permutation(p.begin(),p.end()));
 vector<bool> seen(1u<<21);
 for(unsigned code=0;code<(1u<<21);code++)if(!seen[code]){
   R a(n);
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)
     if(code>>idx[i][j]&1)a[i]|=1u<<j;
     else a[j]|=1u<<i;
   seeds.push_back(a);
   for(size_t k=0;k<bits.size();k++){
    unsigned mask=code, z=flip[k];
    while(mask){z^=bits[k][__builtin_ctz(mask)];mask&=mask-1;}
    seen[z]=true;
   }
 }
 cout<<"n7 isomorphism representatives: "<<seeds.size()<<endl;
 auto scan=[&](int n){
  map<int,int> allReady, sameSite, both, stableSurv,rawSurv,rawFull;
  unsigned considered=0, strong=0;
  for(auto& a:seeds)for(unsigned m=0;m<(n==7?1u:128u);m++){
   R b=a;
   if(n==8){
    b.push_back(m);
    for(int i=0;i<7;i++)if(!(m>>i&1))b[i]|=1u<<7;
   }
   considered++;
   auto f=profile(b); if(!f.strong)continue;strong++;
   auto r=profile(reverse_arcs(b));swap(f.stable,r.stable);
   bool ready=false,same=false;
   for(int i=0;i<n;i++)for(int j=0;j<n;j++)
     if(f.stable[i]&&r.stable[j]&&abs(f.d[i])==1&&f.d[i]+r.d[j]==0)
       {ready=true;same|=i==j;}
   if(ready){
   if(!allReady[f.H]++) {
    cout<<"opposite ready H="<<f.H<<" rows";
    for(auto x:b)cout<<" "<<x;cout<<" stable unit ports";
    for(int i=0;i<n;i++) {
       if(f.stable[i]&&abs(f.d[i])==1)cout<<" F"<<i<<":"<<f.d[i];
       if(r.stable[i]&&abs(r.d[i])==1)cout<<" R"<<i<<":"<<r.d[i];
    }
    cout<<endl;
   }
   }
   if(same)sameSite[f.H]++;
   if(mask(f,r,false)!=3)continue;
   both[f.H]++;int sr=0, rr=0;
   for(int rev=0;rev<=1;rev++)for(int q=0;q<n;q++){
     const auto& s=rev?r:f;
     if(abs(s.d[q])!=1)continue;
     auto x=leaf(b,q,rev);auto xf=profile(x), xr=profile(reverse_arcs(x));swap(xf.stable,xr.stable);
     if(mask(xf,xr,false)==3) rr|=s.d[q]==1?1:2;
     if(s.stable[q]&&mask(xf,xr,true)==3)sr|=s.d[q]==1?1:2;
   }
   if(sr)stableSurv[f.H]++;
   if(rr)rawSurv[f.H]++;
   if(rr==3){
    if(!rawFull[f.H]++) {cout<<"depth-two raw H="<<f.H<<" rows";for(auto x:b)cout<<" "<<x;cout<<endl;}
   }
  }
  cout<<"n="<<n<<" considered="<<considered<<" strong="<<strong<<endl;
  for(auto [h,k]:both)cout<<h<<" both="<<k<<" opposite-ready="<<allReady[h]<<" same-site="<<sameSite[h]<<" stable-survivor="<<stableSurv[h]<<" raw-survivor="<<rawSurv[h]<<" depth-two="<<rawFull[h]<<endl;
 };
 scan(7);scan(8);
}
