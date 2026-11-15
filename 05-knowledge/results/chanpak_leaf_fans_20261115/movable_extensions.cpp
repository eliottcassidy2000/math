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

static pair<int,vector<int>> profile(const R& a){
 int n=a.size(), full=(1<<n)-1, dp[512][9]={}, h[512]={};
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
 return {h[full],d};
}

static R single(const R& a,int q){
 int n=a.size();
 R b=a;
 for(int i=0;i<n;i++)if(i!=q)b[i]|=1u<<n;
 b.push_back(1u<<q); return b;
}
static bool mixed(const vector<int>& d){
 return find(d.begin(),d.end(),1)!=d.end() && find(d.begin(),d.end(),-1)!=d.end();
}
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
  map<int,int> counts, positive, negative, both;
  unsigned considered=0;
  for(auto& a:seeds)for(unsigned m=0;m<(n==7?1u:128u);m++){
   R b=a;
   if(n==8){
    b.push_back(m);
    for(int i=0;i<7;i++)if(!(m>>i&1))b[i]|=1u<<7;
   }
   considered++;
   auto [h,d]=profile(b);
   if(!mixed(d))continue;
   counts[h]++;
   bool pm=false,nm=false;
   for(int q=0;q<n;q++)if(abs(d[q])==1&&mixed(profile(single(b,q)).second)){
     if(d[q]==1)pm=true;else nm=true;
   }
   if(pm){
    cout<<"positive survivor H="<<h<<" rows=";
    for(auto i:b)cout<<i<<",";
    cout<<endl;
    positive[h]++;
   }
   if(nm)negative[h]++;
   if(pm&&nm)both[h]++;
  }
  cout<<"n="<<n<<" considered="<<considered
      <<" (n8 source repeats possible)"<<endl;
  for(auto [h,k]:counts)cout<<h<<": "<<k<<" positive "<<positive[h]
      <<" negative "<<negative[h]<<" both "<<both[h]<<endl;
 };
 scan(7);scan(8);
}
