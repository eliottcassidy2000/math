#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
using Rows=vector<uint16_t>;
struct Hash {size_t operator()(Rows const&r)const{size_t h=0;for(auto v:r)h=(h*1000003)^v;return h;}};
int hp(Rows const&r,int cap=1000000000){int n=r.size(),N=1<<n;vector<array<int,13>> d(N);for(int i=0;i<n;i++)d[1<<i][i]=1;for(int m=1;m<N;m++){int total=0;for(int i=0;i<n;i++)total+=d[m][i];if(total>cap)return cap+1;for(int i=0;i<n;i++)if(d[m][i]){int z=r[i]&~m;while(z){int b=z&-z;z-=b;d[m|b][__builtin_ctz((unsigned)b)]+=d[m][i];}}}return accumulate(d[N-1].begin(),d[N-1].begin()+n,0);}
Rows canonical(Rows const&r){int n=r.size();vector<vector<int>> groups(n);for(int i=0;i<n;i++)groups[__builtin_popcount(r[i])].push_back(i);vector<vector<int>> gs;for(auto g:groups)if(!g.empty())gs.push_back(g);Rows best;vector<int> p;function<void(int)> visit=[&](int k){if(k==(int)gs.size()){vector<int> inv(n);for(int j=0;j<n;j++)inv[p[j]]=j;Rows out(n);for(int j=0;j<n;j++){int z=r[p[j]];while(z){int b=z&-z;z-=b;out[j]|=1<<inv[__builtin_ctz((unsigned)b)];}}if(best.empty()||out<best)best=out;return;}auto g=gs[k];do{for(int v:g)p.push_back(v);visit(k+1);p.resize(p.size()-g.size());}while(next_permutation(g.begin(),g.end()));};visit(0);return best;}
Rows extend(Rows const&r,int mask){Rows s=r;int n=r.size();for(int j=0;j<n;j++)if(!(mask>>j&1))s[j]|=1<<n;s.push_back(mask);return s;}
void printrow(ostream&o,Rows const&r){o<<'[';for(int i=0;i<(int)r.size();i++){if(i)o<<',';o<<r[i];}o<<']';}

int main(int argc,char**argv) {
 if(argc!=2)return 2;
 ofstream bankout(argv[1]); if(!bankout)return 2;
 unordered_map<Rows,int,Hash> bank{{canonical({2,12,145,149,99,143,175,19}),179}};
 cout<<"{\"generation\":[";
 for(int n=8;n<=12;n++) {
  long long candidates=0,accepted=0;
  if(n>8) {
   unordered_map<Rows,int,Hash> next;
   for(auto const&[parent,h]:bank) {
    for(int mask=1;mask<(1<<(n-1))-1;mask++) {
     candidates++; auto rows=extend(parent,mask); int H=hp(rows,541);
     if(H>541)continue;
     accepted++; if(H==541)throw runtime_error("H541 induced strong extension");
     next.emplace(canonical(rows),H);
    }
   }
   bank=move(next);
  }
  if(n>8)cout<<',';
  cout<<"{\"order\":"<<n<<",\"candidates\":"<<candidates<<",\"accepted\":"<<accepted<<",\"strong_classes\":"<<bank.size()<<'}';
  vector<pair<Rows,int>> sorted(bank.begin(),bank.end()); sort(sorted.begin(),sorted.end());
  for(auto const&[rows,H]:sorted) {
   bankout<<"{\"n\":"<<n<<",\"H\":"<<H<<",\"rows\":"; printrow(bankout,rows); bankout<<"}\n";
  }
 }
 cout<<"],\"target_hits\":0}\n";
}
