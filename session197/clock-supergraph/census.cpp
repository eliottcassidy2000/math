#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <unordered_map>
using namespace std;
using Rows=vector<int>;
const int cap=541;
int forward_hp(const Rows&r){int n=r.size(),N=1<<n;vector<array<int,12>>d(N);for(int v=0;v<n;v++)d[1<<v][v]=1;for(int mask=1;mask<N;mask++){int total=0;for(int v=0;v<n;v++)total+=d[mask][v];if(total>cap)return cap+1;for(int v=0;v<n;v++)if(d[mask][v]){int z=r[v]&~mask;while(z){int bit=z&-z;z-=bit;d[mask|bit][__builtin_ctz((unsigned)bit)]+=d[mask][v];}}}int h=0;for(int v=0;v<n;v++)h+=d[N-1][v];return h;}
int backward_hp(const Rows&r){int n=r.size(),N=1<<n;vector<array<int,12>>d(N);vector<int>in(n);for(int u=0;u<n;u++)for(int v=0;v<n;v++)if(r[u]>>v&1)in[v]|=1<<u;for(int mask=1;mask<N;mask++){int total=0;for(int v=0;v<n;v++)if(mask>>v&1){int rem=mask^(1<<v);if(!rem)d[mask][v]=1;else{int z=in[v]&rem;while(z){int bit=z&-z;z-=bit;d[mask][v]+=d[rem][__builtin_ctz((unsigned)bit)];}}total+=d[mask][v];}if(total>cap)return cap+1;}int h=0;for(int v=0;v<n;v++)h+=d[N-1][v];return h;}

struct Hash{size_t operator()(Rows const&r)const{size_t h=0;for(auto x:r)h=(h*1000003)^x;return h;}};
Rows canon(Rows const&r){int n=r.size();vector<int>p(n);iota(p.begin(),p.end(),0);Rows best;do{Rows out(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(r[p[i]]>>p[j]&1)out[i]|=1<<j;if(best.empty()||out<best)best=out;}while(next_permutation(p.begin()+8,p.end()));return best;}
void writebank(ostream&o,unordered_map<Rows,int,Hash>const&bank){vector<pair<Rows,int>>sorted(bank.begin(),bank.end());sort(sorted.begin(),sorted.end());for(auto const&[r,h]:sorted){o<<"{\"n\":"<<r.size()<<",\"H\":"<<h<<",\"rows\":[";for(int i=0;i<r.size();i++){if(i)o<<',';o<<r[i];}o<<"]}\n";}}
int main(int argc,char**argv){if(argc!=2)return 2;ofstream bank(argv[1]);if(!bank)return 2;unordered_map<Rows,int,Hash>cur{{Rows{2,12,145,149,99,143,175,19},179}};writebank(bank,cur);cout<<"{\"generation\":[";for(int n=9;n<=12;n++){unordered_map<Rows,int,Hash>nxt;long long tested=0,accepted=0;for(auto const&[rows,h]:cur)for(int mask=0;mask<(1<<(n-1));mask++){Rows r=rows;r.push_back(mask);for(int j=0;j<n-1;j++)if(!(mask>>j&1))r[j]|=1<<(n-1);tested++;int H=forward_hp(r),H2=backward_hp(r);if(H!=H2)throw runtime_error("forward/backward disagreement");if(H>541)continue;accepted++;if(H==541)throw runtime_error("H541 supergraph found");nxt.emplace(canon(r),H);}if(n>9)cout<<',';cout<<"{\"n\":"<<n<<",\"candidates\":"<<tested<<",\"accepted_before_isomorphism\":"<<accepted<<",\"seed_fixed_classes\":"<<nxt.size()<<'}';cur=move(nxt);writebank(bank,cur);}cout<<"],\"target_hits\":0}\n";}
