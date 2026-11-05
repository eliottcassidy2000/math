#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
using namespace std;
using Rows=vector<int>;
const int cap=541;
int forward_hp(const Rows&r){int n=r.size(),N=1<<n;vector<array<int,12>>d(N);for(int v=0;v<n;v++)d[1<<v][v]=1;for(int mask=1;mask<N;mask++){int total=0;for(int v=0;v<n;v++)total+=d[mask][v];if(total>cap)return cap+1;for(int v=0;v<n;v++)if(d[mask][v]){int z=r[v]&~mask;while(z){int bit=z&-z;z-=bit;d[mask|bit][__builtin_ctz((unsigned)bit)]+=d[mask][v];}}}int h=0;for(int v=0;v<n;v++)h+=d[N-1][v];return h;}
int backward_hp(const Rows&r){int n=r.size(),N=1<<n;vector<array<int,12>>d(N);vector<int>in(n);for(int u=0;u<n;u++)for(int v=0;v<n;v++)if(r[u]>>v&1)in[v]|=1<<u;for(int mask=1;mask<N;mask++){int total=0;for(int v=0;v<n;v++)if(mask>>v&1){int rem=mask^(1<<v);if(!rem)d[mask][v]=1;else{int z=in[v]&rem;while(z){int bit=z&-z;z-=bit;d[mask][v]+=d[rem][__builtin_ctz((unsigned)bit)];}}total+=d[mask][v];}if(total>cap)return cap+1;}int h=0;for(int v=0;v<n;v++)h+=d[N-1][v];return h;}
long long nodes[5]={},counts=0;map<int,long long>hist[5];ofstream bank;
void list(ostream&o,const Rows&r){o<<'[';for(int j=0;j<r.size();j++){if(j)o<<',';o<<r[j];}o<<']';}
void dfs(Rows const&r,Rows const&p,int depth,int h){nodes[depth]++;hist[depth][h]++;if(h==541)throw runtime_error("target541 reached");bank<<"{\"depth\":"<<depth<<",\"H\":"<<h<<",\"ports\":";list(bank,p);bank<<",\"rows\":";list(bank,r);bank<<"}\n";if(depth==4)return;int n=r.size();for(int side=0;side<2;side++)for(int mask=0;mask<8;mask++){Rows rr=r;rr.push_back(0);int sel=0;for(int j=0;j<3;j++)if(mask>>j&1)sel|=1<<p[j];for(int v=0;v<n;v++){bool newwins=((sel>>v)&1)^side;if(newwins)rr[n]|=1<<v;else rr[v]|=1<<n;}int a=forward_hp(rr),b=backward_hp(rr);counts++;if(a!=b)throw runtime_error("independent endpoint recurrences disagree");if(a>cap)continue;for(int drop=0;drop<4;drop++){Rows pp=p;pp.push_back(n);pp.erase(pp.begin()+drop);dfs(rr,pp,depth+1,a);}}}
int main(int argc,char**argv){if(argc!=2)return 2;bank.open(argv[1]);if(!bank)return 2;Rows seed={2,12,145,149,99,143,175,19};if(forward_hp(seed)!=179||backward_hp(seed)!=179)return 3;dfs(seed,{0,1,4},0,179);cout<<"{\"accepted_word_prefixes\":[";for(int d=0;d<=4;d++){if(d)cout<<',';cout<<nodes[d];}cout<<"],\"independent_count_comparisons\":"<<counts<<",\"histograms\":[";for(int d=0;d<=4;d++){if(d)cout<<',';cout<<'{';bool first=true;for(auto[h,c]:hist[d]){if(!first)cout<<',';first=false;cout<<'\"'<<h<<"\":"<<c;}cout<<'}';}cout<<"]}\n";}
