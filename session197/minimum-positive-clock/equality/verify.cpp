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
int matching(Rows const&r,int root){vector<int> O,I;for(int v=0;v<(int)r.size();v++)if(v!=root)((r[root]>>v&1)?O:I).push_back(v);vector<bool>a(1<<I.size());a[0]=true;for(int v:O){auto b=a;for(int m=0;m<(int)a.size();m++)if(a[m])for(int j=0;j<(int)I.size();j++)if(!(m>>j&1)&&(r[v]>>I[j]&1))b[m|1<<j]=true;a=b;}int ans=0;for(int m=0;m<(int)a.size();m++)if(a[m])ans=max(ans,__builtin_popcount((unsigned)m));return ans;}
Rows clone(Rows const&r,int root,int m){vector<int> old;for(int v=0;v<(int)r.size();v++)if(v!=root)old.push_back(v);int n=old.size(),N=n+m;Rows out(N);for(int i=0;i<n;i++){for(int j=0;j<n;j++)if(r[old[i]]>>old[j]&1)out[i]|=1<<j;if(r[old[i]]>>root&1)for(int a=n;a<N;a++)out[i]|=1<<a;}for(int a=n;a<N;a++){for(int j=0;j<n;j++)if(r[root]>>old[j]&1)out[a]|=1<<j;for(int b=a+1;b<N;b++)out[a]|=1<<b;}return out;}
void printrow(ostream&o,Rows const&r){o<<'[';for(int i=0;i<(int)r.size();i++){if(i)o<<',';o<<r[i];}o<<']';}
int main(int argc, char** argv){if(argc!=2){cerr<<"usage: verifier OUTPUT_BANK"<<endl;return 2;}const int cap=217;unordered_map<Rows,int,Hash> bank{{canonical({2,4,1}),3}};ofstream out(argv[1]);if(!out){cerr<<"cannot open bank output"<<endl;return 2;}long long hits=0,tests=0;for(int n=3;n<=10;n++){if(n>3){unordered_map<Rows,int,Hash> next;long long candidates=0,accepted=0;for(auto const&[p,h]:bank){for(int m=1;m<(1<<(n-1))-1;m++){candidates++;auto r=extend(p,m);int H=hp(r,cap);if(H>cap)continue;accepted++;next.emplace(canonical(r),H);}}bank=move(next);cerr<<"order "<<n<<" candidates "<<candidates<<" accepted "<<accepted<<" classes "<<bank.size()<<endl;}vector<pair<Rows,int>> sorted(bank.begin(),bank.end());sort(sorted.begin(),sorted.end());long long localtest=0,localhit=0;for(auto const&[r,H]:sorted){out<<"{\"n\":"<<n<<",\"H\":"<<H<<",\"rows\":";printrow(out,r);out<<"}\n";if(H%6!=1)continue;for(int root=0;root<n;root++){if(matching(r,root)!=2)continue;localtest++;int H2=hp(clone(r,root,2),3*H-2);if(H2!=3*H-2)continue;int H3=hp(clone(r,root,3),9*H-8);if(H3!=9*H-8)continue;cout<<"HIT A="<<(H-1)/3<<" n="<<n<<" root="<<root<<" rows=";printrow(cout,r);cout<<endl;localhit++;}}tests+=localtest;hits+=localhit;cerr<<"checked "<<n<<" rank2 target roots "<<localtest<<" hits "<<localhit<<endl;out.flush();}cout<<"DONE tested="<<tests<<" hits="<<hits<<endl;}
