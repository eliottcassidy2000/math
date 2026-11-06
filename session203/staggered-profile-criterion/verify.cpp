#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;using I=long long;using G=vector<unsigned>;
void need(bool ok,const char*s){if(!ok)throw runtime_error(s);}
void valid(const G&a){int n=a.size();need(n>0&&n<31,"order");for(int i=0;i<n;i++){need(a[i]<(1u<<n)&&!(a[i]>>i&1),"row");for(int j=0;j<i;j++)need(((a[i]>>j&1)+(a[j]>>i&1))==1,"pair");}}
struct Data{vector<I>h;vector<vector<I>>s,e;};
Data profile(const G&a){valid(a);int n=a.size(),N=1<<n;Data z{vector<I>(N),vector<vector<I>>(N,vector<I>(n)),vector<vector<I>>(N,vector<I>(n))};z.h[0]=1;for(int i=0;i<n;i++)z.s[1<<i][i]=z.e[1<<i][i]=1;for(int S=1;S<N;S++)for(int i=0;i<n;i++)if(S>>i&1){int T=S^(1<<i);if(T)for(int j=0;j<n;j++)if(T>>j&1){if(a[i]>>j&1)z.s[S][i]+=z.s[T][j];else z.e[S][i]+=z.e[T][j];}z.h[S]+=z.s[S][i];}return z;}
G child(const G&a,int r){vector<int>order{r};for(int i=0;i<(int)a.size();i++)if(i!=r)order.push_back(i);G u(a.size());for(int i=0;i<(int)a.size();i++)for(int j=0;j<(int)a.size();j++)if(a[order[i]]>>order[j]&1)u[i]|=1u<<j;G out{2,8,3,5};out.resize(a.size()+3);for(int j=1;j<(int)u.size();j++){int x=j+3;for(int i=0;i<4;i++)if(i||u[0]>>j&1)out[i]|=1u<<x;else out[x]|=1u<<i;for(int k=1;k<j;k++){int y=k+3;if(u[j]>>k&1)out[x]|=1u<<y;else out[y]|=1u<<x;}}valid(out);return out;}
int main(){map<int,G>seeds;ifstream sf("SEEDS.tsv");need(bool(sf),"missing seeds");int target,n;while(sf>>target>>n){G a(n);for(auto&r:a)need(bool(sf>>r),"seed row missing");need(seeds.emplace(target,a).second,"duplicate seed");need(profile(a).h.back()==target,"seed count");}need(seeds.size()==194,"seed count records");for(int h=23;h<=409;h+=2)need(seeds.count(h),"seed interval gap");
 ifstream cf("PROFILE_CERTIFICATES.tsv");need(bool(cf),"missing profiles");int B,r,j,root;set<pair<int,int>>keys;I maxcount=0;int graphs=194;while(cf>>B>>r>>j>>root>>n){G a(n);for(auto&row:a)need(bool(cf>>row),"profile row missing");need(B>=137&&B<=997&&B%2==1&&(r==1||r==2)&&j>=0&&j<=8&&root>=0&&root<n,"profile bounds");need(keys.emplace(B,r).second,"duplicate profile");auto z=profile(a);int full=(1<<n)-1;I h=z.h.back(),s=z.s.back()[root],L=0;for(int S=1;S<=full;S++)L+=z.e[S][root]*z.h[full^S];need(h+L-s==B-2*j&&s==r+3*j,"profile equations");auto c=child(a,root);I hc=profile(c).h.back();need(hc==3*B+2*r,"literal compiled count");maxcount=max(maxcount,hc);graphs+=2;}
 need(keys.size()==862,"profile records");for(int B=137;B<=997;B+=2)for(int r:{1,2})need(keys.count({B,r}),"profile interval gap");
 cout<<"{\"status\":\"PASS\",\"seed_interval\":[23,409],\"seed_graphs\":194,\"profile_base_interval\":[137,997],\"profile_base_filter\":\"odd\",\"both_residues_each_base\":true,\"J\":8,\"profile_certificates\":862,\"literal_graph_counts\":"<<graphs<<",\"largest_compiled_target\":"<<maxcount<<",\"all_order_profile_lemma\":\"OPEN\",\"eventual_coverage\":\"CONDITIONAL\"}\n";
}
