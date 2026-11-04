#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>
using namespace std;
void canonical(const vector<unsigned>&r){int n=r.size();vector<vector<int>>groups(n);for(int i=0;i<n;i++)groups[__builtin_popcount(r[i])].push_back(i);vector<int>p;uint32_t best=numeric_limits<uint32_t>::max();uint64_t aut=0;function<void(int)> rec=[&](int d){if(d==n){uint32_t v=0;int bit=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,bit++)if(r[p[i]]>>p[j]&1)v|=1u<<bit;if(v<best){best=v;aut=1;}else if(v==best)aut++;return;}auto g=groups[d];do{p.insert(p.end(),g.begin(),g.end());rec(d+1);p.resize(p.size()-g.size());}while(next_permutation(g.begin(),g.end()));};rec(0);cout<<best<<" "<<aut<<"\n";}
void joint(const vector<unsigned>&r){int n=r.size(),full=(1<<n)-1;for(int s=0;s<n;s++){vector<uint64_t>d((1<<n)*n);d[(1<<s)*n+s]=1;for(int mask=1;mask<=full;mask++){unsigned a=mask;while(a){int v=__builtin_ctz(a);a&=a-1;auto val=d[mask*n+v];if(!val)continue;unsigned b=r[v]&~mask;while(b){int w=__builtin_ctz(b);b&=b-1;d[(mask|(1<<w))*n+w]+=val;}}}for(int t=0;t<n;t++)cout<<d[full*n+t]<<" ";}cout<<"\n";}
int main(int argc,char**argv){bool endpoint=argc>1&&(string(argv[1])=="joint");int k;cin>>k;for(int i=0;i<k;i++){int n;cin>>n;vector<unsigned>r(n);for(auto&v:r)cin>>v;if(endpoint)joint(r);else canonical(r);}}
