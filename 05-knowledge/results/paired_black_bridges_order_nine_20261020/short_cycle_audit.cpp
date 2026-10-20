#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 int n=argc>1?stoi(argv[1]):9;string pre=argc>2?argv[2]:"census_n"+to_string(n);auto st=chrono::steady_clock::now();int m=(n-1)*(n-2)/2;uint32_t full=(1u<<m)-1;vector<int> mp(1u<<m);ifstream f(pre+".map",ios::binary);f.read((char*)mp.data(),mp.size()*4);assert(f.gcount()==(streamsize)mp.size()*4);int N=*max_element(mp.begin(),mp.end())+1;
 int ix[16][16];vector<pair<int,int>>co;for(int y=1;y<n-1;y++)for(int x=n;x>y+1;x--){ix[x][y]=co.size();co.push_back({x,y});}vector<int>sig;for(auto[x,y]:co)sig.push_back(ix[n+1-y][n+1-x]);int s=m/2;vector<uint32_t>l(1u<<s),h(1u<<(m-s));for(uint32_t t=1;t<l.size();t++){int k=__builtin_ctz(t);l[t]=l[t&(t-1)]|(1u<<sig[k]);}for(uint32_t t=1;t<h.size();t++){int k=__builtin_ctz(t);h[t]=h[t&(t-1)]|(1u<<sig[k+s]);}
 vector<uint64_t> keys;keys.reserve(1u<<(m-2));vector<uint64_t>deg(N);uint64_t loops=0;
 for(uint32_t t=0;t<(1u<<(m-1));t++){uint32_t rt=l[t&((1u<<s)-1)]|h[t>>s];if(t==rt||t>rt||t>(full^rt))continue;int a=mp[t],b=mp[full^t];deg[a]++;deg[b]++;if(a==b){loops++;continue;}if(a>b)swap(a,b);keys.push_back((uint64_t(a)<<32)|uint32_t(b));}
 mp.clear();mp.shrink_to_fit();sort(keys.begin(),keys.end());vector<uint64_t> unique;vector<uint32_t>counts;unique.reserve(keys.size());counts.reserve(keys.size());for(auto k:keys){if(unique.empty()||k!=unique.back()){unique.push_back(k);counts.push_back(1);}else counts.back()++;}keys.clear();keys.shrink_to_fit();cerr<<"support pairs="<<unique.size()<<" seconds="<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<endl;
 vector<uint32_t>off(N+1);for(auto k:unique){off[(k>>32)+1]++;off[uint32_t(k)+1]++;}for(int i=1;i<=N;i++)off[i]+=off[i-1];vector<uint32_t>adj(off.back()),cur=off;for(auto k:unique){uint32_t a=k>>32,b=k;adj[cur[a]++]=b;adj[cur[b]++]=a;}for(int i=0;i<N;i++)sort(adj.begin()+off[i],adj.begin()+off[i+1]);
 auto common=[&](uint32_t a,uint32_t b,uint32_t exa,uint32_t exb)->int{uint32_t i=off[a],j=off[b],ie=off[a+1],je=off[b+1];while(i<ie&&j<je){if(adj[i]<adj[j])i++;else if(adj[i]>adj[j])j++;else{uint32_t z=adj[i];if(z!=exa&&z!=exb)return z;i++;j++;}}return -1;};
 uint64_t paral=0,tri=0,square=0,leaf=0,uncovered=0;ofstream quads(pre+"_quadrilaterals.tsv"),fail(pre+"_short_cycle_failures.tsv");quads<<"u\tv\tx\ty\n";fail<<"u\tv\n";
 for(size_t i=0;i<unique.size();i++){uint32_t a=unique[i]>>32,b=unique[i];if(counts[i]>1){paral+=counts[i];continue;}if(min(deg[a],deg[b])==1){leaf++;continue;}if(common(a,b,a,b)>=0){tri++;continue;}bool found=false;if(off[a+1]-off[a]>off[b+1]-off[b])swap(a,b);for(uint32_t j=off[a];j<off[a+1];j++){uint32_t x=adj[j];if(x==b)continue;int y=common(x,b,a,b);if(y>=0){found=true;square++;quads<<a<<'\t'<<b<<'\t'<<x<<'\t'<<y<<'\n';break;}}if(!found){uncovered++;fail<<a<<'\t'<<b<<'\n';}}
 ofstream out(pre+"_short_cycle_audit.json");out<<"{\"n\":"<<n<<",\"support_nonloop_pairs\":"<<unique.size()<<",\"loop_edges\":"<<loops<<",\"digon_edges\":"<<paral<<",\"triangle_edges\":"<<tri<<",\"quadrilateral_edges\":"<<square<<",\"leaf_edges\":"<<leaf<<",\"uncovered_edges\":"<<uncovered<<",\"elapsed_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"}\n";cerr<<"DONE triangle="<<tri<<" square="<<square<<" leaf="<<leaf<<" uncovered="<<uncovered<<" seconds="<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<endl;
}
