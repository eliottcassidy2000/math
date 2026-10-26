#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
struct DSU{vector<uint32_t> p,sz;DSU(uint32_t n):p(n),sz(n,1){iota(p.begin(),p.end(),0);}uint32_t find(uint32_t x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;}void join(uint32_t a,uint32_t b){a=find(a);b=find(b);if(a==b)return;if(sz[a]<sz[b])swap(a,b);p[b]=a;sz[a]+=sz[b];}};
struct Map{uint32_t base=0,ones=0,zeros=0;array<array<uint32_t,256>,4> lut{};uint32_t apply(uint32_t x)const{return base^lut[0][x&255]^lut[1][(x>>8)&255]^lut[2][(x>>16)&255]^lut[3][(x>>24)&255];}};
int main(int argc,char**argv){int maxn=argc>1?stoi(argv[1]):8;need(maxn>=1&&maxn<=8,"limit 1..8");string out=argc>2?argv[2]:"";cout<<"[\n";
for(int n=1;n<=maxn;n++){
 vector<pair<int,int>> edges;int idx[9][9];for(auto &row:idx)fill(begin(row),end(row),-1);for(int i=0;i<n;i++)for(int j=i+2;j<n;j++){idx[i][j]=edges.size();edges.push_back({i,j});}int m=edges.size();uint32_t N=1u<<m;
 auto make_map=[&](vector<int> p,bool converse){Map f;for(int k=0;k<m;k++){auto[i,j]=edges[k];int u=p[i],v=p[j];bool reverse=u>v;if(reverse)swap(u,v);bool flip=reverse^converse;int old=idx[u][v];if(old<0){if(!flip)f.base|=1u<<k;}else{if(flip)f.base|=1u<<k;for(int z=0;z<256;z++)if(z>>(old%8)&1)f.lut[old/8][z]^=1u<<k;}}return f;};
 vector<Map> moves;for(int a=0;a<n;a++)for(int b=a+2;b<n;b++){vector<int>p;for(int i=0;i<a;i++)p.push_back(i);for(int i=a+1;i<=b;i++)p.push_back(i);p.push_back(a);for(int i=b+1;i<n;i++)p.push_back(i);Map f=make_map(p,false);f.zeros=1u<<idx[a][b];if(a)f.ones|=1u<<idx[a-1][a+1];if(b+1<n)f.ones|=1u<<idx[a][b+1];moves.push_back(f);}
 vector<int> rev(n);for(int i=0;i<n;i++)rev[i]=n-1-i;Map reflect=make_map(rev,true);DSU d(N);uint64_t legal=0;for(uint32_t x=0;x<N;x++)for(auto&f:moves)if(!(x&f.zeros)&&(x&f.ones)==f.ones){uint32_t y=f.apply(x);need(y<N,"range");d.join(x,y);legal++;}
 vector<uint32_t> ordinary(N),ordinary_size(N);map<uint32_t,uint32_t> rep;map<uint32_t,uint32_t> hist;for(uint32_t x=0;x<N;x++){ordinary[x]=d.find(x);auto it=rep.find(ordinary[x]);if(it==rep.end())rep[ordinary[x]]=x;}
 for(auto[k,v]:rep){hist[d.sz[k]]++;ordinary_size[k]=d.sz[k];}
 for(uint32_t x=0;x<N;x++){uint32_t y=reflect.apply(x);need(reflect.apply(y)==x,"reflection involution");d.join(x,y);}
 map<uint32_t,uint32_t> merged;uint64_t sc=0;for(auto[k,x]:rep){uint32_t r=d.find(x);if(!merged.count(r))merged[r]=x;if(ordinary[reflect.apply(x)]==k)sc++;}
 need(2*merged.size()==rep.size()+sc,"converse count");
 uint64_t blue=0,black=0,loops=0;map<array<uint32_t,3>,uint32_t> lines;vector<uint64_t> degree(N,0);if(m){for(uint32_t x=0;x<N;x++){uint32_t y=x^(N-1);if(x>y)continue;bool b=reflect.apply(x)==x;need((reflect.apply(y)==y)==b,"colour paired");uint32_t u=merged[d.find(x)],v=merged[d.find(y)];if(u>v)swap(u,v);lines[{u,v,uint32_t(b)}]++;degree[d.find(x)]++;degree[d.find(y)]++;if(b)blue++;else black++;if(u==v)loops++;}for(auto[r,x]:merged)need(degree[r]==d.sz[r],"half-edge degree");need(blue+black==N/2,"line total");}
 if(!out.empty()){
 ofstream f(out+"/nodes_"+to_string(n)+".tsv");f<<"ordinary_representative\ttilings_tc\tmerged_representative\tself_converse\n";for(auto[k,x]:rep)f<<x<<'\t'<<ordinary_size[k]<<'\t'<<merged[d.find(x)]<<'\t'<<(ordinary[reflect.apply(x)]==k)<<'\n';
 if(n<=(argc>3?stoi(argv[3]):7)){ofstream e(out+"/lines_"+to_string(n)+".tsv");e<<"u\tv\tblue\tmultiplicity\n";for(auto[k,c]:lines)e<<k[0]<<'\t'<<k[1]<<'\t'<<k[2]<<'\t'<<c<<'\n';}
 }
 cout<<"{\"n\":"<<n<<",\"tilings\":"<<N<<",\"ordinary_classes\":"<<rep.size()<<",\"merged_classes\":"<<merged.size()<<",\"self_converse\":"<<sc<<",\"legal_right_relocations\":"<<legal<<",\"blue_lines\":"<<blue<<",\"black_lines\":"<<black<<",\"loop_lines\":"<<loops<<",\"distinct_coloured_endpoint_pairs\":"<<lines.size()<<",\"tc_histogram\":{";bool first=true;for(auto[k,c]:hist){if(!first)cout<<',';first=false;cout<<'"'<<k<<"\":"<<c;}cout<<"}}"<<(n<maxn?",":"")<<'\n'<<flush;
}cout<<"]\n";
}
