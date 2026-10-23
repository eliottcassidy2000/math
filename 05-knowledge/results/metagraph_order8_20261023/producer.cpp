// Native order-eight ownership and complete line/port replay.
// Every certificate source is a literal Hamiltonian path of a certified class representative.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr int n=8,C=6880,M=21,N=1<<M,ALL=N-1;
struct Class {uint32_t canon;int H,aut,h,b,conv;array<int,8> adj;};
vector<Class> cs; vector<uint16_t> owner(N,65535);vector<uint8_t> visits(N,0);vector<pair<int,int>> tiles;int sigmaBit[M];uint64_t checks=0,path_total=0;
void need(bool x,const string& why){checks++;if(!x)throw runtime_error(why);}
vector<string> split(const string&s){vector<string>o;string x;istringstream z(s);while(getline(z,x,'\t'))o.push_back(x);return o;}
uint32_t reflect(uint32_t m){uint32_t s=0;for(int i=0;i<M;i++)if(m>>i&1)s|=1u<<sigmaBit[i];return s;}
uint32_t canonical(const Class &c,int &aut){
 array<int,8>s,p;for(int i=0;i<n;i++){s[i]=__builtin_popcount((unsigned)c.adj[i]);p[i]=i;}
 sort(p.begin(),p.end(),[&](int a,int b){return s[a]!=s[b]?s[a]<s[b]:a<b;});vector<pair<int,int>>g;
 for(int a=0;a<n;){int b=a+1;while(b<n&&s[p[b]]==s[p[a]])b++;g.push_back({a,b});a=b;}
 uint32_t best=~0u;aut=0;
 while(true){uint32_t z=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)z=z*2+((c.adj[p[i]]>>p[j])&1);if(z<best){best=z;aut=1;}else if(z==best)aut++;
 int k=(int)g.size()-1;while(k>=0&&!next_permutation(p.begin()+g[k].first,p.begin()+g[k].second))k--;if(k<0)break;}
 return best;
}
int path_dp(const Class&c){int d[256][8]={};for(int i=0;i<n;i++)d[1<<i][i]=1;for(int m=1;m<256;m++)for(int u=0;u<n;u++)if(d[m][u]){int rem=c.adj[u]&~m;while(rem){int v=__builtin_ctz((unsigned)rem);rem&=rem-1;d[m|1<<v][v]+=d[m][u];}}return accumulate(d[255],d[255]+8,0);}
void paths(int ci,int used,int len,array<int,8>&P,int&count){
 if(len==n){uint32_t mask=0;for(int k=0;k<M;k++)if(cs[ci].adj[P[tiles[k].first]]>>P[tiles[k].second]&1)mask|=1u<<k;
  need(owner[mask]==65535||owner[mask]==ci,"source owner collision");owner[mask]=ci;need(visits[mask]<255,"visit counter overflow");visits[mask]++;count++;path_total++;return;}
 int rem=len?cs[ci].adj[P[len-1]]&~used:255;while(rem){int v=__builtin_ctz((unsigned)rem);rem&=rem-1;P[len]=v;paths(ci,used|1<<v,len+1,P,count);}
}
void u16(ofstream&f,unsigned x){f.put(x&255);f.put((x>>8)&255);}void u32(ofstream&f,uint32_t x){for(int j=0;j<4;j++)f.put((x>>(8*j))&255);}void u64(ofstream&f,uint64_t x){for(int j=0;j<8;j++)f.put((x>>(8*j))&255);}
void write_lines(vector<uint64_t>&v,const string&name){
 sort(v.begin(),v.end());ofstream f(name,ios::binary);f<<"source\ttarget\tblue\tblack\n";size_t i=0;
 while(i<v.size()){uint64_t pair=v[i]>>1;unsigned counts[2]={};while(i<v.size()&&(v[i]>>1)==pair){counts[v[i]&1]++;i++;}f<<(pair>>16)<<'\t'<<(pair&65535)<<'\t'<<counts[1]<<'\t'<<counts[0]<<'\n';}
}
int main(int argc,char**argv){try{
 need(argc==3,"usage producer classes.tsv output-directory");string dir=argv[2];ifstream f(argv[1]);need(bool(f),"class input missing");string line;getline(f,line);uint64_t labelled_mass=0;
 while(getline(f,line)){auto r=split(line);need(r.size()>=11,"bad class row");need(stoi(r[0])==(int)cs.size(),"class ids not sequential");Class c{};c.canon=stoul(r[1]);c.H=stoi(r[2]);c.aut=stoi(r[3]);c.h=stoi(r[4]);c.b=stoi(r[5]);c.conv=stoi(r[8]);int k=27;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k--){if(c.canon>>k&1)c.adj[i]|=1<<j;else c.adj[j]|=1<<i;}
 int aut=0;need(canonical(c,aut)==c.canon&&aut==c.aut,"class canonical/aut mismatch");need(path_dp(c)==c.H,"class Hamilton DP mismatch");need(c.h*c.aut==c.H,"class fibre mismatch");labelled_mass+=40320/c.aut;cs.push_back(c);}
 need(cs.size()==C&&labelled_mass==(1u<<28),"class universe mismatch");vector<uint32_t>can;for(auto&c:cs)can.push_back(c.canon);sort(can.begin(),can.end());need(adjacent_find(can.begin(),can.end())==can.end(),"duplicate classes");
 for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)tiles.push_back({i,j});for(int k=0;k<M;k++)sigmaBit[k]=find(tiles.begin(),tiles.end(),make_pair(n-1-tiles[k].second,n-1-tiles[k].first))-tiles.begin();
 for(int ci=0;ci<C;ci++){array<int,8>P{};int count=0;paths(ci,0,0,P,count);need(count==cs[ci].H,"DFS H mismatch");}
 cerr<<"native paths "<<path_total<<" generated\n";
 vector<int>fib(C),bluefib(C),merged_id(C),rep;for(int ci=0;ci<C;ci++)if(ci<=cs[ci].conv){merged_id[ci]=rep.size();rep.push_back(ci);}for(int ci=0;ci<C;ci++)merged_id[ci]=merged_id[min(ci,cs[ci].conv)];
 need(rep.size()==3528,"merged node count");ofstream own(dir+"/owner_map.u16",ios::binary);vector<uint64_t>ordinary,merged,ports;ordinary.reserve(N/2);merged.reserve(N/2);ports.reserve(N/2);uint64_t blue=0,black=0;int apex=find(tiles.begin(),tiles.end(),make_pair(0,7))-tiles.begin();
 for(uint32_t m=0;m<N;m++){int ci=owner[m];need(ci<C,"uncovered tiling");need(visits[m]==cs[ci].aut,"path orbit multiplicity mismatch");u16(own,ci);fib[ci]++;uint32_t sm=reflect(m);bool isblue=m==sm;bluefib[ci]+=isblue;need(owner[sm]==cs[ci].conv,"literal converse reflection mismatch");uint32_t cm=m^ALL;int cj=owner[cm];
  if(m<cm){int a=min(ci,cj),b=max(ci,cj);ordinary.push_back((((uint64_t)a<<16)|b)*2+isblue);a=min(merged_id[ci],merged_id[cj]);b=max(merged_id[ci],merged_id[cj]);merged.push_back((((uint64_t)a<<16)|b)*2+isblue);(isblue?blue:black)++;}
  if(!(m>>apex&1)){int a=1,b=0;for(int k=0;k<M;k++){auto [i,j]=tiles[k];if(i==0&&(m>>k&1))a++;if(j==7&&!(m>>k&1))b++;}uint64_t key=(((uint64_t)ci*C+cj)*8+a)*8+b;key=key*2+isblue;ports.push_back(key);}
 }
 own.close();for(int i=0;i<C;i++)need(fib[i]==cs[i].h&&bluefib[i]==cs[i].b,"native fibre/blue mismatch");need(blue==2048&&black==1046528&&ports.size()==1048576,"line totals");
 write_lines(ordinary,dir+"/ordinary_edges.tsv");write_lines(merged,dir+"/merged_edges.tsv");sort(ports.begin(),ports.end());ofstream pf(dir+"/ports.bin",ios::binary);uint64_t port_groups=0;for(size_t i=0;i<ports.size();){size_t j=i+1;while(j<ports.size()&&ports[j]==ports[i])j++;u64(pf,ports[i]);u32(pf,j-i);i=j;port_groups++;}pf.close();
 ofstream nf(dir+"/merged_nodes.tsv",ios::binary);nf<<"merged_id\tclass_a\tclass_b\tcanon_a\tcanon_b\tfibre\tblue\n";for(int j=0;j<(int)rep.size();j++){int a=rep[j],b=cs[a].conv;nf<<j<<'\t'<<a<<'\t'<<b<<'\t'<<cs[a].canon<<'\t'<<cs[b].canon<<'\t'<<(a==b?fib[a]:fib[a]+fib[b])<<'\t'<<(a==b?bluefib[a]:0)<<'\n';}
 ofstream report(dir+"/producer_summary.json",ios::binary);report<<"{\"checks\":"<<checks<<",\"paths\":"<<path_total<<",\"tilings\":"<<N<<",\"ordinary_classes\":"<<C<<",\"merged_nodes\":"<<rep.size()<<",\"blue_lines\":"<<blue<<",\"black_lines\":"<<black<<",\"port_groups\":"<<port_groups<<"}\n";cout<<"checks "<<checks<<" port_groups "<<port_groups<<"\n";
 }catch(const exception&e){cerr<<e.what()<<'\n';return 1;}return 0;}
