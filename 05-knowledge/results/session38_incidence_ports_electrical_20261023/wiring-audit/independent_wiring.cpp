// Independent Session38 audit. No producer implementation imported.
// Class discovery: score-canonical orbit coverage. Mask classification: exact E
// signature evaluated by direct cyclic triples, unlike producer HP normalization.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
using A=array<uint8_t,8>; using Sig=array<uint64_t,8>;
uint64_t checks=0,canon_perms=0;
void need(bool q,string m){++checks;if(!q)throw runtime_error(m);}
int pc(unsigned x){return __builtin_popcount(x);}
bool arc(const A&a,int i,int j){return(a[i]>>j)&1;}
vector<string> split(string s){vector<string>o;stringstream z(s);string q;while(getline(z,q,'\t'))o.push_back(q);return o;}
A decode(uint32_t word){A a{};int k=27;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++,k--){if((word>>k)&1)a[i]|=1<<j;else a[j]|=1<<i;}return a;}
uint32_t encode(const A&a,const array<int,8>&p){uint32_t z=0;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)z=(z<<1)|arc(a,p[i],p[j]);return z;}
A converse(const A&a){A b;for(int i=0;i<8;i++)b[i]=(~a[i])&255&~(1<<i);return b;}
pair<uint32_t,int> canonical(const A&a){array<int,8>p{0,1,2,3,4,5,6,7};sort(p.begin(),p.end(),[&](int i,int j){return make_pair(pc(a[i]),i)<make_pair(pc(a[j]),j);});vector<pair<int,int>> blocks;for(int s=0,e;s<8;s=e){e=s+1;while(e<8&&pc(a[p[s]])==pc(a[p[e]]))e++;blocks.push_back({s,e});}uint32_t best=~0u;int aut=0;
 auto rec=[&](auto&&self,int b)->void{if(b==(int)blocks.size()){uint32_t q=encode(a,p);canon_perms++;if(q<best){best=q;aut=1;}else if(q==best)aut++;return;}auto[s,e]=blocks[b];sort(p.begin()+s,p.begin()+e);do{self(self,b+1);}while(next_permutation(p.begin()+s,p.begin()+e));};rec(rec,0);return{best,aut};}
struct Inc{array<int,8>v{};array<array<int,8>,8>e{};int triangles=0;};
Inc incidence(const A&a){Inc c;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)for(int k=j+1;k<8;k++){if(arc(a,i,j)==arc(a,j,k)&&arc(a,j,k)==arc(a,k,i)){c.triangles++;c.v[i]++;c.v[j]++;c.v[k]++;c.e[i][j]++;c.e[j][i]++;c.e[i][k]++;c.e[k][i]++;c.e[j][k]++;c.e[k][j]++;}}return c;}
Sig signature(const A&a,const Inc&c){Sig s{};for(int u=0;u<8;u++){array<unsigned,7> q{};int d=0;for(int v=0;v<8;v++)if(arc(a,u,v)){need(c.v[v]<13&&c.e[u][v]<8,"signature width overflow");q[d++]=1+c.v[v]*8+c.e[u][v];}sort(q.begin(),q.begin()+d);uint64_t z=d;for(int i=0;i<d;i++)z|=uint64_t(q[i])<<(3+7*i);s[u]=z;}sort(s.begin(),s.end());return s;}
struct SH{size_t operator()(const Sig&s)const{uint64_t h=1469598103934665603ull;for(auto x:s){h^=x;h*=1099511628211ull;h^=x>>32;h*=1099511628211ull;}return h;}};
int hamiltonian_dp(const A&a){int d[256][8]{};for(int v=0;v<8;v++)d[1<<v][v]=1;for(int s=1;s<256;s++)for(int v=0;v<8;v++)if(d[s][v]){unsigned ext=a[v]&(~s)&255;while(ext){int w=__builtin_ctz(ext);ext&=ext-1;d[s|(1<<w)][w]+=d[s][v];}}return accumulate(d[255],d[255]+8,0);}
struct Cl{int oldid,H,aut,tc,bm,SC,oldconv;uint32_t canon;A a;int conv=-1,node=-1;uint64_t seen=0,blue=0,closed=0;};
int main(int argc,char**argv){try{string dir=argc>1?argv[1]:".";ifstream f(dir+"/classes_input_n8.tsv");need(bool(f),"input absent");string l;getline(f,l);vector<Cl>C;while(getline(f,l)){auto v=split(l);need(v.size()==11,"TSV width");C.push_back({stoi(v[0]),stoi(v[2]),stoi(v[3]),stoi(v[4]),stoi(v[5]),stoi(v[6]),stoi(v[8]),(uint32_t)stoul(v[1])});C.back().a=decode(C.back().canon);}sort(C.begin(),C.end(),[](auto&a,auto&b){return a.canon<b.canon;});need(C.size()==6880,"rep count");map<uint32_t,int> cc;unordered_map<Sig,int,SH> lookup;lookup.reserve(10000);uint64_t mass=0;vector<int> old2new(C.size());for(int i=0;i<(int)C.size();i++){auto&c=C[i];auto[can,aut]=canonical(c.a);need(can==c.canon,"canonical label mismatch");need(aut==c.aut,"automorphism mismatch");need(cc.emplace(can,i).second,"duplicate class");need(40320%aut==0,"orbit integral");mass+=40320/aut;need(hamiltonian_dp(c.a)==c.H,"Hamiltonian count DP mismatch");need(c.H/c.aut==c.tc&&c.H%c.aut==0,"fiber formula");auto inc=incidence(c.a);need(lookup.emplace(signature(c.a,inc),i).second,"E signature collision");old2new[c.oldid]=i;}need(mass==(1ull<<28),"exhaustive labelled orbit mass");uint64_t source_perms=canon_perms;int nodecount=0,selfconverse=0;for(int i=0;i<(int)C.size();i++){auto[can,aut]=canonical(converse(C[i].a));need(cc.count(can),"missing converse");int j=cc[can];C[i].conv=j;need(j==old2new[C[i].oldconv],"input converse ID mismatch");need((j==i)==bool(C[i].SC),"SC mismatch");if(j>=i){C[i].node=C[j].node=nodecount++;}selfconverse+=(j==i);}for(int i=0;i<(int)C.size();i++)need(C[C[i].conv].conv==i,"converse involution");cerr<<"classes certified: "<<C.size()<<" E signatures, labelled orbit mass "<<mass<<", score permutations "<<source_perms<<"\n";
 vector<pair<int,int>> tiles;int index[8][8]{};for(int b=0;b<6;b++)for(int a=b+2;a<8;a++){index[b][a]=tiles.size();tiles.push_back({b,a});}need(tiles.size()==21,"free arcs");array<int,21> refl;int fixed=0;for(int k=0;k<21;k++){auto[b,a]=tiles[k];refl[k]=index[7-a][7-b];fixed+=refl[k]==k;}need(fixed==3,"reflection fixed tiles");const int N=1<<21,F=N-1;vector<uint16_t> cls(N);vector<uint8_t> start(N),end(N),blue(N);uint64_t bcount=0;ofstream masks(dir+"/independent_class_of_mask.u16le",ios::binary);for(int m=0;m<N;m++){A a{};for(int i=1;i<8;i++)a[i]|=1<<(i-1);for(int k=0;k<21;k++){auto[b,z]=tiles[k];if((m>>k)&1)a[b]|=1<<z;else a[z]|=1<<b;}auto inc=incidence(a);auto it=lookup.find(signature(a,inc));need(it!=lookup.end(),"unclassified mask "+to_string(m));int cid=it->second;cls[m]=cid;char bytes[2]={(char)(cid&255),(char)(cid>>8)};masks.write(bytes,2);C[cid].seen++;start[m]=pc(a[7]);end[m]=pc(a[0]);int r=0;for(int k=0;k<21;k++)r|=((m>>k)&1)<<refl[k];bool blu=r==m;blue[m]=blu;bcount+=blu;C[cid].blue+=blu;C[cid].closed+=(m>>index[0][7])&1;
 // Independent exact path-curl inverse: successive adjacent-root triples.
 array<int,8>x{};x[0]=0;x[7]=0;x[6]=1; // root7, path7->6->...->0
 auto sgn=[&](int u,int v){return u==v?0:(arc(a,u,v)?1:-1);};
 for(int i=6;i>0;i--){int kap=sgn(7,i)+sgn(i,i-1)+sgn(i-1,7);x[i-1]=x[i]+1-kap;}
 for(int i=0;i<7;i++)need(x[i]==sgn(7,i),"root star inverse");
 for(int i=0;i<7;i++)for(int j=i+1;j<7;j++){int kap=sgn(7,i)+sgn(i,j)+sgn(j,7);need(kap-x[i]+x[j]==sgn(i,j),"root curl inverse");}
 if((m&262143)==262143)cerr<<"classified "<<m+1<<" masks\n";
 }masks.close();need(bcount==4096,"blue tilings");for(auto&c:C){need(c.seen==(uint64_t)c.tc,"class fiber full mask mismatch");need(c.blue==(uint64_t)c.bm,"blue fiber full mask mismatch");need((c.seen&1)==1,"ordinary fiber odd");}
 using Port=array<int,7>;map<Port,uint64_t> ports;using Edge=array<int,3>;map<Edge,uint64_t> directed,undirected;uint64_t blacklines=0,bluelines=0;for(int m=0;m<N;m++){int r=0;for(int k=0;k<21;k++)r|=((m>>k)&1)<<refl[k];int t=F^m;need(cls[r]==C[cls[m]].conv,"reflection versus converse");need(blue[t]==blue[m],"complement color");need(start[t]==8-start[m]&&end[t]==6-end[m],"endpoint transport");need((F^r)!=m,"complement reflection fixed point");if(!((m>>index[0][7])&1))continue;int u=cls[m],v=cls[t],color=blue[m];ports[{u,start[m],end[m],color,v,start[t],end[t]}]++;int x=C[u].node,y=C[v].node;directed[{x,y,color}]++;if(x>y)swap(x,y);undirected[{x,y,color}]++;if(color)bluelines++;else blacklines++;}need(blacklines==1046528&&bluelines==2048,"line totals");uint64_t blackloops=0,blueloops=0,blackloop_support=0,blueloop_support=0;for(auto&[e,k]:directed)if(!e[2])need(k%2==0,"directed black twin parity including loops");for(auto&[e,k]:undirected){if(!e[2])need(k%2==0,"undirected black twin parity including loops");if(e[0]==e[1]){if(e[2]){blueloops+=k;blueloop_support++;}else{blackloops+=k;blackloop_support++;}}}
 ofstream cl(dir+"/independent_classes.tsv");cl<<"cid\tcanon\tconverse_cid\tnode_id\tH\tAut\tordinary_fiber\tblue_fiber\tclosed_fiber\n";for(int i=0;i<(int)C.size();i++){auto&c=C[i];cl<<i<<'\t'<<c.canon<<'\t'<<c.conv<<'\t'<<c.node<<'\t'<<c.H<<'\t'<<c.aut<<'\t'<<c.seen<<'\t'<<c.blue<<'\t'<<c.closed<<'\n';}
 ofstream pp(dir+"/independent_ports.tsv");pp<<"source_cid\tstart_score\tend_score\tblue\ttarget_cid\ttarget_start_score\ttarget_end_score\tcount\n";for(auto&[p,k]:ports){for(int q:p)pp<<q<<'\t';pp<<k<<'\n';}
 auto dump=[&](string name,const map<Edge,uint64_t>&es){ofstream o(dir+"/"+name);o<<"source_node\ttarget_node\tblue\tliteral_count\tfundamental_count\n";for(auto&[e,k]:es)o<<e[0]<<'\t'<<e[1]<<'\t'<<e[2]<<'\t'<<k<<'\t'<<(e[2]?k:k/2)<<'\n';};dump("independent_directed_edges.tsv",directed);dump("independent_undirected_edges.tsv",undirected);
 ofstream out(dir+"/independent_receipt.json");out<<"{\n\"status\":\"PASS\",\n\"universe\":2097152,\n\"labelled_orbit_mass\":"<<mass<<",\n\"classes\":"<<C.size()<<",\n\"score_canonical_permutations\":"<<source_perms<<",\n\"canonical_permutations_including_converses\":"<<canon_perms<<",\n\"E_signatures\":"<<lookup.size()<<",\n\"merged_nodes\":"<<nodecount<<",\n\"self_converse_classes\":"<<selfconverse<<",\n\"blue_tilings\":"<<bcount<<",\n\"literal_blue_lines\":"<<bluelines<<",\n\"literal_black_lines\":"<<blacklines<<",\n\"fundamental_black_lines\":"<<blacklines/2<<",\n\"port_entries\":"<<ports.size()<<",\n\"directed_edge_entries\":"<<directed.size()<<",\n\"undirected_edge_entries\":"<<undirected.size()<<",\n\"literal_black_loops\":"<<blackloops<<",\n\"blue_loops\":"<<blueloops<<",\n\"black_loop_support\":"<<blackloop_support<<",\n\"blue_loop_support\":"<<blueloop_support<<",\n\"explicit_checks\":"<<checks<<"\n}\n";cerr<<"PASS full independent replay, "<<checks<<" exact checks\n";
 }catch(exception&e){cerr<<"AUDIT FAILURE: "<<e.what()<<"\n";return 1;}}
