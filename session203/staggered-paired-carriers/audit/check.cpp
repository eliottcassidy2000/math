// Independent append-at-end path DP and physical construction controls.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using G=std::vector<unsigned>;using I=std::int64_t;
// h,s,e,L,R,d,C,b,T
using P=std::array<I,9>;
int profiles=0;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
void valid(const G&a){int n=a.size();need(n>0&&n<=16,"order");for(int i=0;i<n;++i){need(a[i]<(1u<<n)&&!(a[i]>>i&1),"range/loop");for(int j=0;j<i;++j)need((a[i]>>j&1)+(a[j]>>i&1)==1,"orientation");}}
P profile(const G&a,int root,bool full=true){valid(a);++profiles;int n=a.size(),N=1<<n,F=N-1;std::vector<I>all(std::size_t(N)*n),marked(std::size_t(N)*n),h(N),s(N),e(N);h[0]=1;for(int i=0;i<n;++i)all[std::size_t(1<<i)*n+i]=1;marked[std::size_t(1<<root)*n+root]=1;
 for(int mask=1;mask<N;++mask)for(unsigned ends=mask;ends;ends&=ends-1){int last=__builtin_ctz(ends);std::size_t at=std::size_t(mask)*n+last;I v=all[at],w=marked[at];h[mask]+=v;s[mask]+=w;if(last==root)e[mask]=v;if(!v&&!w)continue;for(unsigned next=a[last]&~mask;next;next&=next-1){int u=__builtin_ctz(next);std::size_t to=std::size_t(mask|(1<<u))*n+u;all[to]+=v;marked[to]+=w;}}
 I L=0,R=0,C=0,b=0,T=0;for(int mask=0;mask<N;++mask){L+=e[mask]*h[F^mask];R+=s[mask]*h[F^mask];C+=h[mask]*h[F^mask];}int body=F^(1<<root);for(int mask=body;;mask=(mask-1)&body){b+=h[mask]*h[body^mask];if(!mask)break;}
 if(full){std::vector<I>bins(N);for(int mask=body;;mask=(mask-1)&body){for(int cut=mask;;cut=(cut-1)&mask){bins[mask]+=h[cut]*h[mask^cut];if(!cut)break;}T+=s[F^mask]*bins[mask];if(!mask)break;}}
 return {h[F],s[F],e[F],L,R,h[body],C,b,T};}
G move(G a,int root,char op){int n=a.size();a.resize(n+1);for(int i=0;i<n;++i){if(op!='O'&&i==root)a[n]|=1u<<i;else a[i]|=1u<<n;}valid(a);return a;}
G child(const G&a,int root){G c={2,8,3,5};c.resize(a.size()+3);std::vector<int>map(a.size());map[root]=0;int next=4;for(int i=0;i<(int)a.size();++i)if(i!=root)map[i]=next++;for(int i=0;i<(int)a.size();++i)for(int j=0;j<(int)a.size();++j)if(a[i]>>j&1)c[map[i]]|=1u<<map[j];for(int i=1;i<4;++i)for(int j=0;j<(int)a.size();++j)if(j!=root)c[i]|=1u<<map[j];valid(c);return c;}
P transfer(P p,char op,int order){auto[h,s,e,L,R,d,C,b,T]=p;if(op=='O')return{h,s,0,L,2*R,d,2*C,2*b,3*T};if(op=='Z')return{h+R-e,order==1?0:s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b};return{h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C};}
std::array<I,7> seven(P p){return{p[0],p[1],p[2],p[3],p[4],p[5],p[7]};}
I base(P p,int j=0){return p[0]+p[3]-p[1]+2*j;}
std::array<I,4> row(P p,int j){return{base(p,j),p[2],p[4]+p[7],p[5]};}
bool strong(const G&a){int n=a.size();for(int root=0;root<n;++root){unsigned seen=1u<<root,old;do{old=seen;for(int i=0;i<n;++i)if(seen>>i&1)seen|=a[i];}while(old!=seen);if(seen!=(1u<<n)-1)return false;}return true;}
I literal_f(const G&a,int root){int n=a.size();std::vector<int>p;for(int i=0;i<n;++i)if(i!=root)p.push_back(i);I f=0;do{std::vector<int>w{root};w.insert(w.end(),p.begin(),p.end());for(int cut=2;cut<n;++cut){bool ok=bool(a[root]>>w.back()&1);for(int j=0;j<n-1;++j)if(j!=cut-1&&!(a[w[j]]>>w[j+1]&1))ok=false;f+=ok;}}while(std::next_permutation(p.begin(),p.end()));return f;}
G graph(int n,unsigned code){G a(n);int bit=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){if(code>>bit++&1)a[i]|=1u<<j;else a[j]|=1u<<i;}return a;}
int main(){std::array<std::array<G,2>,4> carriers={{{G{2,12,17,101,75,23,167,63},G{2,12,17,165,43,71,31,119}},{G{2,44,49,133,11,89,31,119},G{2,52,25,35,105,197,15,95}},{G{2,148,97,7,45,75,27,125},G{2,148,41,67,45,75,23,125}},{G{4,9,50,85,163,75,151,47},G{4,9,50,85,99,139,167,31}}}};int roots[4]={2,1,0,2};int js[4][2]={{3,3},{4,3},{4,3},{5,5}};I X[4]={677,759,1547,821},coef[4]={1012,898,1582,1046},constant[4]={-313,-139,-27,-225};std::vector<std::vector<std::array<I,7>>>want={
 {{165,10,27,516,237,49,824},{175,11,27,507,249,49,812}},
 {{253,13,29,511,287,29,640},{253,11,29,511,287,29,640}},
 {{259,13,101,1293,267,109,1424},{259,11,101,1293,267,109,1424}},
 {{295,16,33,532,359,33,720},{293,17,33,535,359,33,720}}};int carrier_graphs=0,readouts=0;std::array<std::array<P,2>,4>seed;
 for(int i=0;i<4;++i){for(int r=0;r<2;++r){seed[i][r]=profile(carriers[i][r],roots[i]);need(seven(seed[i][r])==want[i][r],"eight literal seed profiles");need(base(seed[i][r],js[i][r])==X[i]&&seed[i][r][1]==r+1+3*js[i][r],"adjusted seed base/endpoint");}need(row(seed[i][0],js[i][0])==row(seed[i][1],js[i][1]),"coupled seed row");
 for(int len=0;len<=5;++len)for(int code=0;code<(1<<len);++code){P actual[2];for(int r=0;r<2;++r){G a=carriers[i][r];P p=seed[i][r];for(int j=0;j<len;++j){char op=(code>>j&1)?'Z':'O';p=transfer(p,op,a.size());a=move(a,roots[i],op);}actual[r]=profile(a,roots[i],false);need(seven(actual[r])==seven(p),"physical O/Z word");need(base(actual[r],js[i][r])>=137&&(base(actual[r],js[i][r])&1),"odd base floor");++carrier_graphs;if(len<=3){P out=profile(child(a,roots[i]),0,false);need(out[0]==3*base(actual[r],js[i][r])+2*(r+1),"both actual residue readouts");++readouts;}}need(row(actual[0],js[i][0])==row(actual[1],js[i][1]),"paired common-word rows");}
 for(int r=0;r<2;++r){P p=seed[i][r];for(int k=1;k<=20;++k){p=transfer(p,'Z',8+k-1);need(base(p,js[i][r])==coef[i]*(1LL<<k)+constant[i],"pure-Z orbit formula");}need((coef[i]+constant[i]==X[i])==(i==1||i==3),"k=0 exceptional formulas");}}
 need(seed[1][0][6]==4384&&seed[1][1][6]==4552&&seed[1][0][8]==2713&&seed[1][1][8]==2797,"759 C/T coordinates");I hostile[2];for(int r=0;r<2;++r){G a=move(carriers[1][r],roots[1],'J');P p=profile(a,8);hostile[r]=base(p,js[1][r]);need(hostile[r]==(r?4547:4379),"J breaks paired row");}
 int roots_checked=0,face_graphs=0,cover_checks=0;
 for(int n=1;n<=4;++n)for(unsigned code=0;code<(1u<<(n*(n-1)/2));++code){G a=graph(n,code);for(int r=0;r<n;++r){P p=profile(a,r);need(p[2]>=0&&p[5]>=p[2]&&p[4]-p[0]+p[2]-p[5]>=0,"three nonnegative faces");need(literal_f(a,r)==p[4]-p[0]+p[2]-p[5],"direct remaining-cover interpretation");++cover_checks;++roots_checked;for(char op:std::string("OZJ")){G b=move(a,r,op);int rr=op=='J'?n:r;P q=profile(b,rr);need(q==transfer(p,op,n),"all-input nine-coordinate move");if(op=='O')need(q[2]==0,"O face");if(op=='Z')need(q[5]-q[2]==0,"Z face");if(op=='J')need(q[4]-q[0]+q[2]-q[5]==0,"J face");need(literal_f(b,rr)==q[4]-q[0]+q[2]-q[5],"literal child cover interpretation");++cover_checks;++face_graphs;}}}
 G u{2,0};int root=0,prev=-1;std::vector<std::array<I,6>>family;I n3f=-1;
 for(int n=3;n<=14;++n){int oldroot=root;u=move(u,root,'J');prev=oldroot;root=n-1;need(strong(u),"strong U_n");P ur=profile(u,root),up=profile(u,prev);need(ur[1]==1&&up[1]==1&&__builtin_popcount(u[root])==1&&__builtin_popcount(u[prev])==1,"two unit starts/outdegrees");G w=u;w.resize(n+1);w[n]=(1u<<root)|(1u<<prev);for(int j=0;j<n;++j)if(j!=root&&j!=prev)w[j]|=1u<<n;P q=profile(w,n);I e=q[2],gap=q[5]-e,f=q[4]-q[0]+e-q[5];need(q[1]==2&&strong(w),"start-two strong W_n");if(n==3){n3f=f;need(f>=0,"n=3 boundary nonnegative");continue;}need(e>0&&gap>0&&f>0&&q[5]>=n,"strict interior/unbounded d");need(gap==ur[2]+up[2],"deletion-endpoint formula");family.push_back({n+1,q[0],q[1],e,gap,f});}
 std::cout<<"{\"status\":\"PASS\",\"carrier_word_graphs\":"<<carrier_graphs<<",\"literal_four_core_readouts\":"<<readouts<<",\"all_rooted_inputs_through4\":"<<roots_checked<<",\"literal_O_Z_J_children\":"<<face_graphs<<",\"direct_f_cover_checks\":"<<cover_checks<<",\"J_hostile_adjusted_bases\":["<<hostile[0]<<","<<hostile[1]<<"],\"n3_boundary_f\":"<<n3f<<",\"profile_computations\":"<<profiles<<",\"W_family_order_h_s_e_d_minus_e_f\":[";for(int i=0;i<(int)family.size();++i){std::cout<<(i?",":"")<<"[";for(int j=0;j<6;++j)std::cout<<(j?",":"")<<family[i][j];std::cout<<"]";}std::cout<<"],\"pair_completeness_claim\":false,\"F_8_137\":\"OPEN\"}\n";
}
