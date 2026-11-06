#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using I=std::int64_t;using G=std::vector<unsigned>;
void need(bool q,const char*s){if(!q)throw std::runtime_error(s);}
struct Paths{G a;int n,N,F;std::vector<I>H,starts,ends;
 Paths(G aa):a(aa),n(a.size()),N(1<<n),F(N-1),H(N),starts(N*n),ends(N*n){for(int v=0;v<n;++v){need(a[v]<unsigned(N)&&!(a[v]>>v&1),"row");for(int u=0;u<v;++u)need((a[v]>>u&1)+(a[u]>>v&1)==1,"orientation");}H[0]=1;auto walk=[&](auto&&self,int first,int last,unsigned S)->void{++H[S];++starts[S*n+first];++ends[S*n+last];unsigned next=a[last]&~S;while(next){unsigned bit=next&-next;next-=bit;self(self,first,__builtin_ctz(bit),S|bit);}};for(int v=0;v<n;++v)walk(walk,v,v,1u<<v);}
 I L(unsigned S,int r)const{I x=0;for(unsigned A=S;;A=(A-1)&S){x+=ends[A*n+r]*H[S^A];if(!A)break;}return x;}
 I R(unsigned S,int r)const{I x=0;for(unsigned A=S;;A=(A-1)&S){x+=starts[A*n+r]*H[S^A];if(!A)break;}return x;}
 I two(unsigned S)const{I x=0;for(unsigned A=S;;A=(A-1)&S){x+=H[A]*H[S^A];if(!A)break;}return x;}
 std::array<I,9>p(int r)const{unsigned A=F^(1u<<r);I h=H[F],s=starts[F*n+r],e=ends[F*n+r],l=L(F,r),rr=R(F,r),d=H[A],b=two(A),B=h+l-s,D=rr+b-e-B;return{h,s,e,l,rr,d,b,B,D};}
};
bool strong(const G&a){int n=a.size();for(int r=0;r<n;++r){std::vector<int>q{r};unsigned seen=1u<<r;for(int i=0;i<(int)q.size();++i)for(int v=0;v<n;++v)if((a[q[i]]>>v&1)&&!(seen>>v&1)){seen|=1u<<v;q.push_back(v);}if(seen!=(1u<<n)-1)return false;}return true;}
G z(G a,int r){int n=a.size();for(int v=0;v<n;++v)if(v!=r)a[v]|=1u<<n;a.push_back(1u<<r);return a;}
G cyclic(const G&t){G a(t.size()+2);a[0]=2;for(int v=0;v<(int)t.size();++v){a[1]|=1u<<(v+2);a[v+2]=(t[v]<<2)|1u;}return a;}
G graph(int n,unsigned code){G a(n);int bit=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++bit)if(code>>bit&1)a[i]|=1u<<j;else a[j]|=1u<<i;return a;}
int main(){I graphs=0,marked=0,positive=0,eq=0,children=0;std::map<int,I>eqorder;std::set<I>values;for(int n=1;n<=6;++n)for(unsigned bits=0;bits<(1u<<(n*(n-1)/2));++bits){G a=graph(n,bits);Paths p(a);++graphs;for(int r=0;r<n;++r){++marked;auto[h,s,e,L,R,d,b,B,D]=p.p(r);need(s<=d&&L>=h+d-s,"union/deletion inequalities");if(L==h)need(s==d,"correct L equality implication");if(n==1){need(D==0&&s==1,"singleton floor exception");Paths singleton_child(z(a,r));need(singleton_child.p(r)[7]==3&&singleton_child.p(r)[7]-2*B==1,"singleton physical defect exception");continue;}if(!s)continue;++positive;unsigned A=p.F^(1u<<r);I sum=0;for(int v=0;v<n;++v)if(a[r]>>v&1)sum+=p.L(A,v);I f=R-h+e-d,pp=d-e;need(sum==b-L&&f>=0&&pp>=0&&sum>=d,"exact positive pieces");need(D-s==f+2*pp+sum-d&&D>=s,"defect floor");bool criterion=false;int v=-1;if(__builtin_popcount(a[r])==1){v=__builtin_ctz(a[r]);criterion=p.L(A,v)==d&&p.ends[A*n+v]==0;need(s==p.starts[A*n+v]&&pp==p.ends[A*n+v]&&f==0,"unique-outneighbor identities");}need((D==s)==criterion,"iff graph test");if(D==s){++eq;++eqorder[n];values.insert(s);need(n>=3&&v>=0&&p.H[A^(1u<<v)]==s,"nonempty residual count");need(strong(a),"every equality source is strong");}if(n<=5){Paths child(z(a,r));auto out=child.p(r);need(out[1]==s&&out[7]==2*B+D,"physical Z count");++children;}}
 G u=cyclic(a);need(strong(u),"cyclic strong output");Paths pu(u);I t=p.H[p.F],c=p.two(p.F);need(pu.p(0)==std::array<I,9>{t+c,t,t,2*c-t,t+c,t,2*c,3*c-t,t},"all cyclic coordinates");}
 Paths hostile(G{14,4,8,2});need(hostile.H[hostile.F]==3&&hostile.L(hostile.F,0)==3&&hostile.starts[hostile.F*4]==3,"false start-one implication hostile");std::cout<<"{\"status\":\"PASS\",\"literal_tournaments\":"<<graphs<<",\"singleton_formula_D\":0,\"singleton_physical_difference\":1,\"marked_inputs\":"<<marked<<",\"nontrivial_positive_start\":"<<positive<<",\"equality_sources\":"<<eq<<",\"all_equality_sources_strong\":true,\"physical_Z_children\":"<<children<<",\"strong_cyclic_sources\":"<<graphs<<",\"equality_by_order\":{";bool comma=false;for(auto[n,c]:eqorder){std::cout<<(comma?",":"")<<'"'<<n<<"\":"<<c;comma=true;}std::cout<<"},\"equality_start_values\":[";comma=false;for(I v:values){std::cout<<(comma?",":"")<<v;comma=true;}std::cout<<"]}\n";}
