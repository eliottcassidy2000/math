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
I three(const Paths&p,unsigned S){I ans=0;for(unsigned A=S;;A=(A-1)&S){ans+=p.H[A]*p.two(S^A);if(!A)break;}return ans;}
I endthree(const Paths&p,unsigned S,int r){I ans=0;for(unsigned A=S;;A=(A-1)&S){ans+=p.ends[A*p.n+r]*p.two(S^A);if(!A)break;}return ans;}
std::vector<unsigned> components(const G&a){int n=a.size();std::vector<unsigned>out;unsigned used=0;for(int r=0;r<n;r++)if(!(used>>r&1)){unsigned R=0;for(int x=0;x<n;x++){auto reaches=[&](int src,int dst){unsigned seen=1u<<src;std::vector<int>q{src};for(unsigned k=0;k<q.size();k++)for(int v=0;v<n;v++)if((a[q[k]]>>v&1)&&!(seen>>v&1)){seen|=1u<<v;q.push_back(v);}return bool(seen>>dst&1);};if(reaches(r,x)&&reaches(x,r))R|=1u<<x;}out.push_back(R);used|=R;}std::sort(out.begin(),out.end(),[&](unsigned A,unsigned B){return a[__builtin_ctz(A)]>>__builtin_ctz(B)&1;});return out;}
G join(const G&a,const G&b){G g=a;int n=a.size(),m=b.size();for(auto &r:g)r|=((1u<<m)-1)<<n;for(auto r:b)g.push_back(r<<n);return g;}
int main(){G triangle{2,4,1};I cases=0;int maxn=0;std::set<I>counts;for(int sz:{1,3})for(int k=0;k<=2;k++)for(int l=0;l<=2;l++)for(int place=0;place<4;place++){
 G P=triangle,Q=triangle,M=sz==1?G{0}:triangle;for(int i=0;i<k;i++)P=place&1?join(P,G{0}):join(G{0},P);for(int i=0;i<l;i++)Q=place&2?join(Q,G{0}):join(G{0},Q);G T=join(join(P,M),Q);int x=P.size();Paths p(T),pm(M),pc(triangle);G u=cyclic(T);u[1]&=~(1u<<(x+2));u[x+2]|=2;Paths pu(u);I alpha=3*pc.two(pc.F)*pm.two(pm.F)*pc.two(pc.F),beta=2*three(pc,pc.F)*endthree(pm,pm.F,0)*pc.two(pc.F),gamma=3*pc.H[pc.F]*pm.R(pm.F,0)*pc.two(pc.F),pow3=1;for(int i=0;i<k;i++)pow3*=3;I expr=(1ll<<l)*(alpha*(1ll<<k)+beta*pow3-gamma);auto q=pu.p(0);need(q[7]+p.H[p.F]==expr&&q[8]==q[1]&&q[1]==p.H[p.F]&&strong(u),"nonempty-core template");need(alpha>=2*gamma&&beta>0,"template no cancellation");++cases;maxn=std::max(maxn,int(u.size()));counts.insert(q[7]);
 }std::cout<<"{\"status\":\"PASS\",\"literal_template_sources\":"<<cases<<",\"both_unmarked_cores_nonempty\":true,\"max_order\":"<<maxn<<",\"distinct_B\":"<<counts.size()<<"}\n";
}
