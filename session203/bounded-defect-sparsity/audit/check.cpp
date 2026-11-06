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
int main(){I positive=0,multi=0,unique=0,edges=0,patterns=0,runs=0;std::map<I,I>caps;
 for(int n=2;n<=6;n++)for(unsigned code=0;code<(1u<<(n*(n-1)/2));code++){
 G a=graph(n,code);Paths p(a);for(int v=0;v<n;v++){
 auto[h,s,e,L,R,d,b,B,D]=p.p(v);if(!s)continue;positive++;need(D>=s&&D>0,"positive defect");int k=__builtin_popcount(a[v]);need(k>0,"positive start implies successor");unsigned T=p.F^(1u<<v),marks=a[v];
 if(k>=2){need(D>=s+(k-1)*d&&d<=D&&k<=D,"multiple outgoing bound");multi++;}else{if(n==2){need(B==1&&D==3,"rooted edge exception");edges++;continue;}int w=__builtin_ctz(a[v]);T^=1u<<w;marks=T&~a[w];I starts=0,ends=0,covers=0;for(int x=0;x<n;x++)if(marks>>x&1){starts+=p.starts[T*n+x];ends+=p.ends[T*n+x];}for(unsigned A=T;;A=(A-1)&T){I l=0,r=0;for(int x=0;x<n;x++)if(marks>>x&1){l+=p.ends[A*n+x];r+=p.starts[(T^A)*n+x];}covers+=l*r;if(!A)break;}int m=__builtin_popcount(marks);need(D==p.H[T]+2*ends+covers&&s==p.H[T]-starts&&covers>=m*(m-1)/2&&p.H[T]<=D&&m<=D,"exact unique successor reduction");unique++;}
 // Independent SCC extraction on the residual induced graph, preserving original labels.
 std::vector<int>labels;for(int x=0;x<n;x++)if(T>>x&1)labels.push_back(x);G t(labels.size());for(unsigned i=0;i<labels.size();i++)for(unsigned j=0;j<labels.size();j++)if(a[labels[i]]>>labels[j]&1)t[i]|=1u<<j;
 auto cc=components(t);unsigned keep=0;I nontrivialCountProduct=1;int nontrivialVertices=0,markedSingletons=0;std::vector<unsigned>rr;unsigned running=0;
 for(auto C:cc){unsigned old=0;for(unsigned i=0;i<labels.size();i++)if(C>>i&1)old|=1u<<labels[i];bool retained=__builtin_popcount(C)>1||(old&marks);if(__builtin_popcount(C)>1){nontrivialVertices+=__builtin_popcount(C);nontrivialCountProduct*=p.H[old];}else if(old&marks)markedSingletons++;if(retained){if(running){rr.push_back(running);running=0;}keep|=old;}else running|=old;}
 if(running)rr.push_back(running);need(nontrivialCountProduct==p.H[T]&&nontrivialVertices<=D&&markedSingletons<=D&&int(rr.size())<=2*D+1,"bounded fixed core");need(__builtin_popcount(keep)+n-__builtin_popcount(T)+1<=2*D+3,"E core size");
 // Physically add the edge-extension vertex z.
 G ext=a;for(int x=0;x<n;x++)if(x==v)ext[x]|=1u<<n;ext.push_back(((1u<<n)-1)^(1u<<v));Paths pe(ext);need(pe.H[pe.F]==B,"physical edge extension B");
 unsigned unionRuns=0;for(unsigned i=0;i<rr.size();i++){unsigned C=rr[i];unionRuns|=C;for(int x=0;x<=n;x++)if(!(C>>x&1)){unsigned out=ext[x]&C;need(!out||out==C,"run homogeneity including E vertex");}for(unsigned j=i+1;j<rr.size();j++)for(int x=0;x<n;x++)if(C>>x&1)need((ext[x]&rr[j])==rr[j],"ordered runs");for(int x=0;x<n;x++)if(C>>x&1)for(int y=0;y<n;y++)if((C>>y&1)&&(ext[x]>>y&1))for(int z=0;z<n;z++)if((C>>z&1)&&(ext[y]>>z&1))need(ext[x]>>z&1,"transitive run");}
 need((unionRuns|keep)==T&&!(unionRuns&keep),"partition");patterns++;runs+=rr.size();for(int K:{1,2,3,4,5,8,16,26})if(D<=K)caps[K]++;
 }}
 std::cout<<"{\"status\":\"PASS\",\"positive_marked_graphs\":"<<positive<<",\"multiple_successor\":"<<multi<<",\"unique_successor_nonedge\":"<<unique<<",\"rooted_edges\":"<<edges<<",\"extracted_patterns\":"<<patterns<<",\"homogeneous_runs\":"<<runs<<",\"pattern_counts_by_cap\":{";bool comma=false;for(auto[k,n]:caps){std::cout<<(comma?",":"")<<'"'<<k<<"\":"<<n;comma=true;}std::cout<<"}}\n";
}
