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
using Poly=std::map<std::vector<int>,I>;
Poly skeletons(const G&g,int c,int r){Poly poly;std::vector<int>q(r);auto rec=[&](auto&&self,unsigned used,int last,int lastRun)->void{if(used==(1u<<c)-1)poly[q]++;for(int i=0;i<c;i++)if(!(used>>i&1)&&(last<0||(g[last]>>i&1)))self(self,used|(1u<<i),i,-1);for(int j=lastRun+1;j<r;j++)if(last<0||(g[last]>>(c+j)&1)){q[j]++;self(self,used,c+j,j);q[j]--;}};rec(rec,0,-1,-1);return poly;}
I onto(int n,int q){if(q==0)return n==0;I ans=0,choose=1;for(int j=0;j<=q;j++){I pw=1;for(int k=0;k<n;k++)pw*=j;ans+=((q-j)&1?-1:1)*choose*pw;if(j<q)choose=choose*(q-j)/(j+1);}return ans;}
I eval(const Poly&p,const std::vector<int>&n){I ans=0;for(const auto&[q,m]:p){I term=m;for(unsigned i=0;i<n.size();i++)term*=onto(n[i],q[i]);ans+=term;}return ans;}
G expand(const G&g,int c,const std::vector<int>&lens){std::vector<int>type;for(int i=0;i<c;i++)type.push_back(i);for(unsigned i=0;i<lens.size();i++)for(int j=0;j<lens[i];j++)type.push_back(c+i);G out(type.size());for(unsigned i=0;i<type.size();i++)for(unsigned j=i+1;j<type.size();j++){if(type[i]==type[j]||(g[type[i]]>>type[j]&1))out[i]|=1u<<j;else out[j]|=1u<<i;}return out;}
// Direct Hamilton-path enumeration records its actual sequence of core vertices/run segments.
Poly pathHist(const G&g,int c,const std::vector<int>&lens){int n=g.size();std::vector<int>type(n);for(int i=0;i<c;i++)type[i]=i;int v=c;for(unsigned i=0;i<lens.size();i++)for(int j=0;j<lens[i];j++)type[v++]=c+i;Poly out;std::vector<int>q(lens.size());auto rec=[&](auto&&self,unsigned used,int last)->void{if(used==(1u<<n)-1){out[q]++;return;}for(int j=0;j<n;j++)if(!(used>>j&1)&&(last<0||(g[last]>>j&1))){bool newRun=type[j]>=c&&(last<0||type[last]!=type[j]);if(newRun)q[type[j]-c]++;self(self,used|(1u<<j),j);if(newRun)q[type[j]-c]--;}};rec(rec,0,-1);return out;}
int main(){I patterns=0,literal=0,histChecks=0,normalizations=0,normInputs=0;int maxorder=0;for(int c=0;c<=3;c++)for(int r=0;r<=2;r++)for(unsigned code=0;code<(1u<<(c*(c-1)/2));code++)for(unsigned cross=0;cross<(1u<<(c*r));cross++){
 G pat=graph(c,code);pat.resize(c+r);for(int i=0;i<r;i++)for(int j=i+1;j<r;j++)pat[c+i]|=1u<<(c+j);for(int i=0;i<c;i++)for(int j=0;j<r;j++)if(cross>>(i*r+j)&1)pat[i]|=1u<<(c+j);else pat[c+j]|=1u<<i;Poly sk=skeletons(pat,c,r);patterns++;for(const auto&[q,m]:sk){need(m>0,"positive multiplicity");for(int x:q)need(x<=c+1,"occurrence cap");}
 std::vector<int>lens(r);auto literals=[&](auto&&self,int pos)->void{if(pos<r){for(lens[pos]=0;lens[pos]<=3;lens[pos]++)self(self,pos+1);return;}G g=expand(pat,c,lens);Paths p(g);I total=eval(sk,lens);need(total==p.H[p.F],"literal positive skeleton count");Poly observed=pathHist(g,c,lens),expected;for(const auto&[q,m]:sk){I x=m;for(int i=0;i<r;i++)x*=onto(lens[i],q[i]);if(x)expected[q]=x;}need(expected==observed,"literal path-run histograms");literal++;histChecks+=observed.size();maxorder=std::max(maxorder,int(g.size()));};literals(literals,0);
 auto normalized=[&](auto&&self,int pos)->void{if(pos<r){for(lens[pos]=0;lens[pos]<=12;lens[pos]++)self(self,pos+1);return;}I X=eval(sk,lens);need(X>=1,"tournament HP positivity");int lg=0;for(I y=X;y>=2;y/=2)lg++;int N=lg+c+1;auto shorter=lens;bool changed=false;for(int &x:shorter)if(x>N){x=N;changed=true;}need(eval(sk,shorter)==X,"simultaneous normalization");normInputs++;normalizations+=changed;};normalized(normalized,0);
 }
 std::cout<<"{\"status\":\"PASS\",\"patterns\":"<<patterns<<",\"literal_expansions\":"<<literal<<",\"path_run_histogram_entries\":"<<histChecks<<",\"maximum_literal_order\":"<<maxorder<<",\"normalization_inputs\":"<<normInputs<<",\"strict_normalizations\":"<<normalizations<<"}\n";
}
