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
int main(){I graphs=0,roots=0,pairs=0,terminalStrong=0,equalities=0,sharp=0,flips=0,coefficients=0,singleMarks=0,emptyCores=0;std::map<int,I>defectByStart;
 for(int n=1;n<=6;n++)for(unsigned code=0;code<(1u<<(n*(n-1)/2));code++){
 G a=graph(n,code);Paths p(a);graphs++;bool isstrong=strong(a);
 for(int x=0;x<n;x++)for(int y=x+1;y<n;y++){
 I count=0;for(unsigned S=1;S<unsigned(p.F);S++)count+=p.starts[S*n+x]*p.ends[(p.F^S)*n+y]+p.starts[S*n+y]*p.ends[(p.F^S)*n+x];need(count>0,"two-terminal cover");pairs++;
 if(isstrong){need(p.ends[(p.F^(1u<<x))*n+y]+p.ends[(p.F^(1u<<y))*n+x]>0,"strong deletion terminal alternative");terminalStrong++;}}
 auto cc=components(a);for(int v=0;v<n;v++){
 auto[h,s,e,L,R,d,b,B,D]=p.p(v);unsigned A=p.F^(1u<<v),incoming=A&~a[v];I sum=0,twopath=0;for(int x=0;x<n;x++)if(incoming>>x&1)sum+=p.starts[A*n+x];
 for(unsigned S=A;;S=(S-1)&A){I ends=0,starts=0;for(int x=0;x<n;x++)if(incoming>>x&1){ends+=p.ends[S*n+x];starts+=p.starts[(A^S)*n+x];}twopath+=ends*starts;if(!S)break;}
 int deg=__builtin_popcount(incoming);need(L-h==sum+twopath&&twopath>=deg*(deg-1)/2,"exact bad cover split");need((L==h)==(deg==0||(deg==1&&sum==0)),"cover equality iff");if(deg<=1)need(twopath==0,"two-path term vanishes");equalities+=(L==h);
 if(n>1&&s>0&&D==s){sharp++;need(__builtin_popcount(a[v])==1,"unique root successor");int w=__builtin_ctz(a[v]);unsigned T=A^(1u<<w),I=A&~a[w]&~(1u<<w);need(T&&__builtin_popcount(I)<=1&&p.H[T]==s,"residual form and count");if(I){int x=__builtin_ctz(I);need(!p.starts[T*n+x]&&!p.ends[T*n+x],"interior mark");}need(isstrong,"strong equality graph");defectByStart[s]++;}roots++;}
 auto test=[&](int x){G u=cyclic(a);if(x>=0){u[1]&=~(1u<<(x+2));u[x+2]|=2;}Paths c(u);auto q=c.p(0);I t=p.H[p.F],C=p.two(p.F);need(q[1]==t&&q[8]==t&&strong(u),"physical source floor and strongness");if(x<0){need(q[7]+t==3*C,"baseline count");return;}
 I l=p.L(p.F,x),r=p.R(p.F,x),z=endthree(p,p.F,x);need(q[0]==t+C+l-r&&q[2]==t+l&&q[5]==t+l&&q[6]==2*C-2*r+2*z&&q[7]+t==3*C-3*r+2*z,"flip literal graph formulas");flips++;
 int mi=-1;for(unsigned i=0;i<cc.size();i++)if(cc[i]>>x&1)mi=i;need(mi>0&&mi+1<int(cc.size()),"interior SCC");unsigned P=0,M=cc[mi],Q=0,PC=0,QC=0;int k=0,lcount=0;for(int i=0;i<int(cc.size());i++)if(i<mi){P|=cc[i];if(__builtin_popcount(cc[i])==1)k++;else PC|=cc[i];}else if(i>mi){Q|=cc[i];if(__builtin_popcount(cc[i])==1)lcount++;else QC|=cc[i];}
 need(C==p.two(P)*p.two(M)*p.two(Q)&&r==p.H[P]*p.R(M,x)*p.two(Q)&&z==three(p,P)*endthree(p,M,x)*p.two(Q),"component cover factorization");
 I alpha=3*p.two(PC)*p.two(M)*p.two(QC),beta=2*three(p,PC)*endthree(p,M,x)*p.two(QC),gamma=3*p.H[PC]*p.R(M,x)*p.two(QC);I threek=1;for(int i=0;i<k;i++)threek*=3;I expr=(1ll<<lcount)*(alpha*(1ll<<k)+beta*threek-gamma);
 need(alpha>=2*gamma&&beta>0&&gamma>0,"positive constants and no cancellation");need(q[7]+t==expr&&expr>=beta*(1ll<<lcount)*threek,"exact two-index expression/lower bound");coefficients++;singleMarks+=(__builtin_popcount(M)==1);emptyCores+=(!PC||!QC);
 };test(-1);for(int x=0;x<n;x++)if(!p.starts[p.F*n+x]&&!p.ends[p.F*n+x])test(x);
 }
 std::cout<<"{\"status\":\"PASS\",\"literal_graphs\":"<<graphs<<",\"rooted_inputs\":"<<roots<<",\"terminal_pairs\":"<<pairs<<",\"strong_terminal_alternatives\":"<<terminalStrong<<",\"cover_equalities\":"<<equalities<<",\"defect_equalities\":"<<sharp<<",\"baseline_sources\":"<<graphs<<",\"flipped_sources\":"<<flips<<",\"exact_two_parameter_checks\":"<<coefficients<<",\"marked_singleton_cases\":"<<singleMarks<<",\"empty_core_cases\":"<<emptyCores<<"}\n";
}
