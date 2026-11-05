#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using ll=int64_t;
void require(bool b,const string&s){if(!b)throw runtime_error(s);}
using Graph=vector<int>;
Graph birth(Graph r,int mask){int n=r.size();for(int i=0;i<n;i++)if(!(mask>>i&1))r[i]|=1<<n;r.push_back(mask);return r;}
// Independent forward endpoint DP, applied to every candidate including rejects.
ll direct(const Graph&r){int n=r.size(),N=1<<n;vector<ll>d(N*n);for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;for(int m=1;m<N;m++)for(int i=0;i<n;i++){ll z=d[m*n+i];if(!z)continue;int todo=r[i]&~m;while(todo){int b=todo&-todo;todo-=b;int j=__builtin_ctz(b);d[(m|b)*n+j]+=z;}}ll z=0;for(int i=0;i<n;i++)z+=d[(N-1)*n+i];return z;}
bool strong(const Graph&r){int n=r.size();for(int i=0;i<n;i++){int seen=1<<i;while(true){int next=seen;for(int j=0;j<n;j++)if(seen>>j&1)next|=r[j];if(next==seen)break;seen=next;}if(seen!=(1<<n)-1)return false;}return true;}
// Reverse endpoint DP and disjoint path pairs compute every birth response.
struct Response {int n,N;vector<vector<ll>>e,s,c;vector<ll>h;Response(const Graph&r):n(r.size()),N(1<<n),e(N,vector<ll>(n)),s(e),c(n,vector<ll>(n)),h(N){h[0]=1;for(int i=0;i<n;i++)e[1<<i][i]=s[1<<i][i]=1;for(int m=1;m<N;m++)for(int i=0;i<n;i++)if(m>>i&1){int a=m^(1<<i);if(a)for(int j=0;j<n;j++)if(a>>j&1){if(r[j]>>i&1)e[m][i]+=e[a][j];else s[m][i]+=s[a][j];}h[m]+=e[m][i];}for(int m=1;m<N;m++)for(int i=0;i<n;i++)if(m>>i&1)for(int j=0;j<n;j++)if(!(m>>j&1))c[i][j]+=e[m][i]*s[(N-1)^m][j];}
ll value(int mask)const{ll z=0;for(int i=0;i<n;i++)z+=(mask>>i&1)?s[N-1][i]:e[N-1][i];for(int i=0;i<n;i++)if(!(mask>>i&1))for(int j=0;j<n;j++)if(mask>>j&1)z+=c[i][j];return z;}};
void rows(const Graph&r){cout<<"[";for(size_t i=0;i<r.size();i++){if(i)cout<<",";cout<<r[i];}cout<<"]";}
struct Seed{Graph r;int first,second;};
int main(int argc,char**argv){bool full=argc>1&&string(argv[1])=="--full";Graph D={2,12,9,17,7};require(direct(D)==13&&strong(D),"H13 seed");Response p(D);vector<Seed>seeds;ll source1=0,source2=0,source_checked=0;for(int a=0;a<32;a++){Graph q=birth(D,a);ll h=p.value(a);require(direct(q)==h,"source first-birth agreement");source_checked++;if(h==41)source1++;if(h>41)continue;Response t(q);for(int b=0;b<64;b++){auto u=birth(q,b);ll h2=t.value(b);require(direct(u)==h2,"source second-birth agreement");source_checked++;if(h2==41){require(strong(u),"every H41 source is strong");seeds.push_back({u,a,b});source2++;}}}
ll first_tested=0,first_retained=0,second_tested=0,second_retained=0,third_tested=0,hits1=0,hits2=0,hits3=0,final_direct=0;
cout<<"{\"mode\":\""<<(full?"full":"fast")<<"\",\"source_first_hits\":"<<source1<<",\"source_second_hits\":"<<source2<<",\"source_direct_checks\":"<<source_checked<<",\"seeds\":[";
for(size_t si=0;si<seeds.size();si++){auto seed=seeds[si];Response a(seed.r);ll local1=0,local2=0,local3=0;if(si)cout<<",";cout<<"{\"first_mask\":"<<seed.first<<",\"second_mask\":"<<seed.second<<",\"rows\":";rows(seed.r);
 for(int m=0;m<a.N;m++){first_tested++;auto q=birth(seed.r,m);ll h=a.value(m);require(direct(q)==h,"first-birth agreement");if(h==127)hits1++;if(h>127)continue;first_retained++;local1++;Response b(q);
 for(int m2=0;m2<b.N;m2++){second_tested++;auto u=birth(q,m2);ll h2=b.value(m2);require(direct(u)==h2,"second-birth agreement");if(h2==127)hits2++;if(h2>127)continue;second_retained++;local2++;Response c(u);
 for(int m3=0;m3<c.N;m3++){third_tested++;local3++;ll h3=c.value(m3);if(full){require(direct(birth(u,m3))==h3,"third-birth independent agreement");final_direct++;}if(h3==127)hits3++;}}}
 cout<<",\"retained_first\":"<<local1<<",\"retained_second\":"<<local2<<",\"third_masks\":"<<local3<<"}";}
require(source1==0&&source2==10,"complete source inventory");require(hits1==0&&hits2==0&&hits3==0,"no nested H127");
cout<<"],\"first_tested\":"<<first_tested<<",\"first_retained\":"<<first_retained<<",\"second_tested\":"<<second_tested<<",\"second_retained\":"<<second_retained<<",\"third_tested\":"<<third_tested<<",\"final_direct_checks\":"<<final_direct<<",\"target_hits\":["<<hits1<<","<<hits2<<","<<hits3<<"],\"status\":\"PASS\"}\n";
}
