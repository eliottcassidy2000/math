#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
using Word=vector<int>;
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
struct Tour {int n; vector<int> a; vector<vector<Word>> words;
 Tour(int n_,uint64_t code):n(n_),a(n),words(1<<n){int k=0; for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++) if((code>>k)&1)a[i]|=1<<j;else a[j]|=1<<i; Word p;words[0].push_back(p); grow(p,0,-1);}
 void grow(Word&p,int mask,int last){int avail=((1<<n)-1)^mask;if(last>=0)avail&=a[last]; while(avail){int v=__builtin_ctz(avail);avail&=avail-1;p.push_back(v);words[mask|1<<v].push_back(p);grow(p,mask|1<<v,v);p.pop_back();}}
 bool edge(int u,int v)const{return(a[u]>>v)&1;}
 bool strong()const{for(int v=0;v<n;v++){int r=1<<v,old=0;while(r!=old){old=r;for(int i=0;i<n;i++)if((r>>i)&1)r|=a[i];}if(r!=(1<<n)-1)return false;}return true;}
};
bool pow2(int x){return x>0&&(x&(x-1))==0;}
bool defect3(int e){return e==1||(e%3==0&&pow2(e/3))||((e+2)%3==0&&(e+2)/3>=2&&pow2((e+2)/3))||((e+2)%5==0&&(e+2)/5>=2&&pow2((e+2)/5));}
struct Cover{int h,k,l;};
Cover covers(const Tour&t,int U,int O){int k=0,l=0;for(int P=U;;P=(P-1)&U){int Q=U^P;if(Q){int starts=0,ends=0;for(const auto&q:t.words[Q])starts+=(O>>q[0])&1;for(const auto&p:t.words[P])if(!p.empty())ends+=(O>>p.back())&1;k+=starts*ends;l+=starts*t.words[P].size();}if(!P)break;}return{(int)t.words[U].size(),k,l};}
vector<int>reach(const Tour&t,int U){vector<int>r(t.n);for(int i=0;i<t.n;i++)if((U>>i)&1)r[i]=(t.a[i]&U)|(1<<i);for(int k=0;k<t.n;k++)for(int i=0;i<t.n;i++)if((r[i]>>k)&1)r[i]|=r[k];return r;}
int main(){uint64_t roots=0,positive=0,strongs=0,compressed=0,dthree=0;set<int>spect;vector<array<uint64_t,7>>by;for(int n=1;n<=6;n++){uint64_t rr=0,pp=0,ss=0,cc=0,dd=0;for(uint64_t code=0;code<(1ULL<<(n*(n-1)/2));code++){Tour t(n,code);int full=(1<<n)-1,h=t.words[full].size();for(int a=0;a<n;a++)if(__builtin_popcount((unsigned)t.a[a])==1){rr++;int r=__builtin_ctz(t.a[a]),U=full^(1<<a)^(1<<r);if(!U)continue;int O=t.a[r]&U,e=0;for(const auto&p:t.words[full])e+=p.back()==a;auto [hu,k,l]=covers(t,U,O);need(e==l-k,"endpoint identity");int d=h-2*e;need(d==hu+k,"defect identity");need(d>=hu+__builtin_popcount((unsigned)O)*(__builtin_popcount((unsigned)O)-1),"capacity floor");need((k==0)==(__builtin_popcount((unsigned)O)<=1),"zero capacity");if(!e)continue;pp++;if(d==3){dd++;spect.insert(e);need(defect3(e),"defect three classification");}if(!t.strong())continue;ss++;auto re=reach(t,U);int removed=0;for(int z=0;z<n;z++)if((U>>z)&1){bool single=true,later=true;for(int y=0;y<n;y++)if(y!=z&&((re[z]>>y)&1)&&((re[y]>>z)&1))single=false;for(int y=0;y<n;y++)if(((O>>y)&1)&&!t.edge(y,z))later=false;if(single&&later)removed|=1<<z;}int s=__builtin_popcount((unsigned)removed);cc+=s>0;int U0=U^removed;auto c=covers(t,U0,O);need(c.h==hu&&c.k==k&&l==(1<<s)*c.l,"compression");int full0=full^removed;auto re0=reach(t,full0);for(int z=0;z<n;z++)if((full0>>z)&1)need(re0[z]==full0,"compressed strongness");int log=0,b=max(1,d-1);while(b>=2){b/=2;log++;}need(__builtin_popcount((unsigned)full0)<=d+log+4,"skeleton bound");}}
 by.push_back({(uint64_t)n,1ULL<<(n*(n-1)/2),rr,pp,ss,cc,dd});roots+=rr;positive+=pp;strongs+=ss;compressed+=cc;dthree+=dd;}
 cout<<"{\"status\":\"PASS\",\"by_order\":[";bool f=true;for(auto x:by){if(!f)cout<<",";f=false;cout<<"{\"n\":"<<x[0]<<",\"tournaments\":"<<x[1]<<",\"sole_exit_roots\":"<<x[2]<<",\"positive_endpoint_roots\":"<<x[3]<<",\"strong_roots\":"<<x[4]<<",\"nontrivial_compressions\":"<<x[5]<<",\"defect_three_roots\":"<<x[6]<<"}";}cout<<"],\"totals\":{\"sole_exit_roots\":"<<roots<<",\"positive_endpoint_roots\":"<<positive<<",\"strong_roots\":"<<strongs<<",\"nontrivial_compressions\":"<<compressed<<",\"defect_three_roots\":"<<dthree<<"},\"observed_defect_three_endpoints\":[";f=true;for(int x:spect){if(!f)cout<<",";f=false;cout<<x;}cout<<"]}\n";
}
