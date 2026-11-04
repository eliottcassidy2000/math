#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using I=uint64_t;
struct Data { int n,N; vector<I> d,h; };
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
void valid(const vector<int>&r){int n=r.size();for(int i=0;i<n;i++){need(!(r[i]>>n)&&!(r[i]>>i&1),"invalid row");for(int j=0;j<i;j++)need(((r[i]>>j&1)+(r[j]>>i&1))==1,"not tournament");}}
Data table(const vector<int>&r){int n=r.size(),N=1<<n;Data a{n,N,vector<I>(N*n),vector<I>(N)};a.h[0]=1;for(int i=0;i<n;i++)a.d[(1<<i)*n+i]=1;for(int M=1;M<N;M++)for(int i=0;i<n;i++){I z=a.d[M*n+i];if(!z)continue;a.h[M]+=z;int rem=r[i]&~M;while(rem){int b=rem&-rem;rem-=b;a.d[(M|b)*n+__builtin_ctz((unsigned)b)]+=z;}}return a;}
vector<I> response(const Data&a){vector<I> F(a.n);for(int M=1;M<a.N;M++)for(int v=0;v<a.n;v++)F[v]+=a.d[M*a.n+v]*a.h[(a.N-1)^M];return F;}
vector<int> converse(const vector<int>&r){auto a=r;for(int i=0;i<(int)r.size();i++)a[i]^=((1<<r.size())-1)^(1<<i);return a;}
bool strong(const vector<int>&r){int n=r.size(),all=(1<<n)-1;for(int i=0;i<n;i++){int S=1<<i;while(true){int Q=S;for(int j=0;j<n;j++)if(S>>j&1)Q|=r[j];if(Q==S)break;S=Q;}if(S!=all)return false;}return true;}
vector<int> deleted(const vector<int>&r,int v){int n=r.size();vector<int>a(n-1);for(int i=0,ii=0;i<n;i++)if(i!=v){for(int j=0,jj=0;j<n;j++)if(j!=v){a[ii]|=(r[i]>>j&1)<<jj;jj++;}ii++;}return a;}
int main(int argc,char**argv){need(argc>=2,"mode");string mode=argv[1];
if(mode=="small") {set<pair<I,I>> prof;for(int n=3;n<=6;n++){I total=0,st=0,unit=0;for(unsigned m=0;m<(1u<<(n*(n-1)/2));m++){vector<int>r(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(m>>k&1)r[i]|=1<<j;else r[j]|=1<<i;total++;if(!strong(r))continue;st++;auto a=table(r),b=table(converse(r));auto F=response(a);for(int v=0;v<n;v++)if(b.d[(a.N-1)*n+v]==1){unit++;prof.insert({a.h.back(),F[v]});}}cout<<"SMALL "<<n<<" "<<total<<" "<<st<<" "<<unit<<"\n";}I matches=0;for(auto[h,f]:prof)if(prof.count({h+2,f}))matches++;need(matches==0,"unexpected small twin");cout<<"PROFILES "<<prof.size()<<" MATCHES "<<matches<<"\n";return 0;}
need(mode=="bank"||mode=="extend","unknown mode");int n,K;cin>>n>>K;need(n>=3&&n<=9&&K>0,"bank parameters");I tested=0,strongcount=0,low=0,found=0,max7=0;map<I,I> targets;map<I,set<I>> rootF;set<vector<int>> lowrows;
auto check=[&](const vector<int>&r){valid(r);tested++;auto a=table(r);I h=a.h.back();if(r.size()==7)max7=max(max7,h);if(!strong(r))return;strongcount++;if(h<=201&&r.size()==9){low++;lowrows.insert(r);}if(h==191||h==201){targets[h]++;auto F=response(a);for(auto f:F){rootF[h].insert(f);if((h==191&&f==1273)||(h==201&&f==5303))found++;}}};
for(int q=0;q<K;q++){vector<int>r(n);for(int&i:r)cin>>i;valid(r);if(mode=="bank"){check(r);if(n==8)for(int v=0;v<8;v++)check(deleted(r,v));}else for(int M=0;M<(1<<n);M++){auto t=r;for(int v=0;v<n;v++)if(M>>v&1)t[v]|=1<<n;t.push_back(((1<<n)-1)^M);check(t);}}
need(found==0,"target profile was realized");cout<<"SUMMARY "<<mode<<" "<<n<<" "<<K<<" "<<tested<<" "<<strongcount<<" "<<low<<" "<<found<<" "<<max7<<"\n";for(auto[h,count]:targets){cout<<"TARGET "<<h<<" "<<count<<" "<<rootF[h].size();for(auto f:rootF[h])cout<<" "<<f;cout<<"\n";}for(auto&r:lowrows){cout<<"LOW";for(int x:r)cout<<" "<<x;cout<<"\n";}}
