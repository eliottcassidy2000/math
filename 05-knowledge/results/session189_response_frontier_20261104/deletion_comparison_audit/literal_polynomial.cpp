#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <cstdint>
#include <stdexcept>
using namespace std;using I=int64_t;
template<class T>void arr(const vector<T>&a){cout<<'[';for(size_t i=0;i<a.size();i++){if(i)cout<<',';cout<<a[i];}cout<<']';}
int main(){int cases;cin>>cases;cout<<'[';for(int id=0;id<cases;id++){int n;cin>>n;if(n<1||n>9)throw runtime_error("literal bound");vector<int>r(n);for(int&i:r)cin>>i;vector<int>p(n);iota(p.begin(),p.end(),0);vector<int>bs{0};for(int k=1;k<=2;k++)for(int m=1;m<(1<<n);m++)if(__builtin_popcount((unsigned)m)==k)bs.push_back(m);int L=bs.size();vector<int>ix(1<<n,-1);for(int i=0;i<L;i++)ix[bs[i]]=i;vector<I>co(L*L),starts(n),ends(n);
auto term=[&](int ai,int ao,int ci,int coo){for(int u=0;u<=(ao!=0);u++)for(int v=0;v<=(coo!=0);v++){int a=ai|(u?ao:0),c=ci|(v?coo:0);co[ix[a]*L+ix[c]]+=((u+v)%2?-1:1);}};
auto distinct=[&](int a,int c){term(a?1<<p[a-1]:0,a<n?1<<p[a]:0,c?1<<p[c-1]:0,c<n?1<<p[c]:0);};
auto same=[&](int g){term(g?1<<p[g-1]:0,0,0,g<n?1<<p[g]:0);};
do{vector<int>bad;for(int g=1;g<n;g++)if(!(r[p[g-1]]>>p[g]&1))bad.push_back(g);if(bad.empty()){starts[p[0]]++;ends[p[n-1]]++;for(int a=0;a<=n;a++){same(a);for(int c=0;c<=n;c++)if(c!=a)distinct(a,c);}}else if(bad.size()==1){int a=bad[0];same(a);for(int c=0;c<=n;c++)if(c!=a){distinct(a,c);distinct(c,a);}}else if(bad.size()==2){distinct(bad[0],bad[1]);distinct(bad[1],bad[0]);}}while(next_permutation(p.begin(),p.end()));
vector<I>values(L*L);for(int i=0;i<L;i++)for(int j=0;j<L;j++)for(int s=bs[i];;s=(s-1)&bs[i]){for(int t=bs[j];;t=(t-1)&bs[j]){values[i*L+j]+=co[ix[s]*L+ix[t]];if(!t)break;}if(!s)break;}
if(id)cout<<',';cout<<"{\"rows\":";arr(r);cout<<",\"masks\":";arr(bs);cout<<",\"coefficients\":";arr(co);cout<<",\"values\":";arr(values);cout<<",\"starts\":";arr(starts);cout<<",\"ends\":";arr(ends);cout<<'}';}cout<<"]\n";}
