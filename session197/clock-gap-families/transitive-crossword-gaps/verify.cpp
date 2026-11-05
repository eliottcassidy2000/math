#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using R=vector<unsigned>;using Z=unsigned long long;
void need(bool b,const char*s){if(!b)throw runtime_error(s);}int bit(unsigned x){return __builtin_ctz(x);}Z power(Z b,int e){Z x=1;while(e--)x*=b;return x;}
R build(int i,int k,int j,int m){need(i>=0&&k>=0&&j>=0&&m>=0,"parameters");vector<int>w{1};for(int t=0;t<i;t++)w.push_back(0);w.push_back(1);for(int t=0;t<j;t++)w.push_back(0);w.push_back(2);for(int t=0;t<k;t++)w.push_back(0);w.push_back(2);int q=w.size(),n=3+q+m;R r(n);r[0]=(1<<1)|(1<<2);r[2]=1<<1;for(int v=0;v<q;v++){int x=v+3;for(int a=0;a<3;a++){bool out=a<2&&(w[v]>>a&1);if(out)r[a]|=1<<x;else r[x]|=1<<a;}for(int u=v+1;u<q;u++)r[x]|=1<<(u+3);for(int c=3+q;c<n;c++)r[x]|=1<<c;}for(int c=3+q;c<n;c++){r[c]|=7;for(int d=c+1;d<n;d++)r[c]|=1<<d;}return r;}
Z count(R const&r){int n=r.size(),N=1<<n;vector<Z>d(size_t(N)*n);for(int u=0;u<n;u++)d[size_t(1<<u)*n+u]=1;for(int mask=1;mask<N;mask++)for(int u=0;u<n;u++)if(d[size_t(mask)*n+u]){unsigned z=r[u]&~mask;while(z){unsigned b=z&-z;z-=b;d[size_t(mask|b)*n+bit(b)]+=d[size_t(mask)*n+u];}}Z h=0;for(int u=0;u<n;u++)h+=d[size_t(N-1)*n+u];return h;}
Z literal(R const&r){int n=r.size();unsigned full=(1<<n)-1;Z h=0;function<void(int,unsigned)>go=[&](int u,unsigned used){if(used==full){h++;return;}unsigned z=r[u]&~used;while(z){unsigned b=z&-z;z-=b;go(bit(b),used|b);}};for(int u=0;u<n;u++)go(u,1<<u);return h;}
int crossrank(R const&r,int root){int n=r.size();vector<int>mate(n,-1);function<bool(int,unsigned&)>aug=[&](int u,unsigned&seen){for(int v=0;v<n;v++)if(v!=root&&(r[v]>>root&1)&&(r[u]>>v&1)&&!(seen>>v&1)){seen|=1<<v;if(mate[v]<0||aug(mate[v],seen)){mate[v]=u;return true;}}return false;};int ans=0;for(int u=0;u<n;u++)if(r[root]>>u&1){unsigned seen=0;ans+=aug(u,seen);}return ans;}
bool prime(Z n){if(n<2)return false;for(Z p=2;p<=n/p;p++)if(n%p==0)return false;return true;}
void rows(R const&r){cout<<'[';for(int v=0;v<(int)r.size();v++){if(v)cout<<',';cout<<r[v];}cout<<']';}

long long pwr(int b,int e){return static_cast<long long>(power(b,e));}
int main(){int grid=0;cout<<"{\"status\":\"PASS\",\"grid\":[";bool first=true;for(int i=0;i<=2;i++)for(int k=0;k<=3;k++)for(int j=0;j<=2;j++)for(int m=0;m<=3;m++){auto r=build(i,k,j,m);long long U=i?pwr(2,i)-1:1,V=k?pwr(2,k)-1:1,F=pwr(2,i+1)+1;long long A=4*F*(pwr(3,k+1)+pwr(2,k+1))*pwr(2,j);long long D=-2*(F*((long long)pwr(3,k)-(long long)pwr(2,k+1)+1)*pwr(2,j)+((long long)pwr(2,k+1)+1)*(U-1));long long C=U*V-2*(j==0),h=A*pwr(3,m)+D*pwr(2,m)+C;need(h>0&&count(r)==(Z)h&&literal(r)==(Z)h,"general gap tensor count");need(D<=0&&((D==0)==(i<=1&&k<=1)),"finite cancellation control");if(i==1)need(crossrank(r,4)==1,"i module rank one");if(j==1)need(crossrank(r,5+i)==1,"j module rank one");if(k==1)need(crossrank(r,6+i+j)==2,"k module rank two");if(m==1)need(crossrank(r,r.size()-1)==2,"m module rank two");if(!first)cout<<',';first=false;cout<<'['<<i<<','<<k<<','<<j<<','<<m<<','<<r.size()<<','<<A<<','<<D<<','<<C<<','<<h<<']';grid++;}need(grid==144,"complete tensor grid");cout<<"],\"graphs\":"<<grid<<"}\n";}
