#define BIRTH_LIBRARY
#include "two_births.cpp"
int main(){int n;I target;cin>>n>>target;if(n>17)return 2;T t(n);for(auto&x:t)cin>>x;int rewrites=0,checked=0;
 for(int a=0;a<n;a++)for(int b=a+1;b<n;b++){T u=t;u[a]^=1<<b;u[b]^=1<<a;auto counts=births(u);rewrites++;if(counts[0]>target)continue;
 for(int c=0;c<(1<<n);c++)if(counts[c]<=target){auto v=birth(u,c);auto z=births(v);checked++;
 for(int d=0;d<(1<<(n+1));d++)if(z[d]==target&&strong(birth(v,d))){cout<<"{\"found\":true,\"target\":"<<target<<",\"arc\":["<<a<<','<<b<<"],\"cuts\":["<<c<<','<<d<<"],\"counts\":["<<counts[0]<<','<<counts[c]<<','<<target<<"],\"rewrites_checked\":"<<rewrites<<",\"births_checked\":"<<checked<<",\"rows\":";print(birth(v,d));cout<<"}\n";return 0;}}}
 cout<<"{\"found\":false,\"target\":"<<target<<",\"rewrites_checked\":"<<rewrites<<",\"births_checked\":"<<checked<<"}\n";return 0;}
