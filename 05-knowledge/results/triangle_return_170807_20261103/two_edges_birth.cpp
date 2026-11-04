#define BIRTH_LIBRARY
#include "two_births.cpp"
int main(){int n;I target;cin>>n>>target;if(n>17)return 2;T t(n);for(auto&x:t)cin>>x;vector<pair<int,int>>arcs;for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)arcs.push_back({a,b});int checked=0;
 for(int i=0;i<(int)arcs.size();i++)for(int j=0;j<i;j++){auto[a,b]=arcs[i];auto[c,d]=arcs[j];T u=t;u[a]^=1<<b;u[b]^=1<<a;u[c]^=1<<d;u[d]^=1<<c;auto z=births(u);checked++;
 for(int cut=0;cut<(1<<n);cut++)if(z[cut]==target&&strong(birth(u,cut))){cout<<"{\"found\":true,\"target\":"<<target<<",\"arcs\":[["<<a<<','<<b<<"],["<<c<<','<<d<<"]],\"cut\":"<<cut<<",\"intermediate_H\":"<<z[0]<<",\"checked\":"<<checked<<",\"rows\":";print(birth(u,cut));cout<<"}\n";return 0;}}
 cout<<"{\"found\":false,\"target\":"<<target<<",\"checked\":"<<checked<<"}\n";return 0;}
