#define BIRTH_LIBRARY
#include "two_births.cpp"
int main(){int n;I target;while(cin>>n>>target){T t(n);for(auto&x:t)cin>>x;int examined=0,hits=0;cout<<"{\"target\":"<<target<<",\"hits\":[";
 for(int a=0;a<n;a++)for(int b=a+1;b<n;b++){T u=t;u[a]^=1<<b;u[b]^=1<<a;auto counts=births(u);examined++;
 for(int c=0;c<(1<<n);c++)if(counts[c]==target){if(hits++)cout<<',';cout<<"{\"arc\":["<<a<<','<<b<<"],\"cut\":"<<c<<",\"intermediate_H\":"<<counts[0]<<",\"strong\":"<<(strong(birth(u,c))?"true":"false")<<",\"rows\":";print(birth(u,c));cout<<'}'<<flush;}}
 cout<<"],\"examined\":"<<examined<<"}\n"<<flush;}return 0;}
