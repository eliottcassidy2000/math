// Exact finite menu: at most one reversed old arc, arbitrary re-marking,
// then the inherited two-birth reset. All subset counts fit int64 for n<=18:
// E2 is bounded by n*n!, smaller than 19!, and H<=18!.
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;using I=int64_t;using T=vector<unsigned>;
struct Stats{I h;vector<I> E;};
Stats stats(const T&t){int n=t.size(),f=(1<<n)-1;vector<I>e((f+1)*n),h(f+1);h[0]=1;
 for(int i=0;i<n;i++)e[(1<<i)*n+i]=1;
 for(int m=1;m<=f;m++)for(unsigned z=m;z;z&=z-1){int v=countr_zero(z);I a=e[m*n+v];h[m]+=a;
  for(unsigned q=t[v]&~m;q;q&=q-1)e[(m|(q&-q))*n+countr_zero(q)]+=a;}
 vector<I>E(n);for(int m=1;m<=f;m++)for(unsigned z=m;z;z&=z-1){int v=countr_zero(z);E[v]+=e[m*n+v]*h[f^m];}
 return {h[f],E};}
bool strong(const T&t){int n=t.size();unsigned f=(1<<n)-1;for(int rev=0;rev<2;rev++){unsigned r=1,old=0;while(r!=old){old=r;for(int i=0;i<n;i++)if(r>>i&1)r|=rev?(f^t[i]^(1<<i)):t[i];}if(r!=f)return false;}return true;}
void print(const T&t){cout<<'[';for(int i=0;i<(int)t.size();i++){if(i)cout<<',';cout<<t[i];}cout<<']';}
int main(){int n;I target;while(cin>>n>>target){if(n<1||n>18)return 2;T t(n);for(auto&x:t)cin>>x;
 auto base=stats(t);int candidates=0,hits=0;I nearest=INT64_MAX;cout<<"{\"n\":"<<n<<",\"source_H\":"<<base.h<<",\"target\":"<<target<<",\"source_E2\":[";for(int r=0;r<n;r++){if(r)cout<<',';cout<<base.E[r];}cout<<"],\"hits\":[";
 for(int a=-1;a<n;a++)for(int b=a+1;b<n;b++){if(a==-1&&b)continue;T u=t;if(a>=0){u[a]^=1<<b;u[b]^=1<<a;}auto s=a<0?base:stats(u);candidates++;
  for(int r=0;r<n;r++){I z=s.h+2*s.E[r];nearest=min(nearest,z>target?z-target:target-z);if(z==target){if(hits++)cout<<',';cout<<"{\"arc\":["<<a<<','<<b<<"],\"root\":"<<r<<",\"H\":"<<s.h<<",\"E2\":"<<s.E[r]<<",\"strong\":"<<(strong(u)?"true":"false")<<",\"rows\":";print(u);cout<<'}';}}}
 cout<<"],\"candidates\":"<<candidates<<",\"nearest\":"<<nearest<<"}\n"<<flush;
 }return 0;}
