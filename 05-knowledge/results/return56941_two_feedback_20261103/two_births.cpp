#define main inherited_search_main
#include "../adaptive_four_returns_20261103/lift_search.cpp"
#undef main
int main(){int n;I H;cin>>n>>H;T t(n);for(auto&x:t)cin>>x;int f=(1<<n)-1;
 for(int r=0;r<n;r++)for(int side=0;side<2;side++){int first=side?f^(1<<r):1<<r;T u=birth(t,first);auto v=ears(u);
  for(int b:{2,4})for(int cut=1;cut<(int)v.size()-1;cut++)if(v[cut]==3*H+b)cout<<"{\"H\":"<<H<<",\"offset\":"<<b<<",\"first\":"<<first<<",\"second\":"<<cut<<",\"pivot\":"<<r<<",\"side\":"<<side<<"}\n";
  cerr<<"done "<<r<<' '<<side<<'\n';cout.flush();
 }}
