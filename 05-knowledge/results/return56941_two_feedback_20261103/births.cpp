#define main inherited_search_main
#include "../adaptive_four_returns_20261103/lift_search.cpp"
#undef main
int main(){int n;I H;while(cin>>n>>H){T t(n);for(auto&x:t)cin>>x;auto v=ears(t);for(int b:{2,4}){cout<<H<<' '<<b;for(int mask=0;mask<(int)v.size();mask++)if(v[mask]==3*H+b)cout<<' '<<mask;cout<<'\n';}}}
