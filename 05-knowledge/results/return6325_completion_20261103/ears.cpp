#define main inherited_search_main
#include "../adaptive_four_returns_20261103/lift_search.cpp"
#undef main
int main(){int n;cin>>n;T t(n);for(auto&x:t)cin>>x;auto v=ears(t);cout<<'[';for(int i=0;i<(int)v.size();i++){if(i)cout<<',';cout<<v[i];}cout<<"]\n";}
