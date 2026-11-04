#define BIRTH_LIBRARY
#include "two_births.cpp"
// Discovery: stop at first strong hit. Every tested predecessor is named.
int main(){int n;I target;cin>>n>>target;if(n>16)return 2;T t(n);for(auto&x:t)cin>>x;auto a=births(t);int first=0,second=0,weak=0;vector<vector<I>>prefixes;
 for(int m=0;m<(1<<n);m++)if(a[m]<=target){auto u=birth(t,m);auto b=births(u);first++;
 for(int c=0;c<(1<<(n+1));c++)if(b[c]<=target){auto v=birth(u,c);auto d=births(v);second++;prefixes.push_back({m,c,a[m],b[c]});
 for(int e=0;e<(1<<(n+2));e++)if(d[e]==target){if(!strong(birth(v,e))){weak++;continue;}cout<<"{\"found\":true,\"target\":"<<target<<",\"cuts\":["<<m<<','<<c<<','<<e<<"],\"counts\":["<<a[m]<<','<<b[c]<<','<<target<<"],\"first_checked\":"<<first<<",\"second_checked\":"<<second<<",\"rows\":";print(birth(v,e));cout<<"}\n";return 0;}}
 cerr<<"first "<<first<<" second "<<second<<'\n';}
 cout<<"{\"found\":false,\"target\":"<<target<<",\"first_checked\":"<<first<<",\"second_checked\":"<<second<<",\"weak_hits\":"<<weak<<",\"prefixes\":[";for(int i=0;i<(int)prefixes.size();i++){if(i)cout<<',';cout<<'[';for(int j=0;j<4;j++){if(j)cout<<',';cout<<prefixes[i][j];}cout<<']';}cout<<"]}\n";return 0;}
