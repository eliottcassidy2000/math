#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <chrono>
using namespace std;
int n,start;uint64_t full,nodes;vector<uint64_t>adj;vector<int>path;bool first;
void dfs(int v,uint64_t used){
 ++nodes;
 if(used==full){if(start< v){if(!first)cout<<",\n";first=false;cout<<"[";for(int i=0;i<n;i++){if(i)cout<<",";cout<<path[i]+1;}cout<<"]";}return;}
 uint64_t opts=adj[v]&~used;
 while(opts){int w=__builtin_ctzll(opts);opts&=opts-1;path.push_back(w);dfs(w,used|(1ULL<<w));path.pop_back();}
}
int main(int argc,char**argv){n=argc>1?atoi(argv[1]):33;if(n<2||n>60)return 2;full=(1ULL<<n)-1;adj.assign(n,0);for(int a=0;a<n;a++)for(int b=0;b<n;b++){int s=a+b+2,r=sqrt(s);if(a!=b&&r*r==s)adj[a]|=1ULL<<b;}first=true;nodes=0;cout<<"[\n";for(start=0;start<n;start++){path={start};dfs(start,1ULL<<start);}cout<<"\n]\n";cerr<<"n="<<n<<" nodes="<<nodes<<"\n";}
