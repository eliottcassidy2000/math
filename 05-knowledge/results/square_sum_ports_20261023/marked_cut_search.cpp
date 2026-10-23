// G25 with proposed selected edge1--15 and mandatory leaf endpoint18.
// Search every cut that can possibly violate the endpoint/selected-edge budget.
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
using namespace std;
int main(int argc,char**argv){if(argc!=2)return 2;const int n=25;uint64_t full=(1ULL<<n)-1,edge=(1ULL<<0)|(1ULL<<14),leaf=1ULL<<17;vector<uint64_t>a(n),found(n);uint64_t unresolved=full^leaf,tested=0;
 for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)for(int r=2;r*r<=2*n;r++)if(u+v+2==r*r){a[u]|=1ULL<<v;a[v]|=1ULL<<u;}
 int largest=0;for(int k=0;k<=13&&unresolved;k++){
  function<void(int,int,uint64_t)>visit=[&](int start,int left,uint64_t S){if(!unresolved)return;if(left){for(int v=start;v<=n-left;v++)visit(v+1,left-1,S|(1ULL<<v));return;}tested++;uint64_t rem=full^S;int c=0;while(rem){c++;uint64_t front=rem&-rem;rem^=front;while(front){int v=__builtin_ctzll(front);front&=front-1;uint64_t add=a[v]&rem;rem^=add;front|=add;}}int budget=k+1-c;int mandatory=((S&leaf)!=0)+((S&edge)==edge);uint64_t candidates=unresolved;while(candidates){int t=__builtin_ctzll(candidates);candidates&=candidates-1;if(mandatory+((S>>t)&1)>budget){found[t]=S;unresolved&=~(1ULL<<t);}}};visit(0,k,0);largest=k;}
 ofstream out(argv[1]);out<<"{\"n\":25,\"required_edge\":[1,15],\"mandatory_leaf_endpoint\":18,\"searched_through_size\":"<<largest<<",\"cuts_tested\":"<<tested<<",\"all_potentially_violating_sizes\":true,\"endpoints\":[";bool first=true;for(int t=0;t<n;t++)if(t!=17){if(!first)out<<',';first=false;out<<"{\"other_endpoint\":"<<t+1<<",\"ruled_out\":"<<(unresolved>>t&1?"false":"true")<<",\"minimal_violating_cut\":[";bool sep=false;for(int v=0;v<n;v++)if(found[t]>>v&1){if(sep)out<<',';sep=true;out<<v+1;}out<<"]}";}out<<"]}\n";cerr<<"marked-cut tested "<<tested<<" unresolved "<<__builtin_popcountll(unresolved)<<'\n';}
