// Complete reversal-quotiented Hamiltonian paths in the simple square-sum graph.
// No unguarded forcing at degree-two vertices. Every unfinished continuation
// has its current vertex as one endpoint, so it has at most one other leaf.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using namespace std;
int main(int argc,char**argv){try{
 if(argc!=3)throw runtime_error("usage paths maximum_n output.json");int maximum=stoi(argv[1]);if(maximum>40||maximum<1)throw runtime_error("order range1..40");ofstream out(argv[2]);out<<"{\"orders\":[";
 for(int n=1;n<=maximum;n++){
  vector<uint64_t>a(n);vector<pair<int,int>>edges;for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)for(int k=2;k*k<=2*n;k++)if(u+v+2==k*k){a[u]|=1ULL<<v;a[v]|=1ULL<<u;edges.push_back({u,v});}
  vector<int>leaves;bool isolated=false;for(int v=0;v<n;v++){if(__builtin_popcountll(a[v])==1)leaves.push_back(v);if(!a[v])isolated=true;}
  vector<vector<int>>paths;vector<int>p;uint64_t nodes=0,pruned=0;uint64_t full=(1ULL<<n)-1;
  function<void(uint64_t)>dfs=[&](uint64_t used){nodes++;if(used==full){if(!leaves.empty()||n==1||p.front()<p.back()){auto q=p;for(int&x:q)x++;paths.push_back(q);}return;}int last=p.back();uint64_t rem=(full^used)|(1ULL<<last);int otherleaves=0;uint64_t scan=full^used;
   while(scan){int v=__builtin_ctzll(scan);scan&=scan-1;int degree=__builtin_popcountll(a[v]&rem);if(!degree||degree==1&&++otherleaves>1){pruned++;return;}}
   uint64_t reached=1ULL<<last,front=reached;while(front){int v=__builtin_ctzll(front);front&=front-1;uint64_t add=a[v]&rem&~reached;reached|=add;front|=add;}if(reached!=rem){pruned++;return;}
   uint64_t avail=a[last]&~used;while(avail){int v=__builtin_ctzll(avail);avail&=avail-1;p.push_back(v);dfs(used|(1ULL<<v));p.pop_back();}};
  if(n==1){p={0};dfs(1);}else if(!isolated&&leaves.size()<=2){vector<int>starts;if(!leaves.empty())starts={leaves[0]};else{starts.resize(n);iota(starts.begin(),starts.end(),0);}for(int v:starts){p={v};dfs(1ULL<<v);}}
  for(auto&q:paths){auto r=q;reverse(r.begin(),r.end());q=min(q,r);}sort(paths.begin(),paths.end());if(adjacent_find(paths.begin(),paths.end())!=paths.end())throw runtime_error("duplicate reversal class");
  if(n>1)out<<',';out<<"{\"n\":"<<n<<",\"dfs_nodes\":"<<nodes<<",\"pruned\":"<<pruned<<",\"paths\":[";for(size_t i=0;i<paths.size();i++){if(i)out<<',';out<<'[';for(int j=0;j<n;j++){if(j)out<<',';out<<paths[i][j];}out<<']';}out<<"]}";cerr<<n<<" paths "<<paths.size()<<" nodes "<<nodes<<'\n';
 }out<<"]}\n";return 0;
 }catch(exception&e){cerr<<e.what()<<'\n';return 1;}}
