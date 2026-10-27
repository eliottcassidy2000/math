// Phase2 independent class authentication: equitable refinement, then canonical
// maximal arc word and automorphism multiplicity. No producer source imports.
#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
struct Counts {uint64_t records=0,mass=0,permutations=0;set<uint64_t> keys;};
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
int main(){try{
 map<int,Counts> totals;int n;uint64_t code,claimed;
 while(cin>>n>>code>>claimed){
  need(n>=4&&n<=9,"bad n");need(code<(uint64_t(1)<<(n*(n-1)/2)),"extra high bits");
  bool e[9][9]={};vector<int> color(n,0);int bit=n*(n-1)/2-1;
  for(int u=0;u<n;u++)for(int v=u+1;v<n;v++){bool q=(code>>bit--)&1;e[u][v]=q;e[v][u]=!q;color[q?u:v]++;}
  while(true){
   vector<vector<int>> signature(n);set<vector<int>> types;
   for(int u=0;u<n;u++){signature[u]=vector<int>(n+1,0);signature[u][0]=color[u];for(int v=0;v<n;v++)if(e[u][v])signature[u][color[v]+1]++;types.insert(signature[u]);}
   map<vector<int>,int> ids;int k=0;for(auto &s:types)ids[s]=k++;
   vector<int> next(n);for(int u=0;u<n;u++)next[u]=ids.at(signature[u]);
   if(next==color)break;
   color=next;
  }
  vector<vector<int>> cells(*max_element(color.begin(),color.end())+1);
  for(int v=0;v<n;v++)cells[color[v]].push_back(v);
  vector<int> labels;bool first=true;uint64_t best=0,aut=0,perms=0;
  function<void(size_t)> enumerate=[&](size_t c){
   if(c==cells.size()){
    uint64_t word=0;for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)word=(word<<1)|e[labels[u]][labels[v]];
    perms++;if(first||word>best){first=false;best=word;aut=1;}else if(word==best)aut++;return;
   }
   vector<int> cell=cells[c];do{size_t old=labels.size();labels.insert(labels.end(),cell.begin(),cell.end());enumerate(c+1);labels.resize(old);}while(next_permutation(cell.begin(),cell.end()));
  };enumerate(0);
  need(aut==claimed,"automorphism mismatch");auto &t=totals[n];need(t.keys.insert(best).second,"isomorphic duplicate");
  uint64_t fact=1;for(int k=2;k<=n;k++)fact*=k;need(fact%aut==0,"bad orbit size");t.records++;t.mass+=fact/aut;t.permutations+=perms;
 }
 need(cin.eof(),"malformed input");need(totals.size()==6,"missing order");
 cout<<"{\n  \"status\": \"PASS\",\n  \"method\": \"equitable outgoing-color refinement, maximal canonical word, exact automorphism multiplicity and ordinary labelled orbit mass\",\n  \"orders\": {\n";bool comma=false;
 for(auto &[n,t]:totals){uint64_t target=uint64_t(1)<<(n*(n-1)/2);need(t.mass==target,"incomplete orbit mass");if(comma)cout<<",\n";comma=true;cout<<"    \""<<n<<"\": {\"ordinary_classes\": "<<t.records<<", \"labelled_mass\": "<<t.mass<<", \"refined_labelings\": "<<t.permutations<<"}";}
 cout<<"\n  }\n}\n";return 0;
 }catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
