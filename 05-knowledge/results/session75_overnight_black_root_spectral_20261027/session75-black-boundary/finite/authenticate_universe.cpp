// Independent ordinary-class authentication via score-block permutations.
// Input each line: n, upper-triangle MSB code, claimed automorphism order.
// Canonical keys independently reconstructed; orbit mass proves completeness.
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
struct Tot {uint64_t records=0,mass=0,perms=0;set<uint64_t> keys;};
int main(){try{
 map<int,Tot> totals;int n;uint64_t code,claimed;
 while(cin>>n>>code>>claimed){
  if(n<4||n>9)throw runtime_error("order outside gate");
  int e[9][9]={},d[9]={},bit=n*(n-1)/2-1;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,bit--){if((code>>bit)&1)e[i][j]=1,d[i]++;else e[j][i]=1,d[j]++;}
  vector<vector<int>> groups(n);for(int i=0;i<n;i++)groups[d[i]].push_back(i);
  vector<int> p;uint64_t best=UINT64_MAX,aut=0,perms=0;
  function<void(int)> go=[&](int g){
   if(g==n){uint64_t z=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)z=(z<<1)|e[p[i]][p[j]];perms++;if(z<best)best=z,aut=1;else if(z==best)aut++;return;}
   auto a=groups[g];do{for(int v:a)p.push_back(v);go(g+1);p.resize(p.size()-a.size());}while(next_permutation(a.begin(),a.end()));
  };go(0);
  if(aut!=claimed)throw runtime_error("automorphism mismatch n="+to_string(n)+" code="+to_string(code));
  auto &t=totals[n];if(!t.keys.insert(best).second)throw runtime_error("isomorphic duplicate");
  uint64_t fact=1;for(int i=2;i<=n;i++)fact*=i;if(fact%aut)throw runtime_error("nonintegral orbit");t.records++;t.mass+=fact/aut;t.perms+=perms;
 }
 if(totals.size()!=6)throw runtime_error("missing order");
 cout<<"{\n  \"status\": \"PASS\",\n  \"method\": \"independent score-block canonical keys, automorphism counts, exact labelled orbit mass\",\n  \"orders\": {\n";bool first=true;
 for(auto &[n,t]:totals){uint64_t target=uint64_t(1)<<(n*(n-1)/2);if(t.mass!=target)throw runtime_error("incomplete orbit mass");if(!first)cout<<",\n";first=false;cout<<"    \""<<n<<"\": {\"ordinary_classes\": "<<t.records<<", \"labelled_orbit_mass\": "<<t.mass<<", \"score_block_permutations\": "<<t.perms<<"}";}
 cout<<"\n  }\n}\n";
 }catch(const exception&e){cerr<<e.what()<<"\n";return 1;}return 0;}
