// Complete tournament classes through order eight, by literal vertex extension.
// Canonicalization uses all permutations within sorted out-degree cells.
#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
vector<int> canonical(const vector<int>&r){
 int n=r.size();map<int,vector<int>>bydegree;for(int i=0;i<n;i++)bydegree[__builtin_popcount((unsigned)r[i])].push_back(i);
 vector<vector<int>>cells;for(auto&kv:bydegree)cells.push_back(kv.second);
 vector<int>order,best;bool found=false;
 function<void(int)>visit=[&](int g){
  if(g==(int)cells.size()){
   vector<int>out(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(r[order[i]]>>order[j]&1)out[i]|=1<<j;
   if(!found||out<best){best=out;found=true;}return;
  }
  auto c=cells[g];sort(c.begin(),c.end());
  do{order.insert(order.end(),c.begin(),c.end());visit(g+1);order.resize(order.size()-c.size());}while(next_permutation(c.begin(),c.end()));
 };visit(0);return best;
}
int main(){
 vector<vector<int>>level(1), all; int expected[]={1,1,1,2,4,12,56,456,6880};
 for(int n=1;n<=8;n++){
  set<vector<int>>next;
  for(const auto&r:level)for(int mask=0;mask<(1<<(n-1));mask++){
   vector<int>child=r;for(int i=0;i<n-1;i++)if(!(mask>>i&1))child[i]|=1<<(n-1);child.push_back(mask);next.insert(canonical(child));
  }
  level.assign(next.begin(),next.end()); for(const auto& row:level) all.push_back(row); if(level.size()!=(size_t)expected[n])throw runtime_error("class count mismatch");cerr<<n<<" "<<level.size()<<"\n";
 }
 cout << all.size() << "\n"; for (const auto& row: all) { cout << row.size(); for (int x: row) cout << " " << x; cout << "\n"; }
}
