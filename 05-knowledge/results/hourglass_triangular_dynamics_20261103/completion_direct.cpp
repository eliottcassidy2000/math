#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>
using namespace std;
long long census(vector<int> F) {
  int n=F.size(); vector<pair<int,int>> free;
  vector<unsigned> fixed(n,0);
  for(int a=0;a<n;a++) {
    int b=F[a]; if(a==b || F[b]==a) throw runtime_error("invalid partial orientation");
    fixed[a]|=1u<<b;
  }
  for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)
    if(F[a]!=b && F[b]!=a)free.push_back({a,b});
  long long good=0;
  for(unsigned bits=0;bits<(1u<<free.size());bits++) {
    vector<unsigned> out=fixed,in(n,0);
    for(unsigned k=0;k<free.size();k++) {
      auto [a,b]=free[k]; if((bits>>k)&1)out[a]|=1u<<b;else out[b]|=1u<<a;
    }
    for(int a=0;a<n;a++)for(int b=0;b<n;b++)if(out[a]>>b&1)in[b]|=1u<<a;
    auto reaches_all=[&](const vector<unsigned>& adj) {
      unsigned seen=1,old=0;
      while(old!=seen) {
        old=seen;
        for(int a=0;a<n;a++)if(old>>a&1)seen|=adj[a];
      }
      return seen==(1u<<n)-1;
    };
    if(reaches_all(out)&&reaches_all(in))good++;
  }
  return good;
}
int main() {
  cout << "{\"full_strong\":" << census({6,3,1,0,0,1,3,6})
       << ",\"quotient_strong\":" << census({3,3,0,2}) << "}\n";
}
