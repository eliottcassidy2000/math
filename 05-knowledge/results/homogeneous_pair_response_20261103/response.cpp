// Standalone predecessor-subset Hamilton counter for actual port extensions.
#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
U count(const std::vector<unsigned>& rows) {
  const unsigned n=rows.size(), total=1u<<n;
  if(n>16) throw std::runtime_error("declared bound");
  std::vector<unsigned> pred(n);
  for(unsigned u=0;u<n;++u)
    for(unsigned v=0;v<n;++v) if(rows[u]>>v&1u) pred[v]|=1u<<u;
  std::vector<U> dp(static_cast<std::size_t>(total)*n,0);
  for(unsigned v=0;v<n;++v) dp[(1u<<v)*n+v]=1;
  for(unsigned mask=1;mask<total;++mask) {
    if(!(mask&(mask-1))) continue;
    unsigned vertices=mask;
    while(vertices) {
      unsigned v=__builtin_ctz(vertices); vertices&=vertices-1;
      unsigned before=mask^(1u<<v), incoming=pred[v]&before;
      U sum=0;
      while(incoming) {
        unsigned u=__builtin_ctz(incoming); incoming&=incoming-1;
        sum+=dp[static_cast<std::size_t>(before)*n+u];
      }
      dp[static_cast<std::size_t>(mask)*n+v]=sum;
    }
  }
  U result=0;
  for(unsigned v=0;v<n;++v) result+=dp[static_cast<std::size_t>(total-1)*n+v];
  return result;
}
int main() {
  std::ios::sync_with_stdio(false);
  unsigned cases; if(!(std::cin>>cases)) return 1;
  for(unsigned k=0;k<cases;++k) {
    unsigned n,m; std::cin>>n>>m;
    std::vector<unsigned> old(n);
    for(auto& row:old) std::cin>>row;
    for(unsigned j=0;j<m;++j) {
      unsigned a,c;std::cin>>a>>c;
      std::vector<unsigned> rows=old;
      rows.push_back(1u<<(n+1)); rows.push_back(0);
      for(unsigned u=0;u<n;++u) {
        if(a>>u&1u) rows[u]|=1u<<n;else rows[n]|=1u<<u;
        if(c>>u&1u) rows[u]|=1u<<(n+1);else rows[n+1]|=1u<<u;
      }
      std::cout<<count(rows)<<"\n";
    }
  }
}

