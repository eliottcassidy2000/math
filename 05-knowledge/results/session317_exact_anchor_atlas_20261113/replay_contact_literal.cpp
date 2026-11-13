#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
int main() {
  constexpr int n=19;
  constexpr std::uint64_t expected=65869119;
  std::array<int,n> row{};
  auto wins_apex=[](int j){return j==0 || (j>=9 && j!=9);};
  for (int i=0;i<n;++i) for (int j=0;j<n;++j) if (i!=j) {
    bool arc;
    if (i==n-1) arc=wins_apex(j);
    else if (j==n-1) arc=!wins_apex(i);
    else if (i/9!=j/9) arc=i<j;
    else if (i/3!=j/3) arc=(j/3-i/3+3)%3==1;
    else arc=(j-i+3)%3==1;
    if (arc) row[i]|=1<<j;
  }
  std::vector<std::uint64_t> dp((1<<n)*n);
  for (int i=0;i<n;++i) dp[(1<<i)*n+i]=1;
  for (int s=1;s<(1<<n);++s) for (int i=0;i<n;++i) if (s>>i&1) {
    int a=row[i]&~s;
    while (a) {
      int j=__builtin_ctz(a);
      a&=a-1;
      dp[(s|(1<<j))*n+j]+=dp[s*n+i];
    }
  }
  std::uint64_t h=0;
  for (int i=0;i<n;++i) h+=dp[((1<<n)-1)*n+i];
  std::cout << "ternary literal aperture r=2 " << h << '\n';
  return h==expected ? 0 : 1;
}
