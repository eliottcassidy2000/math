// Direct exhaustive simple-cycle parity audit, independent of Walsh transform.
// g++ -O3 -std=c++17 -mpopcnt this_file.cpp -o /tmp/math-d5-sep05-audit
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

void require(bool b, const char* why) { if (!b) throw std::runtime_error(why); }

int main() {
  for (int n=6; n<=8; ++n) {
    int ix[8][8] = {}, dim=0;
    for (int a=1;a<n;++a) for (int b=a+1;b<n;++b) ix[a][b]=ix[b][a]=dim++;
    std::vector<uint32_t> masks;
    std::array<std::vector<uint32_t>,5> layer;
    std::vector<int> perm(n);
    for (int a=0;a<n;++a) perm[a]=a;
    // Generate all ordered distinct k-tuples by full permutations, then dedup.
    std::array<std::set<uint32_t>,5> sets;
    do {
      for (int k=3;k<=std::min(n,7);++k) {
        if (*std::min_element(perm.begin(),perm.begin()+k)!=perm[0] || perm[1]>perm[k-1]) continue;
        uint32_t m=0;
        for (int j=0;j<k;++j) {
          int a=perm[j], b=perm[(j+1)%k];
          if (a && b) m ^= uint32_t(1)<<ix[a][b];
        }
        sets[k-3].insert(m);
      }
    } while (std::next_permutation(perm.begin(),perm.end()));
    int expected_cycles[3][5]={{20,45,72,60,0},{35,105,252,420,360},{56,210,672,1680,2880}};
    for (int j=0;j<5;++j) {
      require(int(sets[j].size())==expected_cycles[n-6][j],"cycle universe");
      layer[j].assign(sets[j].begin(),sets[j].end());
      if(j<4) masks.insert(masks.end(),layer[j].begin(),layer[j].end());
    }
    std::set<uint32_t> expected;
    for (int j=0;j<dim;++j) expected.insert(uint32_t(1)<<j);
    for (int a=1;a<n;++a) {
      uint32_t m=0;
      for (int b=1;b<n;++b) if(a!=b) m |= uint32_t(1)<<ix[a][b];
      expected.insert(m);
    }
    int minimum=std::numeric_limits<int>::max();
    int minimum6=std::numeric_limits<int>::max();
    std::set<uint32_t> minimizers;
    std::set<uint32_t> minimizers6;
    std::map<std::array<int,4>,int> profiles;
    const uint32_t end=uint32_t(1)<<dim;
    for (uint32_t h=0;h<end;++h) {
      int total=0;
      int count7=0;
      if(n<=7) {
        std::array<int,5> c={};
        for(int j=0;j<5;++j) for(uint32_t m:layer[j]) c[j]+=__builtin_parity(h&m);
        if(n==6) {
          require(3*c[3]>=2*c[1],"six-vertex comparison");
          profiles[{c[0],c[1],c[2],c[3]}]++;
        } else require(2*c[4]>=c[1]+c[2],"seven-vertex comparison");
        for(int j=0;j<4;++j) total+=c[j];
        count7=c[4];
      } else {
        for(uint32_t m:masks) total+=__builtin_parity(h&m);
        for(uint32_t m:layer[4]) count7+=__builtin_parity(h&m);
      }
      if(!h) { require(total==0,"balanced hostile"); continue; }
      if(total<minimum) { minimum=total; minimizers.clear(); }
      if(total==minimum) minimizers.insert(h);
      if(n>=7) {
        int total6=total+count7;
        if(total6<minimum6) { minimum6=total6; minimizers6.clear(); }
        if(total6==minimum6) minimizers6.insert(h);
      }
    }
    const int values[3]={64,205,516};
    require(minimum==values[n-6],"D5 minimum");
    require(minimizers==expected,"D5 equality");
    if(n>=7) {
      require(minimum6==(n==7 ? 325 : 1236),"D6 minimum");
      require(minimizers6==expected,"D6 equality");
    }
    std::cout<<"n="<<n<<" classes="<<end<<" minimum="<<minimum<<" minimizers="<<minimizers.size()<<" gap="<<2*minimum<<" PASS\n";
    if(n==6) for(const auto& row:profiles) {
      std::cout<<"profile"; for(int x:row.first) std::cout<<" "<<x;
      std::cout<<" multiplicity="<<row.second<<"\n";
    }
    if(n>=7) std::cout<<"D6 n="<<n<<" minimum="<<minimum6<<" minimizers="<<minimizers6.size()<<" gap="<<2*minimum6<<" PASS\n";
  }
  std::cout<<"PASS: all 2130944 switching classes at n=6,7,8; direct parity only\n";
}
