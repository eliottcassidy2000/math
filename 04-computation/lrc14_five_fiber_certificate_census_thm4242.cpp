// Exact sufficient-certificate census for the 5-fibre deletion lemma on
// S={10,15,20,30,264}, center 50.  This does not replace the full audit.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>

using i64 = std::int64_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int,18> C={
    8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::array<int,12> O={10,15,20,30,60,63,132,170,176,190,264,290};
constexpr std::array<int,7> CD={40,80,85,95,120,145,240};
constexpr std::array<int,11> CN={8,16,42,84,88,126,143,168,193,252,286};
constexpr std::array<std::array<int,5>,11> L={{
    {{10,15,264,0,0}},{{10,20,264,0,0}},{{10,30,264,0,0}},
    {{15,20,264,0,0}},{{15,30,264,0,0}},{{20,30,264,0,0}},
    {{10,15,20,264,0}},{{10,15,30,264,0}},{{10,20,30,264,0}},
    {{15,20,30,264,0}},{{10,15,20,30,264}}
}};
constexpr std::array<int,11> PETALS={3,3,3,3,3,3,4,4,4,4,5};
constexpr std::array<int,11> CORE={6,6,6,6,6,6,5,5,5,5,4};
constexpr std::array<int,11> CUTOFF={434,444,442,465,428,471,399,376,410,389,352};

i64 lcm_exact(i64 a,i64 b) {
    const i128 v=i128(a/std::gcd(a,b))*b;
    if (v>std::numeric_limits<i64>::max()) throw std::runtime_error("lcm");
    return i64(v);
}

bool safe_mid(int s,i64 a,i64 b,i64 d) {
    i128 r=i128(s)*(a+b)%(i128(2)*d);
    if (r<0) r+=i128(2)*d;
    return i128(7)*r>=d && i128(7)*r<=i128(13)*d;
}

std::pair<i64,i64> mass(std::vector<int> speeds) {
    std::sort(speeds.begin(),speeds.end());
    speeds.erase(std::unique(speeds.begin(),speeds.end()),speeds.end());
    i64 d=1;
    for (int s:speeds) d=lcm_exact(d,14LL*s);
    std::vector<i64> walls{0,d};
    for (int s:speeds) {
        const i64 u=d/(14LL*s);
        for (int j=0;j<s;++j) {
            walls.push_back((14LL*j+1)*u);
            walls.push_back((14LL*j+13)*u);
        }
    }
    std::sort(walls.begin(),walls.end());
    walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
    i64 n=0;
    for (std::size_t i=0;i+1<walls.size();++i) {
        bool good=true;
        for (int s:speeds) if (!safe_mid(s,walls[i],walls[i+1],d)) {
            good=false; break;
        }
        if (good) n+=walls[i+1]-walls[i];
    }
    return {n,d};
}

bool excluded(int r) {
    if (r==50) return true;
    return std::find(C.begin(),C.end(),r)!=C.end() ||
           std::find(O.begin(),O.end(),r)!=O.end();
}

i64 choose(int n,int k) {
    if (k<0 || k>n) return 0;
    i64 v=1;
    for (int j=1;j<=k;++j) v=v*(n-j+1)/j;
    return v;
}

} // namespace

int main() {
    std::map<std::vector<int>,std::pair<i64,i64>> cache;
    i64 total=0,certified=0,div5_total=0,div5_certified=0;
    std::array<i64,11> lane_total{},lane_certified{};
    i128 min_positive_gap=-1;
    std::vector<int> min_gap_speeds;
    int min_gap_e=0;

    for (int lane=0;lane<11;++lane) {
        std::vector<int> fixed={10}; // center 50 / 5
        for (int j=0;j<PETALS[lane];++j) {
            if (L[lane][j]%5==0) fixed.push_back(L[lane][j]/5);
        }
        for (int r=1;r<CUTOFF[lane];++r) {
            if (excluded(r)) continue;
            const bool rdiv=r%5==0;
            for (int dmask=0;dmask<(1<<7);++dmask) {
                const int dc=__builtin_popcount(unsigned(dmask));
                const int nc=CORE[lane]-dc;
                const i64 multiplicity=choose(11,nc);
                if (multiplicity==0) continue;
                lane_total[lane]+=multiplicity;
                total+=multiplicity;
                if (rdiv) div5_total+=multiplicity;
                const int ecount=1+(rdiv?0:1)+nc; // 264, maybe r, non-5 core
                if (ecount>=5) continue;
                std::vector<int> speeds=fixed;
                if (rdiv) speeds.push_back(r/5);
                for (int j=0;j<7;++j) if (dmask&(1<<j)) speeds.push_back(CD[j]/5);
                std::sort(speeds.begin(),speeds.end());
                speeds.erase(std::unique(speeds.begin(),speeds.end()),speeds.end());
                auto it=cache.find(speeds);
                if (it==cache.end()) it=cache.emplace(speeds,mass(speeds)).first;
                const auto [n,den]=it->second;
                const i128 gap=i128(63)*(5-ecount)*n-i128(20)*den;
                if (gap<0) continue;
                lane_certified[lane]+=multiplicity;
                certified+=multiplicity;
                if (rdiv) div5_certified+=multiplicity;
                if (min_positive_gap<0 || gap*cache[min_gap_speeds].second <
                    min_positive_gap*den) {
                    min_positive_gap=gap;
                    min_gap_speeds=speeds;
                    min_gap_e=ecount;
                }
            }
        }
    }
    for (int lane=0;lane<11;++lane) {
        std::cout << "LANE labels=";
        for (int j=0;j<PETALS[lane];++j) {
            if (j) std::cout << ',';
            std::cout << L[lane][j];
        }
        std::cout << " total=" << lane_total[lane]
                  << " fiber_certified=" << lane_certified[lane]
                  << " uncertified=" << lane_total[lane]-lane_certified[lane]
                  << '\n';
    }
    std::cout << "SUMMARY total=" << total
              << " fiber_certified=" << certified
              << " uncertified=" << total-certified
              << " div5_total=" << div5_total
              << " div5_certified=" << div5_certified
              << " cached_divisible_systems=" << cache.size()
              << " weakest_nonnegative_scaled_gap=";
    if (min_positive_gap>=0) {
        const i64 den=cache[min_gap_speeds].second;
        std::cout << static_cast<long long>(min_positive_gap) << '/' << den
                  << " E=" << min_gap_e << " D=";
        for (std::size_t i=0;i<min_gap_speeds.size();++i) {
            if (i) std::cout << ',';
            std::cout << min_gap_speeds[i];
        }
    }
    std::cout << " status=FINITE-EXACT-SUFFICIENT-CERTIFICATE-NOT-EQUIVALENCE\n";
}
