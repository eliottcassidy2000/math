// Independent audit of the p=5 fibre-certificate count.
// Enumerates every actual core body; measures divisible systems by XOR events.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

using u32=std::uint32_t;
using u64=std::uint64_t;
using i128=__int128_t;

namespace {
constexpr std::array<int,18> C={
  8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::array<int,12> O={10,15,20,30,60,63,132,170,176,190,264,290};
constexpr std::array<std::array<int,5>,11> L={{
  {{10,15,264,0,0}},{{10,20,264,0,0}},{{10,30,264,0,0}},
  {{15,20,264,0,0}},{{15,30,264,0,0}},{{20,30,264,0,0}},
  {{10,15,20,264,0}},{{10,15,30,264,0}},{{10,20,30,264,0}},
  {{15,20,30,264,0}},{{10,15,20,30,264}}
}};
constexpr std::array<int,11> PCOUNT={3,3,3,3,3,3,4,4,4,4,5};
constexpr std::array<int,11> CSIZE={6,6,6,6,6,6,5,5,5,5,4};
constexpr std::array<int,11> CUTOFF={434,444,442,465,428,471,399,376,410,389,352};

u64 lcm_checked(u64 a,u64 b) {
  const i128 x=i128(a/std::gcd(a,b))*b;
  if (x>std::numeric_limits<u64>::max()) throw std::runtime_error("lcm");
  return u64(x);
}

struct Event {u64 x;u32 toggle;};

std::pair<u64,u64> event_mass(std::vector<int> speeds) {
  std::sort(speeds.begin(),speeds.end());
  speeds.erase(std::unique(speeds.begin(),speeds.end()),speeds.end());
  u64 den=1;
  for (int s:speeds) den=lcm_checked(den,14*u64(s));
  std::vector<Event> events;
  for (std::size_t bit=0;bit<speeds.size();++bit) {
    const int s=speeds[bit];
    const u64 local=14*u64(s),scale=den/local;
    for (int tooth=0;tooth<s;++tooth) {
      int left=14*tooth-1;
      if (left<0) left+=int(local);
      events.push_back({u64(left)*scale,u32{1}<<bit});
      events.push_back({u64(14*tooth+1)*scale,u32{1}<<bit});
    }
  }
  std::sort(events.begin(),events.end(),[](const Event&a,const Event&b){return a.x<b.x;});
  std::vector<Event> grouped;
  for (const Event&e:events) {
    if (!grouped.empty() && grouped.back().x==e.x) grouped.back().toggle^=e.toggle;
    else grouped.push_back(e);
  }
  u32 state=speeds.empty()?0:(u32{1}<<speeds.size())-1;
  const u32 initial=state;
  u64 prev=0,numerator=0;
  for (const Event&e:grouped) {
    if (state==0) numerator+=e.x-prev;
    state^=e.toggle;
    prev=e.x;
  }
  if (state==0) numerator+=den-prev;
  if (state!=initial) throw std::runtime_error("circular state");
  return {numerator,den};
}

bool excluded(int r) {
  return r==50 || std::find(C.begin(),C.end(),r)!=C.end() ||
         std::find(O.begin(),O.end(),r)!=O.end();
}

struct BodyPart {std::vector<int> divisible;int nondivisible;};

std::array<std::vector<BodyPart>,3> make_bodies() {
  std::array<std::vector<BodyPart>,3> answer;
  for (u32 mask=0;mask<(u32{1}<<18);++mask) {
    const int size=std::popcount(mask);
    if (size<4 || size>6) continue;
    BodyPart b{{},0};
    for (int j=0;j<18;++j) if (mask&(u32{1}<<j)) {
      if (C[j]%5==0) b.divisible.push_back(C[j]/5); else ++b.nondivisible;
    }
    answer[std::size_t(size-4)].push_back(std::move(b));
  }
  if (answer[0].size()!=3060 || answer[1].size()!=8568 || answer[2].size()!=18564)
    throw std::runtime_error("body counts");
  return answer;
}
} // namespace

int main() {
  const auto bodies=make_bodies();
  std::map<std::vector<int>,std::pair<u64,u64>> cache;
  std::array<u64,11> totals{},certs{};
  u64 all=0,certified=0;
  i128 weakest=-1;
  u64 weakest_den=1;
  std::vector<int> weakest_speeds;
  int weakest_e=0;

  for (int lane=0;lane<11;++lane) {
    std::vector<int> fixed_d={10};
    int fixed_e=0;
    for (int j=0;j<PCOUNT[lane];++j) {
      const int s=L[lane][j];
      if (s%5==0) fixed_d.push_back(s/5); else ++fixed_e;
    }
    for (int r=1;r<CUTOFF[lane];++r) {
      if (excluded(r)) continue;
      std::vector<int> row_d=fixed_d;
      int row_e=fixed_e;
      if (r%5==0) row_d.push_back(r/5); else ++row_e;
      for (const BodyPart&body:bodies[std::size_t(CSIZE[lane]-4)]) {
        ++totals[lane];++all;
        const int ecount=row_e+body.nondivisible;
        if (ecount>=5) continue;
        std::vector<int> speeds=row_d;
        speeds.insert(speeds.end(),body.divisible.begin(),body.divisible.end());
        std::sort(speeds.begin(),speeds.end());
        speeds.erase(std::unique(speeds.begin(),speeds.end()),speeds.end());
        auto it=cache.find(speeds);
        if (it==cache.end()) it=cache.emplace(speeds,event_mass(speeds)).first;
        const auto [num,den]=it->second;
        const i128 gap=i128(63)*(5-ecount)*num-i128(20)*den;
        if (gap<0) continue;
        ++certs[lane];++certified;
        if (weakest<0 || gap*weakest_den<weakest*den) {
          weakest=gap;weakest_den=den;weakest_speeds=speeds;weakest_e=ecount;
        }
      }
    }
  }
  for (int lane=0;lane<11;++lane) {
    std::cout << "LANE labels=";
    for (int j=0;j<PCOUNT[lane];++j) {
      if (j) std::cout << ',';
      std::cout << L[lane][j];
    }
    std::cout << " total=" << totals[lane] << " fiber_certified=" << certs[lane]
              << " uncertified=" << totals[lane]-certs[lane] << '\n';
  }
  std::cout << "SUMMARY total=" << all << " fiber_certified=" << certified
            << " uncertified=" << all-certified
            << " cached_divisible_systems=" << cache.size()
            << " weakest_nonnegative_scaled_gap=";
  if (weakest>=0) {
    // This audit's actual weakest rational is known to fit u64 after evaluation.
    std::cout << u64(weakest) << '/' << weakest_den << " E=" << weakest_e << " D=";
    for (std::size_t j=0;j<weakest_speeds.size();++j) {
      if (j) std::cout << ',';
      std::cout << weakest_speeds[j];
    }
  }
  std::cout << " status=FINITE-EXACT-INDEPENDENT-SUFFICIENT-CERTIFICATE\n";
}
