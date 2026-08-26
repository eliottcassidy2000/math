#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

constexpr int N = 11;
constexpr int FULL = (1 << N) - 1;

struct OrbitEncoding {
  int variables = 0;
  std::array<int, N * (N - 1) / 2> variable{};
  std::array<int, N * (N - 1) / 2> positive{};
};

int edge_index(int a, int b) {
  if (a > b) std::swap(a, b);
  int z = 0;
  for (int i = 0; i < a; ++i) z += N - 1 - i;
  return z + b - a - 1;
}

std::array<int, N> permutation(const std::string& type) {
  std::array<int, N> p{};
  std::iota(p.begin(), p.end(), 0);
  auto cycle = [&](std::initializer_list<int> c) {
    std::vector<int> v(c);
    for (size_t i = 0; i < v.size(); ++i) p[v[i]] = v[(i + 1) % v.size()];
  };
  if (type == "3x2") {
    cycle({0,1,2}); cycle({3,4,5});
  } else if (type == "3x3") {
    cycle({0,1,2}); cycle({3,4,5}); cycle({6,7,8});
  } else if (type == "5x2") {
    cycle({0,1,2,3,4}); cycle({5,6,7,8,9});
  } else if (type == "11") {
    cycle({0,1,2,3,4,5,6,7,8,9,10});
  } else {
    std::cerr << "unknown type: " << type << '\n';
    std::exit(2);
  }
  return p;
}

OrbitEncoding make_encoding(const std::array<int, N>& p) {
  OrbitEncoding e;
  e.variable.fill(-1);
  e.positive.fill(-1);
  for (int a = 0; a < N; ++a) for (int b = a + 1; b < N; ++b) {
    int ix = edge_index(a,b);
    if (e.variable[ix] >= 0) continue;
    const int variable = e.variables++;
    int x = a, y = b;
    do {
      int u = std::min(x,y), v = std::max(x,y);
      int j = edge_index(u,v);
      int positive = (x < y); // representative a->b maps to lower->higher iff x<y.
      if (e.variable[j] >= 0 &&
          (e.variable[j] != variable || e.positive[j] != positive)) {
        std::cerr << "edge orbit reversal conflict\n";
        std::exit(3);
      }
      e.variable[j] = variable;
      e.positive[j] = positive;
      x = p[x]; y = p[y];
    } while (x != a || y != b);
  }
  return e;
}

std::array<uint16_t, N> decode(uint64_t assignment, const OrbitEncoding& e) {
  std::array<uint16_t, N> out{};
  int ix = 0;
  for (int a = 0; a < N; ++a) for (int b = a + 1; b < N; ++b, ++ix) {
    bool rep = (assignment >> e.variable[ix]) & 1;
    bool lower_wins = (rep == bool(e.positive[ix]));
    if (lower_wins) out[a] |= uint16_t(1u << b);
    else out[b] |= uint16_t(1u << a);
  }
  return out;
}

bool is_prime(const std::array<uint16_t, N>& out) {
  for (int s = 0; s <= FULL; ++s) {
    const int k = __builtin_popcount(unsigned(s));
    if (k < 2 || k == N) continue;
    bool module = true;
    uint16_t outside = uint16_t(FULL ^ s);
    while (outside) {
      int x = __builtin_ctz(outside);
      outside &= uint16_t(outside - 1);
      int toward = __builtin_popcount(unsigned(out[x] & s));
      if (toward != 0 && toward != k) { module = false; break; }
    }
    if (module) return false;
  }
  return true;
}

std::string label(const std::array<uint16_t, N>& out) {
  std::string s;
  s.reserve(N * (N - 1) / 2);
  for (int a = 0; a < N; ++a) for (int b = a + 1; b < N; ++b)
    s.push_back((out[a] & (1u << b)) ? '1' : '0');
  return s;
}

std::string digraph6(const std::array<uint16_t, N>& out) {
  std::string s;
  s.push_back('&');
  s.push_back(char(63 + N));
  int value = 0, used = 0;
  for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
    value = (value << 1) | ((out[i] >> j) & 1u);
    if (++used == 6) {
      s.push_back(char(63 + value));
      value = used = 0;
    }
  }
  if (used) s.push_back(char(63 + (value << (6 - used))));
  return s;
}

std::string score_key(const std::array<uint16_t, N>& out) {
  std::array<int,N> d{};
  for (int i=0;i<N;++i) d[i]=__builtin_popcount(unsigned(out[i]));
  std::sort(d.begin(),d.end());
  std::ostringstream z;
  for (int i=0;i<N;++i) { if(i) z<<','; z<<d[i]; }
  return z.str();
}

} // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
  _setmode(_fileno(stdout), _O_BINARY);
#endif
  std::string type = "3x2";
  bool emit = false, emit_d6 = false;
  uint64_t residue = 0, modulus = 1;
  for (int i=1;i<argc;++i) {
    std::string a=argv[i];
    if(a=="--type" && i+1<argc) type=argv[++i];
    else if(a=="--emit-prime") emit=true;
    else if(a=="--emit-digraph6") emit_d6=true;
    else if(a=="--res" && i+1<argc) residue=std::stoull(argv[++i]);
    else if(a=="--mod" && i+1<argc) modulus=std::stoull(argv[++i]);
    else { std::cerr<<"usage: --type 3x2|3x3|5x2|11 [--emit-prime|--emit-digraph6] [--res r --mod m]\n"; return 2; }
  }
  auto p=permutation(type);
  auto e=make_encoding(p);
  if(e.variables>=63){std::cerr<<"too many variables\n";return 3;}
  const uint64_t total=uint64_t(1)<<e.variables;
  if(!modulus || residue>=modulus){std::cerr<<"bad residue/modulus\n";return 2;}
  uint64_t prime=0,regular=0;
  std::map<std::string,uint64_t> scores;
  for(uint64_t a=residue;a<total;a+=modulus){
    auto out=decode(a,e);
    if(!is_prime(out))continue;
    ++prime;
    auto key=score_key(out);
    ++scores[key];
    if(key=="5,5,5,5,5,5,5,5,5,5,5")++regular;
    if(emit)std::cout<<label(out)<<'\n';
    else if(emit_d6)std::cout<<digraph6(out)<<'\n';
  }
  if(!emit && !emit_d6){
    std::cout<<"type="<<type<<" variables="<<e.variables<<" total="<<total
             <<" residue="<<residue<<" modulus="<<modulus
             <<" prime_rows="<<prime<<" regular_rows="<<regular
             <<" score_keys="<<scores.size()<<'\n';
    for(const auto& [key,count]:scores)std::cout<<"score="<<key<<" rows="<<count<<'\n';
  }
}
