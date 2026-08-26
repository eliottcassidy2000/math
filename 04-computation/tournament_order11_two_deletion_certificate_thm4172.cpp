#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using i64 = long long;
using i128 = __int128_t;

struct Tournament {
  int n = 0;
  std::array<std::uint16_t, 12> out{};
  std::string label;
};

static i64 choose2(i64 x) { return x * (x - 1) / 2; }

static std::string i128str(i128 x) {
  if (x == 0) return "0";
  bool neg = x < 0;
  if (neg) x = -x;
  std::string s;
  while (x) {
    s.push_back(char('0' + x % 10));
    x /= 10;
  }
  if (neg) s.push_back('-');
  std::reverse(s.begin(), s.end());
  return s;
}

static Tournament parse_bits(const std::string& bits) {
  Tournament t;
  int e = int(bits.size());
  while (t.n * (t.n - 1) / 2 < e) ++t.n;
  if (t.n * (t.n - 1) / 2 != e || t.n > 12) {
    throw std::runtime_error("nontriangular/oversize label");
  }
  t.label = bits;
  int k = 0;
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j, ++k) {
      if (bits[k] == '1') t.out[i] |= std::uint16_t(1u << j);
      else if (bits[k] == '0') t.out[j] |= std::uint16_t(1u << i);
      else throw std::runtime_error("nonbinary label");
    }
  }
  return t;
}

static Tournament from_code(int n, std::uint64_t code) {
  Tournament t;
  t.n = n;
  int k = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j, ++k) {
      bool one = ((code >> k) & 1u) != 0;
      t.label.push_back(one ? '1' : '0');
      if (one) t.out[i] |= std::uint16_t(1u << j);
      else t.out[j] |= std::uint16_t(1u << i);
    }
  }
  return t;
}

static Tournament extend(const Tournament& q, unsigned pattern) {
  Tournament t = q;
  int x = q.n;
  t.n = q.n + 1;
  for (int i = 0; i < q.n; ++i) {
    if ((pattern >> i) & 1u) t.out[x] |= std::uint16_t(1u << i);
    else t.out[i] |= std::uint16_t(1u << x);
  }
  t.label.clear();
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j) {
      t.label.push_back(((t.out[i] >> j) & 1u) ? '1' : '0');
    }
  }
  return t;
}

static Tournament induced(const Tournament& t, unsigned deleted, std::vector<int>* old_vertices = nullptr) {
  Tournament q;
  std::array<int, 12> map{};
  map.fill(-1);
  for (int i = 0; i < t.n; ++i) if (((deleted >> i) & 1u) == 0) {
    map[i] = q.n++;
    if (old_vertices) old_vertices->push_back(i);
  }
  for (int i = 0; i < t.n; ++i) if (map[i] >= 0) {
    for (int j = i + 1; j < t.n; ++j) if (map[j] >= 0) {
      if ((t.out[i] >> j) & 1u) q.out[map[i]] |= std::uint16_t(1u << map[j]);
      else q.out[map[j]] |= std::uint16_t(1u << map[i]);
    }
  }
  for (int i = 0; i < q.n; ++i) for (int j = i + 1; j < q.n; ++j)
    q.label.push_back(((q.out[i] >> j) & 1u) ? '1' : '0');
  return q;
}

static bool strong(const Tournament& t) {
  auto reach = [&](bool reverse) {
    std::uint16_t seen = 1, frontier = 1;
    while (frontier) {
      int u = __builtin_ctz(frontier);
      frontier &= std::uint16_t(frontier - 1);
      std::uint16_t nb = 0;
      if (!reverse) nb = t.out[u];
      else {
        for (int v = 0; v < t.n; ++v) if ((t.out[v] >> u) & 1u) nb |= std::uint16_t(1u << v);
      }
      nb &= std::uint16_t((1u << t.n) - 1u);
      nb &= std::uint16_t(~seen);
      seen |= nb;
      frontier |= nb;
    }
    return seen == std::uint16_t((1u << t.n) - 1u);
  };
  return reach(false) && reach(true);
}

static bool prime(const Tournament& t) {
  if (t.n < 3) return false;
  unsigned full = (1u << t.n) - 1u;
  for (unsigned m = 1; m < full; ++m) {
    int sz = __builtin_popcount(m);
    if (sz < 2 || sz == t.n) continue;
    bool module = true;
    unsigned outside = full ^ m;
    while (outside) {
      int v = __builtin_ctz(outside);
      outside &= outside - 1;
      unsigned relation = t.out[v] & m;
      if (relation != 0 && relation != m) { module = false; break; }
    }
    if (module) return false;
  }
  return true;
}

static i64 automorphism_count(const Tournament& t) {
  std::array<int, 12> deg{}, map{}, used{};
  map.fill(-1);
  for (int i = 0; i < t.n; ++i) deg[i] = __builtin_popcount(t.out[i]);
  i64 count = 0;
  auto dfs = [&](auto&& self, int depth) -> void {
    if (depth == t.n) { ++count; return; }
    int u = -1, best = 100;
    for (int x = 0; x < t.n; ++x) if (map[x] < 0) {
      int options = 0;
      for (int y = 0; y < t.n; ++y) if (!used[y] && deg[y] == deg[x]) ++options;
      if (options < best) { best = options; u = x; }
    }
    for (int v = 0; v < t.n; ++v) if (!used[v] && deg[v] == deg[u]) {
      bool ok = true;
      for (int w = 0; w < t.n; ++w) if (map[w] >= 0) {
        if (((t.out[u] >> w) & 1u) != ((t.out[v] >> map[w]) & 1u)) { ok = false; break; }
      }
      if (!ok) continue;
      map[u] = v; used[v] = 1;
      self(self, depth + 1);
      used[v] = 0; map[u] = -1;
    }
  };
  dfs(dfs, 0);
  return count;
}

static int prime_card_count(const Tournament& t) {
  int count = 0;
  for (int v = 0; v < t.n; ++v) count += prime(induced(t, 1u << v));
  return count;
}

struct Capacity {
  i64 H = 0;
  std::array<std::array<i64, 12>, 12> c{};
};

static Capacity capacities(const Tournament& t) {
  const int n = t.n;
  const int size = 1 << n;
  std::vector<std::array<i64, 12>> ending(size), starting(size), before(size), after(size);
  std::array<std::uint16_t, 12> incoming{};
  for (int u = 0; u < n; ++u) {
    unsigned vs = t.out[u];
    while (vs) {
      int v = __builtin_ctz(vs);
      vs &= vs - 1;
      incoming[v] |= std::uint16_t(1u << u);
    }
    ending[1 << u][u] = 1;
    starting[1 << u][u] = 1;
  }
  for (int mask = 1; mask < size; ++mask) {
    unsigned vs = unsigned(mask);
    while (vs) {
      int v = __builtin_ctz(vs);
      vs &= vs - 1;
      int rest = mask ^ (1 << v);
      unsigned ps = unsigned(rest) & incoming[v];
      while (ps) {
        int u = __builtin_ctz(ps);
        ps &= ps - 1;
        ending[mask][v] += ending[rest][u];
      }
      unsigned ss = unsigned(rest) & t.out[v];
      while (ss) {
        int u = __builtin_ctz(ss);
        ss &= ss - 1;
        starting[mask][v] += starting[rest][u];
      }
    }
  }
  for (int b = 0; b < n; ++b) {
    before[0][b] = after[0][b] = 1;
  }
  for (int mask = 1; mask < size; ++mask) {
    for (int b = 0; b < n; ++b) {
      unsigned vs = unsigned(mask) & incoming[b];
      while (vs) {
        int v = __builtin_ctz(vs);
        vs &= vs - 1;
        before[mask][b] += ending[mask][v];
      }
      vs = unsigned(mask) & t.out[b];
      while (vs) {
        int v = __builtin_ctz(vs);
        vs &= vs - 1;
        after[mask][b] += starting[mask][v];
      }
    }
  }
  Capacity ans;
  for (int v = 0; v < n; ++v) ans.H += ending[size - 1][v];
  unsigned full = unsigned(size - 1);
  for (int x = 0; x < n; ++x) for (int y = x + 1; y < n; ++y) {
    unsigned rem = full ^ (1u << x) ^ (1u << y);
    unsigned left = rem;
    i64 value = 0;
    for (;;) {
      unsigned right = rem ^ left;
      value += before[left][x] * after[right][y];
      value += before[left][y] * after[right][x];
      if (left == 0) break;
      left = (left - 1) & rem;
    }
    ans.c[x][y] = ans.c[y][x] = value;
  }
  return ans;
}

struct Packet { i64 C = 0, D = 0; };

static Packet packet(const Tournament& t, const Capacity& cap, unsigned deleted) {
  std::array<i64, 12> d{}, h{};
  i64 total = 0, q2 = 0;
  for (int i = 0; i < t.n; ++i) if (((deleted >> i) & 1u) == 0) {
    for (int j = i + 1; j < t.n; ++j) if (((deleted >> j) & 1u) == 0) {
      i64 z = cap.c[i][j];
      total += z;
      q2 += z * z;
      d[i] += z; d[j] += z;
      if ((t.out[i] >> j) & 1u) { h[i] += z; h[j] -= z; }
      else { h[i] -= z; h[j] += z; }
    }
  }
  Packet p;
  i64 sumd2 = 0;
  for (int i = 0; i < t.n; ++i) if (((deleted >> i) & 1u) == 0) {
    p.C += d[i] * h[i];
    sumd2 += d[i] * d[i];
  }
  i64 numer = total * total + q2 - sumd2;
  if (numer & 1) throw std::runtime_error("odd D numerator");
  p.D = numer / 2;
  return p;
}

struct Summary {
  i64 rows = 0, prime_rows = 0, strong_rows = 0;
  i64 parent_fail = 0, all_restriction_cert_fail = 0, abs_bary_cert_fail = 0;
  i64 prime_parent_fail = 0, strong_parent_fail = 0;
  i64 prime_all_restriction_cert_fail = 0, strong_all_restriction_cert_fail = 0;
  i64 prime_abs_bary_cert_fail = 0, strong_abs_bary_cert_fail = 0;
  i64 restriction_failures = 0, identity_failures = 0, zero_D = 0;
  long double max_restriction_ratio = -1, max_abs_bary_ratio = -1;
  std::string max_restriction_witness, max_abs_witness;
  std::string first_all_fail, first_abs_fail;
};

static void evaluate(const Tournament& t, bool target, Summary& s, const std::string& tag) {
  Capacity cap = capacities(t);
  Packet pp = packet(t, cap, 0);
  ++s.rows;
  bool is_prime = prime(t), is_strong = strong(t);
  if (is_prime) ++s.prime_rows;
  if (is_strong) ++s.strong_rows;
  int parent_a = t.n - 3;
  int parent_b = (t.n % 2 == 0) ? 2 : 4;
  if (i128(parent_a) * std::llabs(pp.C) >= i128(parent_b) * pp.D) {
    ++s.parent_fail;
    s.prime_parent_fail += is_prime;
    s.strong_parent_fail += is_strong;
  }

  if (t.n < 6) return;
  int k = t.n - 2;
  int a = k - 3;
  int b = (k % 2 == 0) ? 2 : 4;
  i128 sumC = 0, sumD = 0, sumAbsC = 0;
  bool all_ok = true;
  for (int u = 0; u < t.n; ++u) for (int v = u + 1; v < t.n; ++v) {
    Packet pr = packet(t, cap, (1u << u) | (1u << v));
    sumC += pr.C; sumD += pr.D; sumAbsC += std::llabs(pr.C);
    if (pr.D <= 0) ++s.zero_D;
    i128 lhs = i128(a) * std::llabs(pr.C), rhs = i128(b) * pr.D;
    long double ratio = pr.D == 0 ? std::numeric_limits<long double>::infinity()
                                         : (long double)a * std::llabs(pr.C) / ((long double)b * pr.D);
    if (ratio > s.max_restriction_ratio) {
      s.max_restriction_ratio = ratio;
      std::ostringstream os;
      os << tag << " delete=" << u << "," << v << " C=" << pr.C << " D=" << pr.D;
      s.max_restriction_witness = os.str();
    }
    if (lhs >= rhs) { all_ok = false; ++s.restriction_failures; }
  }
  i128 expectedC = i128(choose2(t.n - 3)) * pp.C;
  i128 expectedD = i128(choose2(t.n - 4)) * pp.D;
  if (sumC != expectedC || sumD != expectedD) {
    ++s.identity_failures;
    std::cerr << "identity failure " << tag << "\n";
  }
  if (!all_ok) {
    ++s.all_restriction_cert_fail;
    s.prime_all_restriction_cert_fail += is_prime;
    s.strong_all_restriction_cert_fail += is_strong;
    if (s.first_all_fail.empty()) s.first_all_fail = tag;
  }
  long double absratio = sumD == 0 ? std::numeric_limits<long double>::infinity()
      : (long double)a * (long double)sumAbsC / ((long double)b * (long double)sumD);
  if (absratio > s.max_abs_bary_ratio) {
    s.max_abs_bary_ratio = absratio;
    std::ostringstream os;
    os << tag << " sumAbsC=" << i128str(sumAbsC) << " sumD=" << i128str(sumD);
    s.max_abs_witness = os.str();
  }
  if (i128(a) * sumAbsC >= i128(b) * sumD) {
    ++s.abs_bary_cert_fail;
    s.prime_abs_bary_cert_fail += is_prime;
    s.strong_abs_bary_cert_fail += is_strong;
    if (s.first_abs_fail.empty()) s.first_abs_fail = tag;
  }
  (void)target;
}

static void print_summary(const std::string& mode, const Summary& s) {
  std::cout << "mode=" << mode
            << " rows=" << s.rows << " prime_rows=" << s.prime_rows
            << " strong_rows=" << s.strong_rows
            << " parent_fail=" << s.parent_fail
            << " prime_parent_fail=" << s.prime_parent_fail
            << " strong_parent_fail=" << s.strong_parent_fail
            << " all_restriction_cert_fail=" << s.all_restriction_cert_fail
            << " prime_all_restriction_cert_fail=" << s.prime_all_restriction_cert_fail
            << " strong_all_restriction_cert_fail=" << s.strong_all_restriction_cert_fail
            << " abs_bary_cert_fail=" << s.abs_bary_cert_fail
            << " prime_abs_bary_cert_fail=" << s.prime_abs_bary_cert_fail
            << " strong_abs_bary_cert_fail=" << s.strong_abs_bary_cert_fail
            << " restriction_failures=" << s.restriction_failures
            << " identity_failures=" << s.identity_failures
            << " zero_D=" << s.zero_D << "\n";
  std::cout << std::setprecision(18)
            << "max_restriction_ratio=" << s.max_restriction_ratio
            << " witness=" << s.max_restriction_witness << "\n"
            << "max_abs_bary_ratio=" << s.max_abs_bary_ratio
            << " witness=" << s.max_abs_witness << "\n";
  if (!s.first_all_fail.empty()) std::cout << "first_all_fail=" << s.first_all_fail << "\n";
  if (!s.first_abs_fail.empty()) std::cout << "first_abs_fail=" << s.first_abs_fail << "\n";
}

static std::uint64_t splitmix64(std::uint64_t& state) {
  std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

int main(int argc, char** argv) {
  try {
    if (argc >= 2 && std::string(argv[1]) == "--stdin") {
      Summary s;
      std::string bits;
      while (std::cin >> bits) evaluate(parse_bits(bits), true, s, bits);
      print_summary("stdin", s);
      return 0;
    }
    if (argc >= 2 && std::string(argv[1]) == "--cube") {
      Tournament q = parse_bits("111111101111111111111101111110110111110111111");
      Summary all, prime_only;
      for (unsigned z = 0; z < (1u << q.n); ++z) {
        Tournament t = extend(q, z);
        std::string tag = "pattern=" + std::to_string(z) + " label=" + t.label;
        evaluate(t, false, all, tag);
        if (prime(t)) evaluate(t, true, prime_only, tag);
      }
      print_summary("cube_all", all);
      print_summary("cube_prime", prime_only);
      for (unsigned z : {0u, 1023u}) {
        Tournament t = extend(q, z);
        Capacity fullcap = capacities(t);
        int fail = 0, fail_actual_strong = 0;
        long double max_corrected = -1;
        std::string witness;
        for (int u = 0; u < 11; ++u) for (int v = u + 1; v < 11; ++v) {
          unsigned del = (1u << u) | (1u << v);
          Packet corrected = packet(t, fullcap, del);
          if (i128(3) * std::llabs(corrected.C) < i128(2) * corrected.D) continue;
          ++fail;
          std::vector<int> old;
          Tournament card = induced(t, del, &old);
          bool card_strong = strong(card);
          fail_actual_strong += card_strong;
          if (card_strong) {
            Capacity acap = capacities(card);
            Packet actual = packet(card, acap, 0);
            long double ratio = (long double)3 * std::llabs(corrected.C) / ((long double)2 * corrected.D);
            if (ratio > max_corrected) {
              max_corrected = ratio;
              std::ostringstream os;
              os << "pattern=" << z << " pair=" << u << "," << v
                 << " corrected_C=" << corrected.C << " corrected_D=" << corrected.D
                 << " actual_C=" << actual.C << " actual_D=" << actual.D
                 << " actual_ratio=" << (long double)3 * std::llabs(actual.C) / ((long double)2 * actual.D);
              witness = os.str();
            }
          }
        }
        std::cout << "cube_boundary pattern=" << z << " failing_pairs=" << fail
                  << " failing_with_actual_strong_card=" << fail_actual_strong
                  << " max_corrected_on_strong_card=" << max_corrected
                  << " witness=" << witness << "\n";
      }
      return 0;
    }
    if (argc >= 2 && std::string(argv[1]) == "--root-hostile") {
      Tournament q = parse_bits("111111101111111111111101111110110111110111111");
      for (unsigned z : {336u, 432u, 368u, 400u}) {
        Tournament t = extend(q, z);
        std::cout << "pattern=" << z << " prime=" << prime(t)
                  << " aut=" << automorphism_count(t)
                  << " prime_cards=" << prime_card_count(t)
                  << " label=" << t.label << "\n";
      }
      return 0;
    }
    if (argc >= 2 && std::string(argv[1]) == "--thm4133") {
      Tournament t;
      t.n = 12;
      const unsigned masks[12] = {3070,3644,3704,3824,4064,4032,3970,3846,3598,1024,2049,512};
      for (int i = 0; i < 12; ++i) t.out[i] = std::uint16_t(masks[i]);
      for (int i = 0; i < 12; ++i) for (int j = i + 1; j < 12; ++j)
        t.label.push_back(((t.out[i] >> j) & 1u) ? '1' : '0');
      Summary s;
      evaluate(t, true, s, "THM-4133 label=" + t.label);
      print_summary("thm4133", s);
      Capacity fullcap = capacities(t);
      for (int u = 0; u < 12; ++u) for (int v = u + 1; v < 12; ++v) {
        if (!((u == 1 && v == 2) || (u == 9 && v == 10))) continue;
        unsigned del = (1u << u) | (1u << v);
        Packet corrected = packet(t, fullcap, del);
        if (i128(7) * std::llabs(corrected.C) < i128(2) * corrected.D) continue;
        std::vector<int> old;
        Tournament q = induced(t, del, &old);
        Capacity qcap = capacities(q);
        Packet actual = packet(q, qcap, 0);
        i64 delta_sum = 0, delta_max = 0, delta_nonzero = 0;
        for (int i = 0; i < q.n; ++i) for (int j = i + 1; j < q.n; ++j) {
          i64 delta = fullcap.c[old[i]][old[j]] - qcap.c[i][j];
          if (delta < 0) throw std::runtime_error("negative two-deletion defect");
          delta_sum += delta;
          delta_max = std::max(delta_max, delta);
          delta_nonzero += (delta != 0);
        }
        std::cout << "failing_pair=" << u << "," << v
                  << " corrected_C=" << corrected.C << " corrected_D=" << corrected.D
                  << " actual_strong=" << strong(q) << " actual_C=" << actual.C
                  << " actual_D=" << actual.D
                  << " actual_ratio=" << std::setprecision(18)
                  << (long double)7 * std::llabs(actual.C) / ((long double)2 * actual.D)
                  << " delta_nonzero=" << delta_nonzero << " delta_sum=" << delta_sum
                  << " delta_max=" << delta_max << "\n";
      }
      return 0;
    }
    if (argc >= 3 && std::string(argv[1]) == "--exhaust") {
      int n = std::atoi(argv[2]);
      int e = n * (n - 1) / 2;
      if (e >= 63) throw std::runtime_error("exhaust order too large");
      Summary all, strong_only, prime_only;
      std::uint64_t total = 1ull << e;
      for (std::uint64_t code = 0; code < total; ++code) {
        Tournament t = from_code(n, code);
        std::string tag = "code=" + std::to_string(code) + " label=" + t.label;
        evaluate(t, false, all, tag);
        if (strong(t)) evaluate(t, true, strong_only, tag);
        if (prime(t)) evaluate(t, true, prime_only, tag);
      }
      print_summary("exhaust_all_n" + std::to_string(n), all);
      print_summary("exhaust_strong_n" + std::to_string(n), strong_only);
      print_summary("exhaust_prime_n" + std::to_string(n), prime_only);
      return 0;
    }
    if (argc >= 4 && std::string(argv[1]) == "--random") {
      int n = std::atoi(argv[2]);
      i64 count = std::atoll(argv[3]);
      int e = n * (n - 1) / 2;
      if (e > 63 || n < 6 || count < 0) throw std::runtime_error("bad random request");
      const std::uint64_t seed0 = 0x4167416841694123ULL;
      std::uint64_t state = seed0;
      std::uint64_t mask = (e == 64) ? ~0ULL : ((1ULL << e) - 1ULL);
      Summary all;
      for (i64 row = 0; row < count; ++row) {
        std::uint64_t code = splitmix64(state) & mask;
        Tournament t = from_code(n, code);
        evaluate(t, false, all, "row=" + std::to_string(row) + " code=" + std::to_string(code) + " label=" + t.label);
      }
      print_summary("random_n" + std::to_string(n) + "_seed_4167416841694123", all);
      return 0;
    }
    std::cerr << "usage: --stdin | --cube | --root-hostile | --thm4133 | --exhaust n | --random n count\n";
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
}
