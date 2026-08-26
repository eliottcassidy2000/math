#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

[[noreturn]] static void die(const std::string& message) {
  std::cerr << "FAIL " << message << '\n';
  std::exit(1);
}

struct Tournament {
  int n = 0;
  std::array<std::uint16_t, 12> out{};
};

static Tournament parse_label(const std::string& label) {
  int n = 0;
  while (n * (n - 1) / 2 < static_cast<int>(label.size())) ++n;
  if (n * (n - 1) / 2 != static_cast<int>(label.size()) || n > 11)
    die("malformed upper-triangle label");
  Tournament t;
  t.n = n;
  std::size_t at = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const char bit = label[at++];
      if (bit == '1') t.out[i] |= std::uint16_t(1u << j);
      else if (bit == '0') t.out[j] |= std::uint16_t(1u << i);
      else die("non-binary upper-triangle label");
    }
  }
  return t;
}

static std::string label_of(const Tournament& t) {
  std::string answer;
  for (int i = 0; i < t.n; ++i)
    for (int j = i + 1; j < t.n; ++j)
      answer.push_back(((t.out[i] >> j) & 1u) ? '1' : '0');
  return answer;
}

static Tournament from_code(int n, u64 code) {
  Tournament t;
  t.n = n;
  int at = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j, ++at) {
      if ((code >> at) & 1u) t.out[i] |= std::uint16_t(1u << j);
      else t.out[j] |= std::uint16_t(1u << i);
    }
  }
  return t;
}

static Tournament extend(const Tournament& q, int pattern) {
  Tournament t = q;
  const int x = q.n;
  t.n = q.n + 1;
  for (int i = 0; i < q.n; ++i) {
    if ((pattern >> i) & 1) t.out[x] |= std::uint16_t(1u << i);
    else t.out[i] |= std::uint16_t(1u << x);
  }
  return t;
}

static Tournament delete_vertex(const Tournament& t, int removed) {
  Tournament card;
  card.n = t.n - 1;
  std::array<int, 12> image{};
  int next = 0;
  for (int v = 0; v < t.n; ++v) image[v] = (v == removed) ? -1 : next++;
  for (int u = 0; u < t.n; ++u) if (u != removed) {
    for (int v = 0; v < t.n; ++v) if (v != removed && ((t.out[u] >> v) & 1u))
      card.out[image[u]] |= std::uint16_t(1u << image[v]);
  }
  return card;
}

static bool strong(const Tournament& t) {
  const std::uint16_t full = std::uint16_t((1u << t.n) - 1u);
  for (int reverse = 0; reverse < 2; ++reverse) {
    std::uint16_t reached = 1u, frontier = 1u;
    while (frontier) {
      const std::uint16_t bit = std::uint16_t(frontier & std::uint16_t(-frontier));
      frontier ^= bit;
      const int u = __builtin_ctz(bit);
      std::uint16_t neighbors = 0;
      if (!reverse) neighbors = t.out[u];
      else {
        for (int v = 0; v < t.n; ++v)
          if ((t.out[v] >> u) & 1u) neighbors |= std::uint16_t(1u << v);
      }
      const std::uint16_t fresh = std::uint16_t(neighbors & full & ~reached);
      reached |= fresh;
      frontier |= fresh;
    }
    if (reached != full) return false;
  }
  return true;
}

// Slow direct definition, used as an independent control at small orders.
static bool prime_by_subsets(const Tournament& t) {
  if (t.n < 3) return false;
  const int full = (1 << t.n) - 1;
  for (int module = 1; module < full; ++module) {
    const int size = __builtin_popcount(static_cast<unsigned>(module));
    if (size < 2) continue;
    bool homogeneous = true;
    for (int v = 0; v < t.n && homogeneous; ++v) {
      if ((module >> v) & 1) continue;
      const int wins = __builtin_popcount(static_cast<unsigned>(t.out[v] & module));
      if (wins != 0 && wins != size) homogeneous = false;
    }
    if (homogeneous) return false;
  }
  return true;
}

// Pair-closure criterion, used for the 9.7-million-class stream.
static bool prime_by_pair_closure(const Tournament& t) {
  if (t.n < 3) return false;
  const std::uint16_t full = std::uint16_t((1u << t.n) - 1u);
  for (int a = 0; a < t.n; ++a) {
    for (int b = a + 1; b < t.n; ++b) {
      std::uint16_t module = std::uint16_t((1u << a) | (1u << b));
      while (module != full) {
        std::uint16_t add = 0;
        std::uint16_t outside = std::uint16_t(full ^ module);
        while (outside) {
          const std::uint16_t bit = std::uint16_t(outside & std::uint16_t(-outside));
          outside ^= bit;
          const int v = __builtin_ctz(bit);
          const std::uint16_t rel = std::uint16_t(t.out[v] & module);
          if (rel != 0 && rel != module) add |= bit;
        }
        if (!add) return false;
        module |= add;
      }
    }
  }
  return true;
}

static u64 hamiltonian_paths(const Tournament& t) {
  const int size = 1 << t.n;
  std::vector<std::array<u64, 12>> dp(size);
  for (int v = 0; v < t.n; ++v) dp[1 << v][v] = 1;
  for (int mask = 1; mask < size; ++mask) {
    int vertices = mask;
    while (vertices) {
      const int bit = vertices & -vertices;
      vertices ^= bit;
      const int v = __builtin_ctz(static_cast<unsigned>(bit));
      const int rest = mask ^ bit;
      int pred = rest;
      while (pred) {
        const int pbit = pred & -pred;
        pred ^= pbit;
        const int p = __builtin_ctz(static_cast<unsigned>(pbit));
        if ((t.out[p] >> v) & 1u) dp[mask][v] += dp[rest][p];
      }
    }
  }
  u64 answer = 0;
  for (int v = 0; v < t.n; ++v) answer += dp[size - 1][v];
  return answer;
}

struct CapacityPacket {
  u64 H = 0;
  std::array<std::array<i64, 12>, 12> cap{};
  i64 C = 0;
  i64 D = 0;
};

// This deliberately avoids exposed-block endpoint convolution.  It computes
// the coefficient from four literal one-ear children: F(S)=H(T+x_S), where
// x_S -> v iff v is in S.  Since x_empty is a sink, F(empty)=H(T).
static CapacityPacket literal_ear_packet(const Tournament& t) {
  CapacityPacket p;
  p.H = hamiltonian_paths(t);
  std::array<u64, 12> singleton{};
  for (int i = 0; i < t.n; ++i)
    singleton[i] = hamiltonian_paths(extend(t, 1 << i));
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j) {
      const u64 pair = hamiltonian_paths(extend(t, (1 << i) | (1 << j)));
      const i64 c = i64(singleton[i]) + i64(singleton[j]) - i64(p.H) - i64(pair);
      // At order three a transitive edge can have zero exposed capacity.
      // The positive-capacity theorem used for Johnson denominators starts at
      // order four; negativity would signal an ear-coordinate/sign error.
      if (c < 0) die("negative direct-switch capacity");
      p.cap[i][j] = p.cap[j][i] = c;
    }
  }
  std::array<i64, 12> degree{}, signed_degree{};
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j) {
      degree[i] += p.cap[i][j];
      degree[j] += p.cap[i][j];
      const i64 sign = ((t.out[i] >> j) & 1u) ? 1 : -1;
      signed_degree[i] += sign * p.cap[i][j];
      signed_degree[j] -= sign * p.cap[i][j];
    }
  }
  for (int i = 0; i < t.n; ++i) p.C += degree[i] * signed_degree[i];
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j) {
      for (int u = 0; u < t.n; ++u) {
        for (int v = u + 1; v < t.n; ++v) {
          if (std::tie(i, j) >= std::tie(u, v)) continue;
          if (i != u && i != v && j != u && j != v)
            p.D += p.cap[i][j] * p.cap[u][v];
        }
      }
    }
  }
  return p;
}

static i64 ceil_ratio(i128 numerator, i128 denominator) {
  if (denominator <= 0) die("bad ceil denominator");
  if (numerator >= 0) return i64((numerator + denominator - 1) / denominator);
  return i64(numerator / denominator);
}

struct ProfileSummary {
  i64 central_coset = std::numeric_limits<i64>::min();
  i64 best_outer_coset = std::numeric_limits<i64>::min();
  i64 central_actual = std::numeric_limits<i64>::min();
  i64 best_outer_actual = std::numeric_limits<i64>::min();
  std::vector<int> rational_t, coset_t, actual_t;
};

static ProfileSummary exact_profile(const Tournament& t, const CapacityPacket& p) {
  const int size = 1 << t.n;
  std::vector<i64> field(size, i64(p.H));
  for (int mask = 0; mask < size; ++mask) {
    for (int i = 0; i < t.n; ++i) if ((mask >> i) & 1) {
      for (int j = 0; j < t.n; ++j) if (!((mask >> j) & 1)) {
        if ((t.out[i] >> j) & 1u) field[mask] += p.cap[i][j];
      }
    }
  }
  struct Layer { int imbalance; i128 jnum, jden; i64 coset, actual; };
  std::vector<Layer> layers;
  for (int m = 1; m < t.n; ++m) {
    std::vector<i64> values;
    for (int mask = 0; mask < size; ++mask)
      if (__builtin_popcount(static_cast<unsigned>(mask)) == m) values.push_back(field[mask]);
    i128 total = 0, square = 0;
    for (i64 value : values) { total += value; square += i128(value) * value; }
    const i128 jnum = square - i128(p.H) * total;
    const i128 jden = total - i128(values.size()) * p.H;
    if (jden <= 0) die("nonpositive Johnson denominator");
    const i64 anchor = values.front();
    i64 lattice = 0;
    for (i64 value : values) lattice = std::gcd(lattice, std::llabs(value - anchor));
    i64 coset = anchor;
    if (lattice) {
      const i64 k = ceil_ratio(jnum - i128(anchor) * jden, jden * lattice);
      coset += k * lattice;
    }
    const i64 actual = *std::max_element(values.begin(), values.end());
    layers.push_back({t.n - 2 * m, jnum, jden, coset, actual});
  }
  auto better_j = [](const Layer& a, const Layer& b) {
    return a.jnum * b.jden > b.jnum * a.jden;
  };
  const Layer* bestj = &layers.front();
  i64 bestc = layers.front().coset, besta = layers.front().actual;
  for (const Layer& layer : layers) {
    if (better_j(layer, *bestj)) bestj = &layer;
    bestc = std::max(bestc, layer.coset);
    besta = std::max(besta, layer.actual);
  }
  ProfileSummary s;
  for (const Layer& layer : layers) {
    if (layer.jnum * bestj->jden == bestj->jnum * layer.jden) s.rational_t.push_back(layer.imbalance);
    if (layer.coset == bestc) s.coset_t.push_back(layer.imbalance);
    if (layer.actual == besta) s.actual_t.push_back(layer.imbalance);
    const bool central = std::abs(layer.imbalance) == 1;
    if (central) {
      s.central_coset = std::max(s.central_coset, layer.coset);
      s.central_actual = std::max(s.central_actual, layer.actual);
    } else {
      s.best_outer_coset = std::max(s.best_outer_coset, layer.coset);
      s.best_outer_actual = std::max(s.best_outer_actual, layer.actual);
    }
  }
  return s;
}

static std::vector<i64> layer_lattices(const Tournament& t, const CapacityPacket& p) {
  const int size = 1 << t.n;
  std::vector<i64> field(size, i64(p.H));
  for (int mask = 0; mask < size; ++mask) {
    for (int i = 0; i < t.n; ++i) if ((mask >> i) & 1) {
      for (int j = 0; j < t.n; ++j) if (!((mask >> j) & 1) && ((t.out[i] >> j) & 1u))
        field[mask] += p.cap[i][j];
    }
  }
  std::vector<i64> answer;
  for (int m = 1; m < t.n; ++m) {
    i64 anchor = -1, lattice = 0;
    for (int mask = 0; mask < size; ++mask) {
      if (__builtin_popcount(static_cast<unsigned>(mask)) != m) continue;
      if (anchor < 0) anchor = field[mask];
      else lattice = std::gcd(lattice, std::llabs(field[mask] - anchor));
    }
    answer.push_back(lattice);
  }
  return answer;
}

static Tournament critical_t(int n) {
  Tournament t;
  t.n = n;
  const int k = (n - 1) / 2;
  for (int i = 0; i < n; ++i)
    for (int d = 1; d <= k; ++d) t.out[i] |= std::uint16_t(1u << ((i + d) % n));
  return t;
}

static Tournament critical_u(int n) {
  Tournament t = critical_t(n);
  const int k = (n - 1) / 2;
  for (int a = k + 1; a < n; ++a) {
    for (int b = a + 1; b < n; ++b) {
      if ((t.out[a] >> b) & 1u) {
        t.out[a] ^= std::uint16_t(1u << b);
        t.out[b] |= std::uint16_t(1u << a);
      } else {
        t.out[b] ^= std::uint16_t(1u << a);
        t.out[a] |= std::uint16_t(1u << b);
      }
    }
  }
  return t;
}

static Tournament critical_w(int n) {
  Tournament t;
  t.n = n;
  const int x = n - 1;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n - 1; ++j) t.out[i] |= std::uint16_t(1u << j);
    if ((i & 1) == 0) t.out[x] |= std::uint16_t(1u << i);
    else t.out[i] |= std::uint16_t(1u << x);
  }
  return t;
}

static int clone_vertices(const Tournament& q, int pattern) {
  int count = 0;
  for (int v = 0; v < q.n; ++v) {
    bool clone = true;
    for (int y = 0; y < q.n; ++y) {
      if (y == v) continue;
      const bool xy = ((pattern >> y) & 1) != 0;
      const bool vy = ((q.out[v] >> y) & 1u) != 0;
      if (xy != vy) { clone = false; break; }
    }
    count += clone;
  }
  return count;
}

struct ExtensionSummary { int nonstrong = 0, decomposable_strong = 0, prime = 0; };

static ExtensionSummary classify_extensions(const Tournament& q, bool direct_prime) {
  ExtensionSummary s;
  const int all = (1 << q.n) - 1;
  for (int pattern = 0; pattern <= all; ++pattern) {
    const Tournament t = extend(q, pattern);
    const bool st = strong(t);
    const bool pr = direct_prime ? prime_by_subsets(t) : prime_by_pair_closure(t);
    if (!st) {
      ++s.nonstrong;
      if (pr || (pattern != 0 && pattern != all)) die("unexpected nonstrong extension type");
    } else if (!pr) {
      ++s.decomposable_strong;
      if (clone_vertices(q, pattern) != 1) die("strong decomposable extension is not a unique clone pattern");
    } else {
      ++s.prime;
      if (clone_vertices(q, pattern) != 0) die("prime extension has clone module");
    }
  }
  return s;
}

struct AutomorphismData {
  u64 count = 0;
  i128 burnside_sum = 0;
};

static AutomorphismData automorphisms(const Tournament& t) {
  std::array<int, 12> degree{}, map{}, inverse{};
  map.fill(-1);
  inverse.fill(-1);
  for (int i = 0; i < t.n; ++i) degree[i] = __builtin_popcount(t.out[i]);
  AutomorphismData answer;
  auto dfs = [&](auto&& self, int assigned) -> void {
    if (assigned == t.n) {
      ++answer.count;
      int cycles = 0, fixed = 0;
      std::array<bool, 12> seen{};
      for (int i = 0; i < t.n; ++i) {
        fixed += map[i] == i;
        if (!seen[i]) {
          ++cycles;
          int u = i;
          while (!seen[u]) { seen[u] = true; u = map[u]; }
        }
      }
      answer.burnside_sum += (i128(1) << cycles) - 2 - 2 * fixed;
      return;
    }
    int chosen = -1;
    std::vector<int> chosen_candidates;
    for (int u = 0; u < t.n; ++u) if (map[u] < 0) {
      std::vector<int> candidates;
      for (int v = 0; v < t.n; ++v) if (inverse[v] < 0 && degree[u] == degree[v]) {
        bool ok = true;
        for (int w = 0; w < t.n && ok; ++w) if (map[w] >= 0) {
          const bool uw = ((t.out[u] >> w) & 1u) != 0;
          const bool image = ((t.out[v] >> map[w]) & 1u) != 0;
          ok = uw == image;
        }
        if (ok) candidates.push_back(v);
      }
      if (chosen < 0 || candidates.size() < chosen_candidates.size()) {
        chosen = u;
        chosen_candidates = std::move(candidates);
      }
    }
    for (int v : chosen_candidates) {
      map[chosen] = v;
      inverse[v] = chosen;
      self(self, assigned + 1);
      inverse[v] = -1;
      map[chosen] = -1;
    }
  };
  dfs(dfs, 0);
  return answer;
}

using Polynomial = std::vector<i64>;

static Polynomial poly_add(const Polynomial& a, const Polynomial& b, i64 scale_b = 1) {
  Polynomial c(a.size());
  for (std::size_t i = 0; i < a.size(); ++i) c[i] = a[i] + scale_b * b[i];
  return c;
}

static Polynomial poly_mul(const Polynomial& a, const Polynomial& b, int max_degree) {
  Polynomial c(a.size());
  for (std::size_t x = 0; x < a.size(); ++x) if (a[x]) {
    for (std::size_t y = 0; y < b.size(); ++y) if (b[y]) {
      const std::size_t z = x | y;
      if (__builtin_popcount(static_cast<unsigned>(z)) <= max_degree) c[z] += a[x] * b[y];
    }
  }
  return c;
}

static int polynomial_degree(const Polynomial& p) {
  int degree = -1;
  for (std::size_t m = 0; m < p.size(); ++m) if (p[m])
    degree = std::max(degree, __builtin_popcount(static_cast<unsigned>(m)));
  return degree;
}

static std::vector<i64> zeta_values(const Polynomial& p, int q) {
  std::vector<i64> values = p;
  for (int i = 0; i < q; ++i)
    for (int mask = 0; mask < (1 << q); ++mask)
      if ((mask >> i) & 1) values[mask] += values[mask ^ (1 << i)];
  return values;
}

struct PolynomialAuditResult {
  int max_H_degree = -1, max_old_degree = -1, max_new_degree = -1;
  int max_packet_degree = -1;
  int literal_patterns = 0;
  int strong_patterns = 0, strong_rational_failures = 0;
  int H_terms = 0, C_terms = 0, D_terms = 0, plus_terms = 0, minus_terms = 0;
};

// Interpolate H and every capacity from literal values at weight <=2.  The
// general degree theorem is proved separately; higher-weight literal patterns
// are independent checks, not assumptions injected into the construction.
static PolynomialAuditResult polynomial_audit(const Tournament& q,
                                               const std::vector<int>& check_patterns) {
  const int cube = 1 << q.n;
  const int n = q.n + 1;
  std::vector<std::pair<int,int>> edges;
  for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) edges.emplace_back(i,j);
  std::vector<int> interpolation_patterns{0};
  for (int i = 0; i < q.n; ++i) interpolation_patterns.push_back(1 << i);
  for (int i = 0; i < q.n; ++i)
    for (int j = i + 1; j < q.n; ++j) interpolation_patterns.push_back((1 << i) | (1 << j));
  std::vector<i64> H_values(cube);
  std::vector<std::vector<i64>> values(edges.size(), std::vector<i64>(cube));
  for (int pattern : interpolation_patterns) {
    const CapacityPacket packet = literal_ear_packet(extend(q, pattern));
    H_values[pattern] = i64(packet.H);
    for (std::size_t e = 0; e < edges.size(); ++e)
      values[e][pattern] = packet.cap[edges[e].first][edges[e].second];
  }
  Polynomial H(cube);
  H[0] = H_values[0];
  for (int i = 0; i < q.n; ++i) H[1 << i] = H_values[1 << i] - H_values[0];
  for (int i = 0; i < q.n; ++i) for (int j = i + 1; j < q.n; ++j) {
    const int ij = (1 << i) | (1 << j);
    H[ij] = H_values[ij] - H_values[1 << i] - H_values[1 << j] + H_values[0];
  }
  const CapacityPacket parent_packet = literal_ear_packet(q);
  Polynomial H_formula(cube);
  H_formula[0] = i64(parent_packet.H);
  for (int i = 0; i < q.n; ++i) for (int j = i + 1; j < q.n; ++j) {
    const bool i_to_j = ((q.out[i] >> j) & 1u) != 0;
    const int tail = i_to_j ? i : j;
    H_formula[1 << tail] += parent_packet.cap[i][j];
    H_formula[(1 << i) | (1 << j)] -= parent_packet.cap[i][j];
  }
  if (H != H_formula) die("interpolated H disagrees with the parent cut formula");
  std::vector<Polynomial> cap(edges.size(), Polynomial(cube));
  for (std::size_t e = 0; e < edges.size(); ++e) {
    cap[e][0] = values[e][0];
    for (int i = 0; i < q.n; ++i) cap[e][1 << i] = values[e][1 << i] - values[e][0];
    for (int i = 0; i < q.n; ++i) for (int j = i + 1; j < q.n; ++j) {
      const int ij = (1 << i) | (1 << j);
      cap[e][ij] = values[e][ij] - values[e][1 << i] - values[e][1 << j] + values[e][0];
    }
  }
  PolynomialAuditResult result;
  result.max_H_degree = polynomial_degree(H);
  std::vector<Polynomial> degree(n, Polynomial(cube)), signed_degree(n, Polynomial(cube));
  for (std::size_t e = 0; e < edges.size(); ++e) {
    const auto [i,j] = edges[e];
    degree[i] = poly_add(degree[i], cap[e]);
    degree[j] = poly_add(degree[j], cap[e]);
    if (j == q.n) {
      if (cap[e][1 << i] != 0) die("new-edge capacity depends on its mutual arc bit");
      for (int a = 0; a < q.n; ++a) for (int b = a + 1; b < q.n; ++b)
        if (cap[e][(1 << a) | (1 << b)] != 0) die("new-edge capacity is not affine");
      result.max_new_degree = std::max(result.max_new_degree, polynomial_degree(cap[e]));
      Polynomial sign(cube); sign[0] = 1; sign[1 << i] = -2;
      Polynomial term = poly_mul(sign, cap[e], 2);
      signed_degree[i] = poly_add(signed_degree[i], term);
      signed_degree[j] = poly_add(signed_degree[j], term, -1);
    } else {
      result.max_old_degree = std::max(result.max_old_degree, polynomial_degree(cap[e]));
      const i64 sign = ((q.out[i] >> j) & 1u) ? 1 : -1;
      signed_degree[i] = poly_add(signed_degree[i], cap[e], sign);
      signed_degree[j] = poly_add(signed_degree[j], cap[e], -sign);
    }
  }
  Polynomial C(cube), D(cube);
  for (int i = 0; i < n; ++i) C = poly_add(C, poly_mul(degree[i], signed_degree[i], 4));
  for (std::size_t a = 0; a < edges.size(); ++a) {
    for (std::size_t b = a + 1; b < edges.size(); ++b) {
      const auto [i,j] = edges[a]; const auto [u,v] = edges[b];
      if (i != u && i != v && j != u && j != v)
        D = poly_add(D, poly_mul(cap[a], cap[b], 4));
    }
  }
  const Polynomial plus = poly_add(D, C, 2), minus = poly_add(D, C, -2);
  result.max_packet_degree = std::max({polynomial_degree(C), polynomial_degree(D),
                                       polynomial_degree(plus), polynomial_degree(minus)});
  auto nonzero = [](const Polynomial& p) {
    return int(std::count_if(p.begin(), p.end(), [](i64 x) { return x != 0; }));
  };
  result.H_terms = nonzero(H); result.C_terms = nonzero(C); result.D_terms = nonzero(D);
  result.plus_terms = nonzero(plus); result.minus_terms = nonzero(minus);
  const auto evaluated_H = zeta_values(H, q.n);
  const auto Cvalues = zeta_values(C, q.n), Dvalues = zeta_values(D, q.n);
  for (int pattern = 0; pattern < cube; ++pattern) {
    if (!strong(extend(q, pattern))) continue;
    ++result.strong_patterns;
    result.strong_rational_failures +=
        Dvalues[pattern] <= 2 * std::llabs(Cvalues[pattern]);
  }
  std::vector<std::vector<i64>> cap_values;
  cap_values.reserve(cap.size());
  for (const Polynomial& p : cap) cap_values.push_back(zeta_values(p, q.n));
  for (int pattern : check_patterns) {
    if (pattern < 0 || pattern >= cube) die("bad polynomial check pattern");
    const CapacityPacket packet = literal_ear_packet(extend(q, pattern));
    if (evaluated_H[pattern] != i64(packet.H))
      die("higher-weight literal H disagrees with quadratic zeta evaluation");
    for (std::size_t e = 0; e < edges.size(); ++e) {
      if (cap_values[e][pattern] != packet.cap[edges[e].first][edges[e].second])
        die("higher-weight literal capacity disagrees with degree-two interpolation");
    }
    if (Cvalues[pattern] != packet.C || Dvalues[pattern] != packet.D)
      die("higher-weight literal packet disagrees with quartic zeta evaluation");
    ++result.literal_patterns;
  }
  return result;
}

static std::string d6_of(const Tournament& t) {
  if (t.n > 62) die("long digraph6 order unsupported");
  std::string bits;
  bits.reserve(t.n * t.n);
  for (int i = 0; i < t.n; ++i)
    for (int j = 0; j < t.n; ++j)
      bits.push_back((i != j && ((t.out[i] >> j) & 1u)) ? '1' : '0');
  while (bits.size() % 6) bits.push_back('0');
  std::string out;
  out.push_back('&'); out.push_back(char(t.n + 63));
  for (std::size_t at = 0; at < bits.size(); at += 6) {
    int value = 0;
    for (int k = 0; k < 6; ++k) value = 2 * value + (bits[at + k] - '0');
    out.push_back(char(value + 63));
  }
  return out;
}

static void run_selftest() {
  u64 small_prime_parents = 0, small_patterns = 0;
  for (int q = 3; q <= 5; ++q) {
    const int edges = q * (q - 1) / 2;
    u64 parents = 0;
    for (u64 code = 0; code < (u64(1) << edges); ++code) {
      const Tournament parent = from_code(q, code);
      const bool a = prime_by_subsets(parent), b = prime_by_pair_closure(parent);
      if (a != b) die("prime checkers disagree in small exhaustive universe");
      if (!a) continue;
      ++parents;
      const ExtensionSummary s = classify_extensions(parent, true);
      if (s.nonstrong != 2 || s.decomposable_strong != 2*q ||
          s.prime != (1 << q) - 2 - 2*q) die("small extension classification count mismatch");
      small_patterns += 1u << q;
    }
    small_prime_parents += parents;
    std::cout << "small_extension_q=" << q << " prime_labelled_parents=" << parents
              << " status=PASS\n";
  }
  std::cout << "small_extension_total_parents=" << small_prime_parents
            << " patterns=" << small_patterns << '\n';

  u64 q6_prime = 0;
  for (u64 code = 0; code < (u64(1) << 15); ++code) {
    const Tournament t = from_code(6, code);
    const bool direct = prime_by_subsets(t), closure = prime_by_pair_closure(t);
    if (direct != closure) die("prime checkers disagree at order six");
    q6_prime += direct;
  }
  std::cout << "prime_checker_crosscheck_q6_labelled=32768 prime=" << q6_prime
            << " status=PASS\n";

  u64 polynomial_parents = 0, polynomial_checks = 0;
  int max_H = -1, max_old = -1, max_new = -1, max_packet = -1;
  for (int q = 2; q <= 4; ++q) {
    const int edges = q * (q - 1) / 2;
    for (u64 code = 0; code < (u64(1) << edges); ++code) {
      const Tournament parent = from_code(q, code);
      std::vector<int> checks;
      const int full = (1 << q) - 1;
      for (int pattern = 0; pattern <= full; ++pattern)
        if (__builtin_popcount(static_cast<unsigned>(pattern)) >= 3) checks.push_back(pattern);
      const auto result = polynomial_audit(parent, checks);
      max_H = std::max(max_H, result.max_H_degree);
      max_old = std::max(max_old, result.max_old_degree);
      max_new = std::max(max_new, result.max_new_degree);
      max_packet = std::max(max_packet, result.max_packet_degree);
      polynomial_checks += result.literal_patterns;
      ++polynomial_parents;
    }
  }
  std::cout << "small_polynomial_parents=" << polynomial_parents
            << " higher_weight_literal_checks=" << polynomial_checks
            << " max_H_degree=" << max_H
            << " max_old_capacity_degree=" << max_old
            << " max_new_capacity_degree=" << max_new
            << " max_packet_degree=" << max_packet << " status=PASS\n";

  for (const auto& item : std::vector<std::pair<std::string,Tournament>>{
           {"T11",critical_t(11)}, {"U11",critical_u(11)}, {"W11",critical_w(11)}}) {
    const CapacityPacket packet = literal_ear_packet(item.second);
    const ProfileSummary profile = exact_profile(item.second, packet);
    if (!prime_by_subsets(item.second) || !strong(item.second)) die("critical family structural failure");
    int prime_cards = 0, strong_cards = 0;
    for (int v = 0; v < item.second.n; ++v) {
      const Tournament card = delete_vertex(item.second, v);
      prime_cards += prime_by_subsets(card);
      strong_cards += strong(card);
    }
    std::cout << item.first << " label=" << label_of(item.second)
              << " H=" << packet.H << " C=" << packet.C << " D=" << packet.D
              << " prime_cards=" << prime_cards << " strong_cards=" << strong_cards
              << " coset_margin=" << (profile.central_coset - profile.best_outer_coset)
              << " actual_margin=" << (profile.central_actual - profile.best_outer_actual)
              << " rational_t=";
    for (int x : profile.rational_t) std::cout << x << ',';
    std::cout << " coset_t=";
    for (int x : profile.coset_t) std::cout << x << ',';
    std::cout << " actual_t=";
    for (int x : profile.actual_t) std::cout << x << ',';
    std::cout << " status=PASS\n";
  }

  const std::string hostile_label = "111111101111111111111101111110110111110111111";
  const Tournament hostile = parse_label(hostile_label);
  if (hostile.n != 10 || !prime_by_subsets(hostile) || !strong(hostile))
    die("asymmetric hostile is not a prime strong order-ten tournament");
  const AutomorphismData aut = automorphisms(hostile);
  if (aut.count == 0 || aut.burnside_sum % aut.count != 0) die("bad Burnside arithmetic");
  const ExtensionSummary ext = classify_extensions(hostile, true);
  std::cout << "asymmetric_parent label=" << hostile_label
            << " aut=" << aut.count
            << " burnside_rooted_prime_orbits=" << static_cast<long long>(aut.burnside_sum / aut.count)
            << " nonstrong=" << ext.nonstrong
            << " strong_decomposable=" << ext.decomposable_strong
            << " prime_patterns=" << ext.prime << " status=PASS\n";

  const Tournament lattice_child = extend(hostile, 308);
  const CapacityPacket lattice_packet = literal_ear_packet(lattice_child);
  const std::vector<i64> lattices = layer_lattices(lattice_child, lattice_packet);
  if (!prime_by_subsets(lattice_child) || lattices != std::vector<i64>({4,2,2,2,2,2,2,2,2,2}))
    die("pattern-308 prime nonparity hostile mismatch");
  std::cout << "prime_nonparity_hostile pattern=308 layer_lattices=";
  for (i64 value : lattices) std::cout << value << ',';
  std::cout << " status=PASS\n";

  const CapacityPacket sink_packet = literal_ear_packet(extend(hostile, 0));
  const CapacityPacket source_packet = literal_ear_packet(extend(hostile, 1023));
  const i64 sink_margin = sink_packet.D - 2 * std::llabs(sink_packet.C);
  const i64 source_twice_abs_c = 2 * std::llabs(source_packet.C);
  const i64 source_margin = source_packet.D - source_twice_abs_c;
  if (strong(extend(hostile, 0)) || strong(extend(hostile, 1023)) ||
      sink_margin != 209423320 || source_twice_abs_c != 3818407904LL ||
      source_packet.D != 2735733720LL || source_margin != -1082674184LL)
    die("uniform nonstrong attachment boundary mismatch");
  std::cout << "uniform_boundary sink_pattern=0 sink_margin=" << sink_margin
            << " source_pattern=1023 source_twice_abs_C=" << source_twice_abs_c
            << " source_D=" << source_packet.D
            << " source_margin=" << source_margin << " status=PASS\n";

  const std::vector<int> checks{7, 31, 63, 127, 255, 341, 511, 682, 767, 895, 1023};
  const auto poly = polynomial_audit(hostile, checks);
  std::cout << "q10_polynomial literal_higher_weight_checks=" << poly.literal_patterns
            << " max_H_degree=" << poly.max_H_degree
            << " max_old_capacity_degree=" << poly.max_old_degree
            << " max_new_capacity_degree=" << poly.max_new_degree
            << " max_packet_degree=" << poly.max_packet_degree
            << " terms_H_C_D_plus_minus=" << poly.H_terms << ',' << poly.C_terms << ',' << poly.D_terms << ','
            << poly.plus_terms << ',' << poly.minus_terms
            << " strong_patterns=" << poly.strong_patterns
            << " strong_rational_failures=" << poly.strong_rational_failures
            << " basis_bound=" << (1 + 10 + 45 + 120 + 210) << " status=PASS\n";
  std::cout << "status=PASS\n";
}

int main(int argc, char** argv) {
  if (argc == 1 || (argc == 2 && std::string(argv[1]) == "--selftest")) {
    run_selftest();
    return 0;
  }
  if (argc == 3 && std::string(argv[1]) == "--stream-prime") {
    const int order = std::atoi(argv[2]);
    u64 total = 0, prime = 0;
    std::string label;
    while (std::cin >> label) {
      const Tournament t = parse_label(label);
      if (t.n != order) die("stream order mismatch");
      ++total;
      prime += prime_by_pair_closure(t);
    }
    std::cout << "order=" << order << " classes=" << total << " prime_classes=" << prime << '\n';
    return 0;
  }
  if (argc == 3 && std::string(argv[1]) == "--emit-prime-extensions") {
    const Tournament q = parse_label(argv[2]);
    int count = 0;
    for (int pattern = 0; pattern < (1 << q.n); ++pattern) {
      const Tournament t = extend(q, pattern);
      if (prime_by_subsets(t)) { std::cout << d6_of(t) << '\n'; ++count; }
    }
    if (count != (1 << q.n) - 2 - 2*q.n) die("prime extension emission count mismatch");
    return 0;
  }
  if (argc == 3 && std::string(argv[1]) == "--list-prime-patterns") {
    const Tournament q = parse_label(argv[2]);
    int index = 0;
    for (int pattern = 0; pattern < (1 << q.n); ++pattern) {
      if (prime_by_subsets(extend(q, pattern)))
        std::cout << ++index << ' ' << pattern << '\n';
    }
    if (index != (1 << q.n) - 2 - 2*q.n) die("prime pattern list count mismatch");
    return 0;
  }
  if (argc == 3 && std::string(argv[1]) == "--full-q10-poly") {
    const Tournament q = parse_label(argv[2]);
    if (q.n != 10) die("full polynomial audit requires an order-ten parent");
    std::vector<int> checks;
    for (int pattern = 0; pattern < (1 << 10); ++pattern)
      if (__builtin_popcount(static_cast<unsigned>(pattern)) >= 3) checks.push_back(pattern);
    const auto result = polynomial_audit(q, checks);
    std::cout << "q10_full_polynomial parent=" << argv[2]
              << " higher_weight_literal_checks=" << result.literal_patterns
              << " max_H_degree=" << result.max_H_degree
              << " max_old_capacity_degree=" << result.max_old_degree
              << " max_new_capacity_degree=" << result.max_new_degree
              << " max_packet_degree=" << result.max_packet_degree
              << " terms_H_C_D_plus_minus=" << result.H_terms << ',' << result.C_terms << ',' << result.D_terms << ','
              << result.plus_terms << ',' << result.minus_terms
              << " strong_patterns=" << result.strong_patterns
              << " strong_rational_failures=" << result.strong_rational_failures
              << " status=PASS\n";
    return 0;
  }
  std::cerr << "usage: audit_thm4169_independent [--selftest]\n"
            << "       audit_thm4169_independent --stream-prime ORDER\n"
            << "       audit_thm4169_independent --emit-prime-extensions LABEL\n"
            << "       audit_thm4169_independent --list-prime-patterns LABEL\n"
            << "       audit_thm4169_independent --full-q10-poly LABEL\n";
  return 2;
}
