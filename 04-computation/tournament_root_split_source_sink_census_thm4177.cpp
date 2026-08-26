#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using i64 = long long;
using u16 = std::uint16_t;

struct Tournament {
  int n = 0;
  std::array<u16, 12> out{};
  std::string label;
};

struct Exposure {
  i64 H = 0;
  std::array<std::array<i64, 12>, 12> c{};
  std::array<std::array<i64, 12>, 12> Q{};
  std::array<i64, 12> start{};
  std::array<i64, 12> end{};
};

struct Packet { i64 C = 0, D = 0; };

static Tournament parse_bits(const std::string& bits) {
  Tournament t;
  const int edges = static_cast<int>(bits.size());
  while (t.n * (t.n - 1) / 2 < edges) ++t.n;
  if (t.n * (t.n - 1) / 2 != edges || t.n > 11) {
    throw std::runtime_error("nontriangular or oversized tournament label");
  }
  t.label = bits;
  int k = 0;
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j, ++k) {
      if (bits[k] == '1') t.out[i] |= static_cast<u16>(1u << j);
      else if (bits[k] == '0') t.out[j] |= static_cast<u16>(1u << i);
      else throw std::runtime_error("nonbinary tournament label");
    }
  }
  return t;
}

static Tournament extend(const Tournament& q, unsigned z) {
  Tournament t = q;
  const int x = q.n;
  ++t.n;
  for (int i = 0; i < q.n; ++i) {
    if ((z >> i) & 1u) t.out[x] |= static_cast<u16>(1u << i);
    else t.out[i] |= static_cast<u16>(1u << x);
  }
  t.label.clear();
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j) {
      t.label.push_back(((t.out[i] >> j) & 1u) ? '1' : '0');
    }
  }
  return t;
}

static bool has_sink(const Tournament& t) {
  for (int i = 0; i < t.n; ++i) if (t.out[i] == 0) return true;
  return false;
}

static bool has_source(const Tournament& t) {
  for (int i = 0; i < t.n; ++i) {
    if (__builtin_popcount(t.out[i]) == t.n - 1) return true;
  }
  return false;
}

static bool strong(const Tournament& t) {
  auto reaches_all = [&](bool reverse) {
    u16 seen = 1, frontier = 1;
    while (frontier) {
      const int u = __builtin_ctz(frontier);
      frontier &= static_cast<u16>(frontier - 1);
      u16 next = 0;
      if (!reverse) next = t.out[u];
      else {
        for (int v = 0; v < t.n; ++v) {
          if ((t.out[v] >> u) & 1u) next |= static_cast<u16>(1u << v);
        }
      }
      next &= static_cast<u16>((1u << t.n) - 1u);
      next &= static_cast<u16>(~seen);
      seen |= next;
      frontier |= next;
    }
    return seen == static_cast<u16>((1u << t.n) - 1u);
  };
  return reaches_all(false) && reaches_all(true);
}

static int minimum_module_size(const Tournament& t) {
  const unsigned full = (1u << t.n) - 1u;
  int answer = t.n + 1;
  for (unsigned m = 1; m < full; ++m) {
    const int size = __builtin_popcount(m);
    if (size < 2 || size >= answer) continue;
    bool module = true;
    unsigned outside = full ^ m;
    while (outside) {
      const int v = __builtin_ctz(outside);
      outside &= outside - 1;
      const unsigned relation = t.out[v] & m;
      if (relation != 0 && relation != m) { module = false; break; }
    }
    if (module) answer = size;
  }
  return answer == t.n + 1 ? 0 : answer;
}

static bool prime(const Tournament& t) {
  return t.n >= 3 && minimum_module_size(t) == 0;
}

static Exposure exposures(const Tournament& t) {
  const int n = t.n;
  const int count = 1 << n;
  std::vector<std::array<i64, 12>> ending(count), starting(count);
  std::vector<std::array<i64, 12>> before(count), after(count);
  std::array<u16, 12> incoming{};
  for (int u = 0; u < n; ++u) {
    unsigned vs = t.out[u];
    while (vs) {
      const int v = __builtin_ctz(vs);
      vs &= vs - 1;
      incoming[v] |= static_cast<u16>(1u << u);
    }
    ending[1 << u][u] = 1;
    starting[1 << u][u] = 1;
  }
  for (int mask = 1; mask < count; ++mask) {
    unsigned vertices = static_cast<unsigned>(mask);
    while (vertices) {
      const int v = __builtin_ctz(vertices);
      vertices &= vertices - 1;
      const int rest = mask ^ (1 << v);
      unsigned pred = static_cast<unsigned>(rest) & incoming[v];
      while (pred) {
        const int u = __builtin_ctz(pred);
        pred &= pred - 1;
        ending[mask][v] += ending[rest][u];
      }
      unsigned succ = static_cast<unsigned>(rest) & t.out[v];
      while (succ) {
        const int u = __builtin_ctz(succ);
        succ &= succ - 1;
        starting[mask][v] += starting[rest][u];
      }
    }
  }
  for (int v = 0; v < n; ++v) before[0][v] = after[0][v] = 1;
  for (int mask = 1; mask < count; ++mask) {
    for (int v = 0; v < n; ++v) {
      unsigned pred = static_cast<unsigned>(mask) & incoming[v];
      while (pred) {
        const int u = __builtin_ctz(pred);
        pred &= pred - 1;
        before[mask][v] += ending[mask][u];
      }
      unsigned succ = static_cast<unsigned>(mask) & t.out[v];
      while (succ) {
        const int u = __builtin_ctz(succ);
        succ &= succ - 1;
        after[mask][v] += starting[mask][u];
      }
    }
  }
  Exposure ans;
  const unsigned full = static_cast<unsigned>(count - 1);
  for (int v = 0; v < n; ++v) {
    ans.end[v] = ending[count - 1][v];
    ans.start[v] = starting[count - 1][v];
    ans.H += ans.end[v];
  }
  for (int x = 0; x < n; ++x) {
    for (int y = 0; y < n; ++y) if (x != y) {
      const unsigned rem = full ^ (1u << x) ^ (1u << y);
      unsigned left = rem;
      for (;;) {
        ans.Q[x][y] += before[left][x] * after[rem ^ left][y];
        if (left == 0) break;
        left = (left - 1) & rem;
      }
    }
  }
  for (int x = 0; x < n; ++x) for (int y = x + 1; y < n; ++y) {
    ans.c[x][y] = ans.c[y][x] = ans.Q[x][y] + ans.Q[y][x];
  }
  return ans;
}

static Packet packet(const Tournament& t, const Exposure& cap, int active_n = -1) {
  if (active_n < 0) active_n = t.n;
  std::array<i64, 12> d{}, h{};
  i64 W = 0, squares = 0;
  for (int i = 0; i < active_n; ++i) for (int j = i + 1; j < active_n; ++j) {
    const i64 value = cap.c[i][j];
    W += value;
    squares += value * value;
    d[i] += value;
    d[j] += value;
    if ((t.out[i] >> j) & 1u) { h[i] += value; h[j] -= value; }
    else { h[i] -= value; h[j] += value; }
  }
  Packet ans;
  i64 degree_squares = 0;
  for (int i = 0; i < active_n; ++i) {
    ans.C += d[i] * h[i];
    degree_squares += d[i] * d[i];
  }
  const i64 numerator = W * W + squares - degree_squares;
  if (numerator & 1) throw std::runtime_error("odd disjoint-energy numerator");
  ans.D = numerator / 2;
  return ans;
}

static i64 gate(const Packet& p, int eta) { return p.D + 2 * eta * p.C; }

static std::pair<i64, i64> root_split(const Tournament& q,
                                      const Tournament& child,
                                      const Exposure& cap) {
  const int n = q.n, x = n;
  Exposure old;
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) old.c[i][j] = cap.c[i][j];
  const Packet base = packet(q, old);
  std::array<i64, 12> d{}, h{}, a{};
  i64 W = 0;
  for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
    const i64 value = old.c[i][j];
    W += value; d[i] += value; d[j] += value;
    if ((q.out[i] >> j) & 1u) { h[i] += value; h[j] -= value; }
    else { h[i] -= value; h[j] += value; }
  }
  i64 D = base.D, C = base.C, A = 0, S = 0;
  for (int i = 0; i < n; ++i) {
    a[i] = cap.c[i][x];
    const int sigma = ((child.out[i] >> x) & 1u) ? 1 : -1;
    D += a[i] * (W - d[i]);
    C += a[i] * h[i] + sigma * a[i] * d[i] + sigma * a[i] * a[i];
    A += a[i]; S += sigma * a[i];
  }
  C -= A * S;
  return {C, D};
}

static std::vector<i64> all_hamilton_counts(const Tournament& t) {
  const int count = 1 << t.n;
  std::vector<std::array<i64, 12>> ending(count);
  std::vector<i64> H(count);
  for (int i = 0; i < t.n; ++i) ending[1 << i][i] = 1;
  for (int mask = 1; mask < count; ++mask) {
    for (int v = 0; v < t.n; ++v) if ((mask >> v) & 1) {
      const int rest = mask ^ (1 << v);
      for (int u = 0; u < t.n; ++u) {
        if (((rest >> u) & 1) && ((t.out[u] >> v) & 1u)) ending[mask][v] += ending[rest][u];
      }
      H[mask] += ending[mask][v];
    }
  }
  H[0] = 1;
  return H;
}

static std::vector<Exposure> odd_path_layers(const Tournament& t) {
  const int layer_count = t.n / 2;
  const int full = (1 << t.n) - 1;
  const auto H = all_hamilton_counts(t);
  std::vector<Exposure> layers(layer_count);
  for (int start = 0; start < t.n; ++start) {
    auto dfs = [&](auto&& self, int v, int mask, int length) -> void {
      if (length > 0 && (length & 1)) {
        const int layer = (length - 1) / 2;
        const int i = std::min(start, v), j = std::max(start, v);
        layers[layer].c[i][j] += 2 * H[full ^ mask];
        layers[layer].c[j][i] = layers[layer].c[i][j];
      }
      if (length == t.n - 1) return;
      unsigned next = t.out[v] & static_cast<unsigned>(~mask);
      while (next) {
        const int u = __builtin_ctz(next);
        next &= next - 1;
        self(self, u, mask | (1 << u), length + 1);
      }
    };
    dfs(dfs, start, 1 << start, 0);
  }
  return layers;
}

static Exposure tensor_sum(const Exposure& a, const Exposure& b, int n) {
  Exposure ans;
  for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
    ans.c[i][j] = ans.c[j][i] = a.c[i][j] + b.c[i][j];
  }
  return ans;
}

struct LayerStats {
  i64 rows = 0, path_edge_checks = 0, path_edge_failures = 0;
  i64 plus_self_checks = 0, plus_self_negative = 0, plus_self_zero = 0;
  i64 plus_cross_checks = 0, plus_cross_negative = 0, plus_cross_zero = 0;
  i64 minus_self_checks = 0, minus_self_negative = 0, minus_self_zero = 0;
  i64 minus_cross_checks = 0, minus_cross_negative = 0, minus_cross_zero = 0;
  i64 min_plus_self = std::numeric_limits<i64>::max();
  i64 min_plus_cross = std::numeric_limits<i64>::max();
  i64 min_minus_self = std::numeric_limits<i64>::max();
  i64 min_minus_cross = std::numeric_limits<i64>::max();
  std::vector<std::string> negative_plus_rows;
};

struct Census {
  int q = 0;
  i64 parents = 0, presentations = 0, strong_parents = 0, prime_parents = 0;
  i64 sink_parents = 0, source_parents = 0;
  i64 root_formula_checks = 0, root_formula_failures = 0;
  i64 root_split_checks = 0, root_split_failures = 0;
  i64 plus_negative = 0, plus_zero = 0, minus_negative = 0, minus_zero = 0;
  i64 child_sinks = 0, child_sources = 0;
  i64 sink_plus_negative = 0, sink_plus_zero = 0, sink_plus_positive = 0;
  i64 source_minus_negative = 0, source_minus_zero = 0, source_minus_positive = 0;
  i64 nosink_plus_nonpositive = 0, nosource_minus_nonpositive = 0;
  i64 strong_children = 0, prime_children = 0;
  i64 strong_child_gate_failures = 0, prime_child_gate_failures = 0;
  i64 strong_parent_nonuniform = 0, strong_parent_nonuniform_failures = 0;
  i64 prime_parent_nonuniform = 0, prime_parent_nonuniform_failures = 0;
  i64 prime_parent_prime_children = 0;
  i64 plus_parent_existence_disagree = 0, minus_parent_existence_disagree = 0;
  i64 plus_corner_min_disagree = 0, minus_corner_min_disagree = 0;
  i64 min_strong_nonuniform_gate = std::numeric_limits<i64>::max();
  i64 min_prime_nonuniform_gate = std::numeric_limits<i64>::max();
  LayerStats layers;
};

static void update_min(i64 value, i64& target) { if (value < target) target = value; }

static void audit_layers(const Tournament& q, const Exposure& cap, Census& s) {
  const auto layers = odd_path_layers(q);
  const auto subset_H = all_hamilton_counts(q);
  const int full = (1 << q.n) - 1;
  ++s.layers.rows;
  for (int i = 0; i < q.n; ++i) for (int j = i + 1; j < q.n; ++j) {
    i64 sum = 0;
    for (const auto& layer : layers) sum += layer.c[i][j];
    ++s.layers.path_edge_checks;
    if (sum != cap.c[i][j]) ++s.layers.path_edge_failures;
    if (!layers.empty()) {
      const i64 direct = 2 * subset_H[full ^ (1 << i) ^ (1 << j)];
      if (layers[0].c[i][j] != direct) ++s.layers.path_edge_failures;
    }
  }
  const bool no_sink = !has_sink(q), no_source = !has_source(q);
  std::vector<i64> plus_self(layers.size()), minus_self(layers.size());
  const i64 full_plus = gate(packet(q, cap), +1);
  for (std::size_t i = 0; i < layers.size(); ++i) {
    const Packet p = packet(q, layers[i]);
    plus_self[i] = gate(p, +1); minus_self[i] = gate(p, -1);
    if (no_sink) {
      ++s.layers.plus_self_checks;
      s.layers.plus_self_negative += plus_self[i] < 0;
      s.layers.plus_self_zero += plus_self[i] == 0;
      update_min(plus_self[i], s.layers.min_plus_self);
    }
    if (no_source) {
      ++s.layers.minus_self_checks;
      s.layers.minus_self_negative += minus_self[i] < 0;
      s.layers.minus_self_zero += minus_self[i] == 0;
      update_min(minus_self[i], s.layers.min_minus_self);
    }
  }
  std::vector<std::vector<i64>> plus_cross(layers.size(), std::vector<i64>(layers.size()));
  for (std::size_t i = 0; i < layers.size(); ++i) for (std::size_t j = i + 1; j < layers.size(); ++j) {
    const Packet p = packet(q, tensor_sum(layers[i], layers[j], q.n));
    plus_cross[i][j] = gate(p, +1) - plus_self[i] - plus_self[j];
    const i64 minus_cross = gate(p, -1) - minus_self[i] - minus_self[j];
    if (no_sink) {
      ++s.layers.plus_cross_checks;
      s.layers.plus_cross_negative += plus_cross[i][j] < 0;
      s.layers.plus_cross_zero += plus_cross[i][j] == 0;
      update_min(plus_cross[i][j], s.layers.min_plus_cross);
    }
    if (no_source) {
      ++s.layers.minus_cross_checks;
      s.layers.minus_cross_negative += minus_cross < 0;
      s.layers.minus_cross_zero += minus_cross == 0;
      update_min(minus_cross, s.layers.min_minus_cross);
    }
  }
  if (no_sink) {
    for (std::size_t i = 0; i < layers.size(); ++i) if (plus_self[i] < 0) {
      i64 other_self = 0, cross = 0;
      std::ostringstream detail;
      for (std::size_t j = 0; j < layers.size(); ++j) if (j != i) other_self += plus_self[j];
      for (std::size_t a = 0; a < layers.size(); ++a) for (std::size_t b = a + 1; b < layers.size(); ++b) {
        cross += plus_cross[a][b];
      }
      detail << "q=" << q.n << " label=" << q.label << " ell=" << (2 * i + 1)
             << " self=" << plus_self[i] << " full=" << full_plus
             << " other_self=" << other_self << " all_cross=" << cross
             << " strong=" << strong(q) << " prime=" << prime(q)
             << " source=" << has_source(q) << " min_module=" << minimum_module_size(q)
             << " incident_cross=";
      bool first = true;
      for (std::size_t j = 0; j < layers.size(); ++j) if (j != i) {
        if (!first) detail << ',';
        first = false;
        const i64 value = j < i ? plus_cross[j][i] : plus_cross[i][j];
        detail << (2 * j + 1) << ':' << value;
      }
      s.layers.negative_plus_rows.push_back(detail.str());
    }
  }
}

static void audit_parent(const Tournament& q, Census& s) {
  const Exposure parent = exposures(q);
  const bool parent_strong = strong(q), parent_prime = prime(q);
  const bool parent_sink = has_sink(q), parent_source = has_source(q);
  ++s.parents;
  s.strong_parents += parent_strong;
  s.prime_parents += parent_prime;
  s.sink_parents += parent_sink;
  s.source_parents += parent_source;
  audit_layers(q, parent, s);

  const int patterns = 1 << q.n, full = patterns - 1;
  std::vector<i64> plus(patterns), minus(patterns);
  for (int z = 0; z < patterns; ++z) {
    const Tournament child = extend(q, static_cast<unsigned>(z));
    const Exposure cap = exposures(child);
    const Packet p = packet(child, cap);
    plus[z] = gate(p, +1); minus[z] = gate(p, -1);
    ++s.presentations;
    s.plus_negative += plus[z] < 0; s.plus_zero += plus[z] == 0;
    s.minus_negative += minus[z] < 0; s.minus_zero += minus[z] == 0;
    const bool sink = has_sink(child), source = has_source(child);
    s.child_sinks += sink; s.child_sources += source;
    if (sink) {
      s.sink_plus_negative += plus[z] < 0;
      s.sink_plus_zero += plus[z] == 0;
      s.sink_plus_positive += plus[z] > 0;
    } else s.nosink_plus_nonpositive += plus[z] <= 0;
    if (source) {
      s.source_minus_negative += minus[z] < 0;
      s.source_minus_zero += minus[z] == 0;
      s.source_minus_positive += minus[z] > 0;
    } else s.nosource_minus_nonpositive += minus[z] <= 0;
    const bool child_strong = strong(child), child_prime = prime(child);
    s.strong_children += child_strong; s.prime_children += child_prime;
    s.strong_child_gate_failures += child_strong && (plus[z] <= 0 || minus[z] <= 0);
    s.prime_child_gate_failures += child_prime && (plus[z] <= 0 || minus[z] <= 0);
    if (parent_prime && child_prime) ++s.prime_parent_prime_children;
    if (z != 0 && z != full && parent_strong) {
      ++s.strong_parent_nonuniform;
      s.strong_parent_nonuniform_failures += plus[z] <= 0 || minus[z] <= 0;
      update_min(std::min(plus[z], minus[z]), s.min_strong_nonuniform_gate);
    }
    if (z != 0 && z != full && parent_prime) {
      ++s.prime_parent_nonuniform;
      s.prime_parent_nonuniform_failures += plus[z] <= 0 || minus[z] <= 0;
      update_min(std::min(plus[z], minus[z]), s.min_prime_nonuniform_gate);
    }

    const auto split = root_split(q, child, cap);
    ++s.root_split_checks;
    if (split.first != p.C || split.second != p.D) ++s.root_split_failures;
    for (int i = 0; i < q.n; ++i) {
      i64 predicted = parent.start[i] + parent.end[i];
      for (int j = 0; j < q.n; ++j) if (j != i) {
        predicted += ((z >> j) & 1) ? parent.Q[i][j] : parent.Q[j][i];
      }
      ++s.root_formula_checks;
      if (predicted != cap.c[i][q.n]) ++s.root_formula_failures;
    }
  }
  i64 plus_best = std::numeric_limits<i64>::max();
  i64 minus_best = std::numeric_limits<i64>::max();
  for (int z = 1; z < patterns; ++z) update_min(plus[z], plus_best);
  for (int z = 0; z < full; ++z) update_min(minus[z], minus_best);
  s.plus_parent_existence_disagree += ((plus_best <= 0) != parent_sink);
  s.minus_parent_existence_disagree += ((minus_best <= 0) != parent_source);
  s.plus_corner_min_disagree += ((plus_best < plus[0]) != parent_sink);
  s.minus_corner_min_disagree += ((minus_best < minus[full]) != parent_source);
}

static void print_optional_min(i64 value) {
  if (value == std::numeric_limits<i64>::max()) std::cout << "NA";
  else std::cout << value;
}

static void print_census(const Census& s) {
  std::cout << "q=" << s.q
            << " parents=" << s.parents << " presentations=" << s.presentations
            << " strong_parents=" << s.strong_parents << " prime_parents=" << s.prime_parents
            << " sink_parents=" << s.sink_parents << " source_parents=" << s.source_parents << '\n';
  std::cout << "  algebra root_formula=" << s.root_formula_checks << '/' << s.root_formula_failures
            << " root_split=" << s.root_split_checks << '/' << s.root_split_failures
            << " path_edges=" << s.layers.path_edge_checks << '/' << s.layers.path_edge_failures << '\n';
  std::cout << "  intrinsic plus_neg=" << s.plus_negative << " plus_zero=" << s.plus_zero
            << " sink_children=" << s.child_sinks
            << " sink_plus_neg_zero_pos=" << s.sink_plus_negative << ',' << s.sink_plus_zero << ',' << s.sink_plus_positive
            << " nosink_plus_nonpositive=" << s.nosink_plus_nonpositive
            << " minus_neg=" << s.minus_negative << " minus_zero=" << s.minus_zero
            << " source_children=" << s.child_sources
            << " source_minus_neg_zero_pos=" << s.source_minus_negative << ',' << s.source_minus_zero << ',' << s.source_minus_positive
            << " nosource_minus_nonpositive=" << s.nosource_minus_nonpositive << '\n';
  std::cout << "  child_controls strong=" << s.strong_children << " strong_gate_fail=" << s.strong_child_gate_failures
            << " prime=" << s.prime_children << " prime_gate_fail=" << s.prime_child_gate_failures << '\n';
  std::cout << "  parent_controls strong_nonuniform=" << s.strong_parent_nonuniform
            << " strong_fail=" << s.strong_parent_nonuniform_failures << " strong_min=";
  print_optional_min(s.min_strong_nonuniform_gate);
  std::cout << " prime_nonuniform=" << s.prime_parent_nonuniform
            << " prime_fail=" << s.prime_parent_nonuniform_failures << " prime_min=";
  print_optional_min(s.min_prime_nonuniform_gate);
  std::cout << " prime_parent_prime_children=" << s.prime_parent_prime_children << '\n';
  std::cout << "  parent_equivalences plus_nonzero_nonpositive_iff_sink_disagree=" << s.plus_parent_existence_disagree
            << " minus_nonfull_nonpositive_iff_source_disagree=" << s.minus_parent_existence_disagree
            << " plus_corner_min_iff_sink_disagree=" << s.plus_corner_min_disagree
            << " minus_corner_min_iff_source_disagree=" << s.minus_corner_min_disagree << '\n';
  std::cout << "  layers plus_self=" << s.layers.plus_self_checks << ',' << s.layers.plus_self_negative << ',' << s.layers.plus_self_zero << ',';
  print_optional_min(s.layers.min_plus_self);
  std::cout << " plus_cross=" << s.layers.plus_cross_checks << ',' << s.layers.plus_cross_negative << ',' << s.layers.plus_cross_zero << ',';
  print_optional_min(s.layers.min_plus_cross);
  std::cout << " minus_self=" << s.layers.minus_self_checks << ',' << s.layers.minus_self_negative << ',' << s.layers.minus_self_zero << ',';
  print_optional_min(s.layers.min_minus_self);
  std::cout << " minus_cross=" << s.layers.minus_cross_checks << ',' << s.layers.minus_cross_negative << ',' << s.layers.minus_cross_zero << ',';
  print_optional_min(s.layers.min_minus_cross);
  std::cout << '\n';
  for (const auto& row : s.layers.negative_plus_rows) std::cout << "  negative_plus_layer " << row << '\n';
}

int main() {
  std::array<Census, 9> census;
  for (int q = 0; q <= 8; ++q) census[q].q = q;
  std::string bits;
  while (std::cin >> bits) {
    const Tournament q = parse_bits(bits);
    if (q.n < 3 || q.n > 8) throw std::runtime_error("expected parent orders 3 through 8");
    audit_parent(q, census[q.n]);
  }
  std::cout << "THM4177_PRIMARY_V1\n";
  for (int q = 3; q <= 8; ++q) print_census(census[q]);
  return 0;
}
