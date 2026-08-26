#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <omp.h>
#include <openssl/sha.h>

// Independent complete audit for THM-4137.
//
// Input is the ordinary upper-triangle ASCII stream from
//
//     gentourng -q 10
//
// (all unlabeled order-ten tournaments, not the generator's strong-only
// stream).  This program performs its own strong-connectivity filter.  For
// every surviving record it computes the complete 1024-entry response table
// by a contracted adjacent-block engine.  The two extremal records are then
// replayed by a second, literal 11-vertex Hamilton-path DP.  No theorem-4137
// primary artifact is read or linked.

namespace {

constexpr int N = 10;
constexpr int SZ = 1 << N;
constexpr int FULL = SZ - 1;
constexpr uint16_t CENTRAL_MASK = uint16_t(1u << (N / 2));
constexpr uint64_t EXPECTED_ALL = 9733056;
constexpr uint64_t EXPECTED_STRONG = 9355949;
constexpr const char* EXPECTED_ALL_RAW_SHA256 =
    "af2873154068897522bc15477d989b0577877d2bbebc08aea3082353e5378b67";
constexpr const char* EXPECTED_STRONG_RAW_SHA256 =
    "47bcaa3ef6272261dee3092735b47b3d2154d882aae6a8420c964cd3ef7289b7";
constexpr const char* EXPECTED_PROFILE_SHA256 =
    "8b2d93cd83ccafa86315b540cbeba5338812ad4d6393712bd403163b954ebab9";
constexpr const char* EXPECTED_SEMANTIC_SHA256 =
    "2984d66aad38b7b7a5fa4de0ec14e6e08ff383bf277eaa2ae4290dca1841f045";

struct Rat {
  int64_t n = 0;
  int64_t d = 1;
};

int64_t gcd64(int64_t a, int64_t b) {
  if (a < 0) a = -a;
  if (b < 0) b = -b;
  return std::gcd(a, b);
}

Rat make_rat(int64_t n, int64_t d) {
  if (!d) throw std::runtime_error("zero rational denominator");
  if (d < 0) {
    n = -n;
    d = -d;
  }
  const int64_t g = gcd64(n, d);
  return {n / g, d / g};
}

int compare(Rat a, Rat b) {
  const __int128 left = static_cast<__int128>(a.n) * b.d;
  const __int128 right = static_cast<__int128>(b.n) * a.d;
  return (left > right) - (left < right);
}

std::string rational_string(Rat a) {
  return std::to_string(a.n) + "/" + std::to_string(a.d);
}

int64_t ceil_quotient(__int128 numerator, __int128 denominator) {
  if (denominator <= 0) throw std::runtime_error("bad ceil denominator");
  __int128 q = numerator / denominator;
  const __int128 r = numerator % denominator;
  if (r > 0) ++q;
  if (q < std::numeric_limits<int64_t>::min() ||
      q > std::numeric_limits<int64_t>::max())
    throw std::runtime_error("ceil quotient overflow");
  return static_cast<int64_t>(q);
}

struct Adjacency {
  std::array<uint16_t, N> out{};
  std::array<uint16_t, N> in{};
};

Adjacency decode(const std::string& label) {
  if (label.size() != N * (N - 1) / 2)
    throw std::runtime_error("input record is not a 45-bit tournament label");
  Adjacency a;
  size_t cursor = 0;
  for (int u = 0; u < N; ++u) {
    for (int v = u + 1; v < N; ++v) {
      const char bit = label[cursor++];
      if (bit == '1') a.out[u] |= uint16_t(1u << v);
      else if (bit == '0') a.out[v] |= uint16_t(1u << u);
      else throw std::runtime_error("nonbinary tournament label");
    }
  }
  for (int u = 0; u < N; ++u) {
    for (int v = 0; v < N; ++v) {
      if (a.out[u] & (1u << v)) a.in[v] |= uint16_t(1u << u);
    }
  }
  return a;
}

bool locally_strong(const Adjacency& a) {
  auto all_reached = [&](bool reverse) {
    uint16_t reached = 1;
    uint16_t frontier = 1;
    while (frontier) {
      const int v = __builtin_ctz(frontier);
      frontier &= uint16_t(frontier - 1);
      const uint16_t next = reverse ? a.in[v] : a.out[v];
      const uint16_t fresh = next & uint16_t(FULL ^ reached);
      reached |= fresh;
      frontier |= fresh;
    }
    return reached == FULL;
  };
  // Every vertex is reachable from 0, and every vertex reaches 0.
  return all_reached(false) && all_reached(true);
}

struct Layer {
  Rat J;
  int64_t floor = 0;
  int64_t maximum = 0;
  int64_t lattice = 0;
};

struct Profile {
  std::string label;
  int64_t H = 0;
  int64_t W2 = 0;       // 2 W
  int64_t D4x4 = 0;     // 4 D_4
  int64_t Chdx4 = 0;    // 4 C_hd
  Rat theta;
  Rat rho;
  uint16_t rational_mask = 0;
  uint16_t coset_mask = 0;
  uint16_t actual_mask = 0;
  std::array<Layer, N - 1> layers{};
  std::array<uint32_t, SZ> values{};
};

Profile contracted_profile(const std::string& label) {
  const Adjacency a = decode(label);
  if (!locally_strong(a))
    throw std::runtime_error("contracted evaluator received nonstrong input");

  Profile p;
  p.label = label;

  // Hamilton paths of each induced subset, split by their final/initial
  // vertex.  This recurrence is independent of the primary Python route.
  uint32_t ending[SZ][N] = {};
  uint32_t starting[SZ][N] = {};
  for (int v = 0; v < N; ++v) {
    ending[1 << v][v] = 1;
    starting[1 << v][v] = 1;
  }
  for (int mask = 1; mask < SZ; ++mask) {
    uint16_t vertices = uint16_t(mask);
    while (vertices) {
      const int v = __builtin_ctz(vertices);
      vertices &= uint16_t(vertices - 1);
      const int rest = mask ^ (1 << v);
      uint16_t possible = uint16_t(rest) & a.in[v];
      while (possible) {
        const int u = __builtin_ctz(possible);
        possible &= uint16_t(possible - 1);
        ending[mask][v] += ending[rest][u];
      }
      possible = uint16_t(rest) & a.out[v];
      while (possible) {
        const int u = __builtin_ctz(possible);
        possible &= uint16_t(possible - 1);
        starting[mask][v] += starting[rest][u];
      }
    }
  }
  for (int v = 0; v < N; ++v) p.H += ending[FULL][v];
  if (p.H <= 0) throw std::runtime_error("nonpositive Hamilton count");

  // before[L][v]: Hamilton orders of L whose last vertex points to v.
  // after[R][v]:  Hamilton orders of R whose first vertex is pointed to by v.
  uint32_t before[SZ][N] = {};
  uint32_t after[SZ][N] = {};
  for (int v = 0; v < N; ++v) before[0][v] = after[0][v] = 1;
  for (int mask = 1; mask < SZ; ++mask) {
    for (int boundary = 0; boundary < N; ++boundary) {
      uint16_t choices = uint16_t(mask) & a.in[boundary];
      while (choices) {
        const int u = __builtin_ctz(choices);
        choices &= uint16_t(choices - 1);
        before[mask][boundary] += ending[mask][u];
      }
      choices = uint16_t(mask) & a.out[boundary];
      while (choices) {
        const int u = __builtin_ctz(choices);
        choices &= uint16_t(choices - 1);
        after[mask][boundary] += starting[mask][u];
      }
    }
  }

  uint32_t edge_capacity[N][N] = {};
  int64_t outgoing_capacity[N] = {};
  struct UndirectedEdge { int u; int v; uint32_t capacity; };
  std::vector<UndirectedEdge> edges;
  edges.reserve(N * (N - 1) / 2);
  for (int x = 0; x < N; ++x) {
    for (int y = x + 1; y < N; ++y) {
      int u = x, v = y;
      if (!(a.out[u] & (1u << v))) std::swap(u, v);
      const int remaining = FULL ^ (1 << u) ^ (1 << v);
      uint64_t capacity = 0;
      for (int left = remaining;; left = (left - 1) & remaining) {
        const int right = remaining ^ left;
        // [u,v] is the good adjacent block; [v,u] is its reversed block.
        capacity += uint64_t(before[left][u]) * after[right][v];
        capacity += uint64_t(before[left][v]) * after[right][u];
        if (!left) break;
      }
      if (capacity > std::numeric_limits<uint32_t>::max())
        throw std::runtime_error("edge capacity overflow");
      const uint32_t c = static_cast<uint32_t>(capacity);
      edge_capacity[u][v] = edge_capacity[v][u] = c;
      outgoing_capacity[u] += c;
      p.W2 += c;
      edges.push_back({u, v, c});
    }
  }

  // Compute the exact scaled packet by a deliberately direct edge-pair loop.
  int64_t d2[N] = {};
  int64_t h2[N] = {};
  for (int v = 0; v < N; ++v) {
    for (int u = 0; u < N; ++u) if (u != v) {
      const int64_t c = edge_capacity[v][u];
      d2[v] += c;
      h2[v] += (a.out[v] & (1u << u)) ? c : -c;
    }
    p.Chdx4 += d2[v] * h2[v];
  }
  for (size_t i = 0; i < edges.size(); ++i) {
    for (size_t j = i + 1; j < edges.size(); ++j) {
      const auto& e = edges[i];
      const auto& f = edges[j];
      if (e.u != f.u && e.u != f.v && e.v != f.u && e.v != f.v)
        p.D4x4 += int64_t(e.capacity) * f.capacity;
    }
  }
  if (p.D4x4 <= 0) throw std::runtime_error("nonpositive disjoint-edge term");
  p.theta = make_rat((N - 3) * p.Chdx4, 2 * p.D4x4);
  p.rho = make_rat((N - 3) * std::llabs(p.Chdx4), 2 * p.D4x4);

  // The response is H plus the capacity of the directed cut S -> S^c.
  int64_t incident_to_subset[N][SZ] = {};
  for (int v = 0; v < N; ++v) {
    for (int mask = 1; mask < SZ; ++mask) {
      const int u = __builtin_ctz(mask);
      incident_to_subset[v][mask] =
          incident_to_subset[v][mask ^ (1 << u)] + edge_capacity[v][u];
    }
  }
  int64_t cut[SZ] = {};
  p.values[0] = static_cast<uint32_t>(p.H);
  for (int mask = 1; mask < SZ; ++mask) {
    const int v = __builtin_ctz(mask);
    const int rest = mask ^ (1 << v);
    cut[mask] = cut[rest] + outgoing_capacity[v]
              - incident_to_subset[v][rest];
    const int64_t value = p.H + cut[mask];
    if (value < p.H || value > std::numeric_limits<uint32_t>::max())
      throw std::runtime_error("response outside uint32 support");
    if (!(value & 1)) throw std::runtime_error("Redei parity failure");
    p.values[mask] = static_cast<uint32_t>(value);
  }
  if (p.values[FULL] != p.H)
    throw std::runtime_error("constant-ear endpoint mismatch");

  Rat best_J{std::numeric_limits<int64_t>::min(), 1};
  int64_t best_floor = std::numeric_limits<int64_t>::min();
  int64_t best_actual = std::numeric_limits<int64_t>::min();
  for (int m = 1; m < N; ++m) {
    int64_t count = 0;
    int64_t sum = 0;
    int64_t sum_of_squares = 0;
    int64_t anchor = -1;
    int64_t lattice = 0;
    int64_t maximum = 0;
    for (int mask = 0; mask < SZ; ++mask) {
      if (__builtin_popcount(mask) != m) continue;
      const int64_t value = p.values[mask];
      if (anchor < 0) anchor = value;
      ++count;
      sum += value;
      sum_of_squares += value * value;
      lattice = gcd64(lattice, value - anchor);
      maximum = std::max(maximum, value);
    }
    const int64_t mean_gap_numerator = sum - count * p.H;
    if (mean_gap_numerator <= 0)
      throw std::runtime_error("nonpositive layer mean gap");
    const Rat J = make_rat(sum_of_squares - sum * p.H,
                           mean_gap_numerator);
    int64_t support_floor = anchor;
    if (lattice) {
      const __int128 offset_numerator =
          static_cast<__int128>(J.n) - static_cast<__int128>(anchor) * J.d;
      support_floor += lattice * ceil_quotient(
          offset_numerator, static_cast<__int128>(J.d) * lattice);
    } else if (compare(J, Rat{anchor, 1}) != 0) {
      throw std::runtime_error("constant layer has nonconstant ratio");
    }
    if (support_floor > maximum)
      throw std::runtime_error("coset floor above actual support");
    p.layers[m - 1] = {J, support_floor, maximum, lattice};

    const int jc = compare(J, best_J);
    if (jc > 0) {
      best_J = J;
      p.rational_mask = uint16_t(1u << m);
    } else if (!jc) {
      p.rational_mask |= uint16_t(1u << m);
    }
    if (support_floor > best_floor) {
      best_floor = support_floor;
      p.coset_mask = uint16_t(1u << m);
    } else if (support_floor == best_floor) {
      p.coset_mask |= uint16_t(1u << m);
    }
    if (maximum > best_actual) {
      best_actual = maximum;
      p.actual_mask = uint16_t(1u << m);
    } else if (maximum == best_actual) {
      p.actual_mask |= uint16_t(1u << m);
    }
  }

  // A second exact expression checks the rational optimizer computation:
  // the minimizing integer t=n-2m must be nearest to theta.
  Rat best_distance{std::numeric_limits<int64_t>::max(), 1};
  uint16_t nearest_mask = 0;
  for (int m = 1; m < N; ++m) {
    const int t = N - 2 * m;
    const Rat distance = make_rat(
        std::llabs(int64_t(t) * p.theta.d - p.theta.n), p.theta.d);
    const int c = compare(distance, best_distance);
    if (c < 0) {
      best_distance = distance;
      nearest_mask = uint16_t(1u << m);
    } else if (!c) {
      nearest_mask |= uint16_t(1u << m);
    }
  }
  if (nearest_mask != p.rational_mask)
    throw std::runtime_error("rational layer / tilt nearest-grid mismatch");
  return p;
}

std::string t_tuple(uint16_t mask, int order = N) {
  std::ostringstream out;
  out << '(';
  bool first = true;
  for (int m = 1; m < order; ++m) {
    if (!(mask & (1u << m))) continue;
    if (!first) out << ',';
    first = false;
    out << order - 2 * m;
  }
  out << ')';
  return out.str();
}

void sha_update(SHA256_CTX& ctx, const void* data, size_t size) {
  SHA256_Update(&ctx, data, size);
}

std::string sha_finish(SHA256_CTX& ctx) {
  unsigned char bytes[SHA256_DIGEST_LENGTH];
  SHA256_Final(bytes, &ctx);
  std::ostringstream out;
  for (unsigned char byte : bytes)
    out << std::hex << std::setw(2) << std::setfill('0') << int(byte);
  return out.str();
}

std::string sha_string(const std::string& text) {
  SHA256_CTX ctx;
  SHA256_Init(&ctx);
  sha_update(ctx, text.data(), text.size());
  return sha_finish(ctx);
}

template <typename U>
void append_unsigned_le(std::vector<unsigned char>& bytes, U value) {
  static_assert(std::is_unsigned<U>::value, "unsigned serializer only");
  for (size_t i = 0; i < sizeof(U); ++i) {
    bytes.push_back(static_cast<unsigned char>(value & 0xffu));
    value >>= 8;
  }
}

void append_i64(std::vector<unsigned char>& bytes, int64_t value) {
  append_unsigned_le<uint64_t>(bytes, static_cast<uint64_t>(value));
}

void hash_profile(SHA256_CTX& ctx, const Profile& p) {
  std::vector<unsigned char> bytes;
  bytes.reserve(4300);
  bytes.insert(bytes.end(), p.label.begin(), p.label.end());
  append_unsigned_le<uint32_t>(bytes, static_cast<uint32_t>(p.H));
  append_i64(bytes, p.W2);
  append_i64(bytes, p.D4x4);
  append_i64(bytes, p.Chdx4);
  append_i64(bytes, p.theta.n);
  append_i64(bytes, p.theta.d);
  append_i64(bytes, p.rho.n);
  append_i64(bytes, p.rho.d);
  append_unsigned_le<uint16_t>(bytes, p.rational_mask);
  append_unsigned_le<uint16_t>(bytes, p.coset_mask);
  append_unsigned_le<uint16_t>(bytes, p.actual_mask);
  for (const Layer& layer : p.layers) {
    append_i64(bytes, layer.J.n);
    append_i64(bytes, layer.J.d);
    append_i64(bytes, layer.floor);
    append_i64(bytes, layer.maximum);
    append_i64(bytes, layer.lattice);
  }
  for (uint32_t value : p.values) append_unsigned_le<uint32_t>(bytes, value);
  sha_update(ctx, bytes.data(), bytes.size());
}

std::string profile_packet(const Profile& p, bool include_layers = true) {
  std::ostringstream out;
  out << "label=" << p.label
      << ",H=" << p.H
      << ",W=" << rational_string(make_rat(p.W2, 2))
      << ",D4=" << rational_string(make_rat(p.D4x4, 4))
      << ",Chd=" << rational_string(make_rat(p.Chdx4, 4))
      << ",theta=" << rational_string(p.theta)
      << ",rho=" << rational_string(p.rho)
      << ",rational_t=" << t_tuple(p.rational_mask)
      << ",coset_t=" << t_tuple(p.coset_mask)
      << ",actual_t=" << t_tuple(p.actual_mask);
  if (include_layers) {
    out << ",layers=[";
    for (int m = 1; m < N; ++m) {
      if (m > 1) out << ';';
      const Layer& layer = p.layers[m - 1];
      out << '(' << m << ',' << (N - 2 * m) << ','
          << rational_string(layer.J) << ',' << layer.floor << ','
          << layer.maximum << ',' << layer.lattice << ')';
    }
    out << ']';
  }
  return out.str();
}

// Literal second evaluator: build each T+x_S and count its Hamilton paths
// without using capacities, packets, or the cut recurrence.
std::array<uint32_t, SZ> literal_child_values(const std::string& label) {
  constexpr int M = N + 1;
  constexpr int MSZ = 1 << M;
  const Adjacency parent = decode(label);
  std::array<uint32_t, SZ> answer{};
  for (int cut = 0; cut < SZ; ++cut) {
    uint16_t child_out[M] = {};
    for (int v = 0; v < N; ++v) child_out[v] = parent.out[v];
    for (int v = 0; v < N; ++v) {
      if (cut & (1 << v)) child_out[N] |= uint16_t(1u << v);
      else child_out[v] |= uint16_t(1u << N);
    }
    uint32_t ending[MSZ][M] = {};
    for (int v = 0; v < M; ++v) ending[1 << v][v] = 1;
    for (int mask = 1; mask < MSZ; ++mask) {
      for (int v = 0; v < M; ++v) if (mask & (1 << v)) {
        const int rest = mask ^ (1 << v);
        for (int u = 0; u < M; ++u) {
          if ((rest & (1 << u)) && (child_out[u] & (1 << v)))
            ending[mask][v] += ending[rest][u];
        }
      }
    }
    for (int v = 0; v < M; ++v) answer[cut] += ending[MSZ - 1][v];
  }
  return answer;
}

// Find an explicit isomorphism reverse(A) -> B.  Degree filtering makes the
// ten-vertex backtracking negligible and avoids relying on nauty a second time.
bool reversal_isomorphism(const std::string& left_label,
                          const std::string& right_label,
                          std::array<int, N>* mapping_out) {
  const Adjacency left = decode(left_label);
  const Adjacency right = decode(right_label);
  std::array<std::vector<int>, N> candidates;
  std::array<int, N> order{};
  for (int v = 0; v < N; ++v) {
    order[v] = v;
    const int reversed_degree = N - 1 - __builtin_popcount(left.out[v]);
    for (int w = 0; w < N; ++w) {
      if (__builtin_popcount(right.out[w]) == reversed_degree)
        candidates[v].push_back(w);
    }
  }
  std::stable_sort(order.begin(), order.end(), [&](int u, int v) {
    return candidates[u].size() < candidates[v].size();
  });
  std::array<int, N> mapping;
  mapping.fill(-1);
  uint16_t used = 0;
  auto search = [&](auto&& self, int depth) -> bool {
    if (depth == N) return true;
    const int v = order[depth];
    for (int w : candidates[v]) {
      if (used & (1u << w)) continue;
      bool compatible = true;
      for (int k = 0; k < depth; ++k) {
        const int u = order[k];
        const int z = mapping[u];
        // reverse(left) has v -> u exactly when left has u -> v.
        const bool source_arc = left.out[u] & (1u << v);
        const bool target_arc = right.out[w] & (1u << z);
        if (source_arc != target_arc) {
          compatible = false;
          break;
        }
      }
      if (!compatible) continue;
      mapping[v] = w;
      used |= uint16_t(1u << w);
      if (self(self, depth + 1)) return true;
      used &= uint16_t(~(1u << w));
      mapping[v] = -1;
    }
    return false;
  };
  if (!search(search, 0)) return false;
  *mapping_out = mapping;
  return true;
}

// Literal n=6 hostile controls.  These deliberately do not invoke the main
// contracted engine and keep rational, exact-coset, and actual maxima separate.
std::string hostile_control(uint64_t code) {
  constexpr int P = 6;
  constexpr int CHILD = P + 1;
  constexpr int PSZ = 1 << P;
  constexpr int CSZ = 1 << CHILD;
  std::string label;
  for (int cursor = 0; cursor < P * (P - 1) / 2; ++cursor)
    label.push_back((code & (uint64_t(1) << cursor)) ? '1' : '0');
  uint16_t parent_out[P] = {};
  int cursor = 0;
  for (int u = 0; u < P; ++u) for (int v = u + 1; v < P; ++v) {
    if (label[cursor++] == '1') parent_out[u] |= uint16_t(1u << v);
    else parent_out[v] |= uint16_t(1u << u);
  }
  uint32_t values[PSZ] = {};
  for (int cut = 0; cut < PSZ; ++cut) {
    uint16_t child_out[CHILD] = {};
    for (int v = 0; v < P; ++v) child_out[v] = parent_out[v];
    for (int v = 0; v < P; ++v) {
      if (cut & (1 << v)) child_out[P] |= uint16_t(1u << v);
      else child_out[v] |= uint16_t(1u << P);
    }
    uint32_t ending[CSZ][CHILD] = {};
    for (int v = 0; v < CHILD; ++v) ending[1 << v][v] = 1;
    for (int mask = 1; mask < CSZ; ++mask) {
      for (int v = 0; v < CHILD; ++v) if (mask & (1 << v)) {
        const int rest = mask ^ (1 << v);
        for (int u = 0; u < CHILD; ++u)
          if ((rest & (1 << u)) && (child_out[u] & (1 << v)))
            ending[mask][v] += ending[rest][u];
      }
    }
    for (int v = 0; v < CHILD; ++v) values[cut] += ending[CSZ - 1][v];
  }
  const int64_t H = values[0];
  Rat best_J{std::numeric_limits<int64_t>::min(), 1};
  int64_t best_floor = std::numeric_limits<int64_t>::min();
  int64_t best_actual = std::numeric_limits<int64_t>::min();
  uint16_t rational_mask = 0, coset_mask = 0, actual_mask = 0;
  struct SmallLayer { Rat J; int64_t floor; int64_t maximum; int64_t lattice; };
  SmallLayer layers[P - 1] = {};
  for (int m = 1; m < P; ++m) {
    int64_t count = 0, sum = 0, sumsq = 0, anchor = -1, lattice = 0, maximum = 0;
    for (int mask = 0; mask < PSZ; ++mask) if (__builtin_popcount(mask) == m) {
      const int64_t value = values[mask];
      if (anchor < 0) anchor = value;
      ++count;
      sum += value;
      sumsq += value * value;
      lattice = gcd64(lattice, value - anchor);
      maximum = std::max(maximum, value);
    }
    const Rat J = make_rat(sumsq - sum * H, sum - count * H);
    int64_t floor = anchor;
    if (lattice) {
      floor += lattice * ceil_quotient(
          static_cast<__int128>(J.n) - static_cast<__int128>(anchor) * J.d,
          static_cast<__int128>(J.d) * lattice);
    }
    layers[m - 1] = {J, floor, maximum, lattice};
    const int c = compare(J, best_J);
    if (c > 0) { best_J = J; rational_mask = uint16_t(1u << m); }
    else if (!c) rational_mask |= uint16_t(1u << m);
    if (floor > best_floor) { best_floor = floor; coset_mask = uint16_t(1u << m); }
    else if (floor == best_floor) coset_mask |= uint16_t(1u << m);
    if (maximum > best_actual) { best_actual = maximum; actual_mask = uint16_t(1u << m); }
    else if (maximum == best_actual) actual_mask |= uint16_t(1u << m);
  }
  std::ostringstream out;
  out << label << "|H=" << H
      << "|rational_t=" << t_tuple(rational_mask, P)
      << "|coset_t=" << t_tuple(coset_mask, P)
      << "|actual_t=" << t_tuple(actual_mask, P)
      << "|layers=[";
  for (int m = 1; m < P; ++m) {
    if (m > 1) out << ';';
    const auto& layer = layers[m - 1];
    out << '(' << m << ',' << (P - 2 * m) << ','
        << rational_string(layer.J) << ',' << layer.floor << ','
        << layer.maximum << ',' << layer.lattice << ')';
  }
  out << "]|literal_child_dp=PASS";
  return out.str();
}

std::string histogram_string(const std::map<std::string, uint64_t>& histogram) {
  std::ostringstream out;
  bool first = true;
  for (const auto& item : histogram) {
    if (!first) out << ',';
    first = false;
    out << item.first << ':' << item.second;
  }
  return out.str();
}

}  // namespace

int main() {
  try {
    uint64_t input_limit = std::numeric_limits<uint64_t>::max();
    if (const char* value = std::getenv("AUDIT_LIMIT"))
      input_limit = std::stoull(value);

    SHA256_CTX all_raw_ctx, strong_raw_ctx, profile_ctx;
    SHA256_Init(&all_raw_ctx);
    SHA256_Init(&strong_raw_ctx);
    SHA256_Init(&profile_ctx);
    constexpr char PROFILE_HEADER[] = "thm4137-n10-independent-profile-v1\n";
    sha_update(profile_ctx, PROFILE_HEADER, sizeof(PROFILE_HEADER) - 1);

    uint64_t all_count = 0, strong_count = 0;
    uint64_t rational_failures = 0, coset_failures = 0;
    uint64_t theta_zero = 0, coset_reorders = 0, actual_noncentral_only = 0;
    int64_t minimum_margin = std::numeric_limits<int64_t>::max();
    Rat worst_rho{-1, 1};
    std::vector<Profile> worst_profiles;
    std::vector<Profile> minimum_margin_profiles;
    std::map<std::string, uint64_t> rational_histogram;
    std::map<std::string, uint64_t> coset_histogram;
    std::map<std::string, uint64_t> actual_histogram;

    constexpr size_t BATCH = 512;
    std::vector<std::string> batch_labels;
    batch_labels.reserve(BATCH);

    auto process_batch = [&]() {
      if (batch_labels.empty()) return;
      std::vector<Profile> profiles(batch_labels.size());
      #pragma omp parallel for schedule(static)
      for (int64_t i = 0; i < static_cast<int64_t>(batch_labels.size()); ++i)
        profiles[i] = contracted_profile(batch_labels[size_t(i)]);
      for (Profile& p : profiles) {
        ++strong_count;
        hash_profile(profile_ctx, p);
        if (p.rational_mask & ~CENTRAL_MASK) ++rational_failures;
        if (p.coset_mask & ~CENTRAL_MASK) ++coset_failures;
        if (p.rational_mask != p.coset_mask) ++coset_reorders;
        if (!(p.actual_mask & CENTRAL_MASK)) ++actual_noncentral_only;
        if (!p.theta.n) ++theta_zero;
        ++rational_histogram[t_tuple(p.rational_mask)];
        ++coset_histogram[t_tuple(p.coset_mask)];
        ++actual_histogram[t_tuple(p.actual_mask)];

        int64_t central_floor = p.layers[N / 2 - 1].floor;
        int64_t outer_floor = std::numeric_limits<int64_t>::min();
        for (int m = 1; m < N; ++m) {
          if (m != N / 2) outer_floor = std::max(outer_floor, p.layers[m - 1].floor);
        }
        const int64_t margin = central_floor - outer_floor;
        if (margin < minimum_margin) {
          minimum_margin = margin;
          minimum_margin_profiles.assign(1, p);
        } else if (margin == minimum_margin) {
          minimum_margin_profiles.push_back(p);
        }
        const int wc = compare(p.rho, worst_rho);
        if (wc > 0) {
          worst_rho = p.rho;
          worst_profiles.assign(1, p);
        } else if (!wc) {
          worst_profiles.push_back(p);
        }
      }
      batch_labels.clear();
    };

    std::string label;
    const char newline = '\n';
    while (all_count < input_limit && (std::cin >> label)) {
      ++all_count;
      sha_update(all_raw_ctx, label.data(), label.size());
      sha_update(all_raw_ctx, &newline, 1);
      const Adjacency a = decode(label);
      if (locally_strong(a)) {
        sha_update(strong_raw_ctx, label.data(), label.size());
        sha_update(strong_raw_ctx, &newline, 1);
        batch_labels.push_back(label);
        if (batch_labels.size() == BATCH) process_batch();
      }
      if (all_count % 131072 == 0)
        std::cerr << "input_records=" << all_count << ",strong_processed="
                  << strong_count << '\n';
    }
    process_batch();
    if (!all_count || !strong_count) throw std::runtime_error("empty audit universe");

    const std::string all_raw_hash = sha_finish(all_raw_ctx);
    const std::string strong_raw_hash = sha_finish(strong_raw_ctx);
    const std::string profile_hash = sha_finish(profile_ctx);

    std::string worst_literal = "NOT_RUN_PARTIAL";
    std::string reversal_pair = "NOT_RUN_PARTIAL";
    std::string reversal_mapping = "()";
    if (all_count == EXPECTED_ALL) {
      if (worst_profiles.size() != 2)
        throw std::runtime_error("full universe did not have two worst profiles");
      for (const Profile& p : worst_profiles) {
        if (literal_child_values(p.label) != p.values)
          throw std::runtime_error("literal child replay disagrees at worst profile");
      }
      worst_literal = "PASS_2_OF_2_ALL_1024_CUTS";
      std::array<int, N> mapping{};
      if (!reversal_isomorphism(worst_profiles[0].label, worst_profiles[1].label,
                                &mapping))
        throw std::runtime_error("worst profiles are not reversal-isomorphic");
      reversal_pair = "PASS";
      std::ostringstream map_out;
      map_out << '(';
      for (int v = 0; v < N; ++v) {
        if (v) map_out << ',';
        map_out << v << "->" << mapping[v];
      }
      map_out << ')';
      reversal_mapping = map_out.str();
    }

    const std::string control_2 = hostile_control(2);
    const std::string control_140 = hostile_control(140);
    const std::string control_20 = hostile_control(20);

    std::vector<std::string> facts;
    auto fact = [&](const std::string& key, const std::string& value) {
      facts.push_back(key + "=" + value);
    };
    fact("all_unlabeled_classes", std::to_string(all_count));
    fact("locally_strong_classes", std::to_string(strong_count));
    fact("locally_nonstrong_classes", std::to_string(all_count - strong_count));
    fact("all_raw_stream_sha256", all_raw_hash);
    fact("locally_filtered_strong_stream_sha256", strong_raw_hash);
    fact("ordered_full_profile_sha256", profile_hash);
    fact("rational_central_failures", std::to_string(rational_failures));
    fact("coset_central_failures", std::to_string(coset_failures));
    fact("coset_reorders_rational", std::to_string(coset_reorders));
    fact("minimum_strict_coset_margin", std::to_string(minimum_margin));
    fact("minimum_margin_multiplicity", std::to_string(minimum_margin_profiles.size()));
    for (size_t i = 0; i < minimum_margin_profiles.size(); ++i)
      fact("minimum_margin_packet_" + std::to_string(i + 1),
           profile_packet(minimum_margin_profiles[i], false));
    fact("theta_zero", std::to_string(theta_zero));
    fact("worst_rho", rational_string(worst_rho));
    fact("worst_multiplicity", std::to_string(worst_profiles.size()));
    for (size_t i = 0; i < worst_profiles.size(); ++i)
      fact("worst_packet_" + std::to_string(i + 1), profile_packet(worst_profiles[i]));
    fact("worst_literal_child_replay", worst_literal);
    fact("worst_reversal_isomorphism", reversal_pair);
    fact("worst_reversal_mapping", reversal_mapping);
    fact("rational_optimizer_histogram", histogram_string(rational_histogram));
    fact("coset_optimizer_histogram", histogram_string(coset_histogram));
    fact("actual_noncentral_only", std::to_string(actual_noncentral_only));
    fact("actual_optimizer_histogram", histogram_string(actual_histogram));
    fact("control_code_2", control_2);
    fact("control_code_140", control_140);
    fact("control_code_20", control_20);

    std::ostringstream semantic_payload;
    semantic_payload << "thm4137-n10-independent-semantic-v1\n";
    for (const std::string& line : facts) semantic_payload << line << '\n';

    std::cout << "implementation=all-class gentourng stream + local two-way reachability + independent contracted good/reversed-block engine\n";
    std::cout << "literal_crosscheck=two worst profiles replayed as 1024 explicit 11-vertex child tournaments\n";
    std::cout << "input_contract=45 ASCII upper-triangle bits per gentourng record; all order-ten classes, not -c\n";
    std::cout << "isomorph_rejection_contract=gentourng one representative per unlabeled class; checked by frozen all-class count and raw digest\n";
    std::cout << "profile_serialization=ASCII header thm4137-n10-independent-profile-v1\\n; for each locally strong record in input order: 45 ASCII bits; H:u32le; W2,D4x4,Chdx4,theta_num,theta_den,rho_num,rho_den:i64le; rational,coset,actual cardinality masks:u16le; m=1..9 J_num,J_den,floor,actual_max,lattice:i64le; F(mask):u32le for masks 0..1023\n";
    for (const std::string& line : facts) std::cout << line << '\n';
    const std::string semantic_hash = sha_string(semantic_payload.str());
    std::cout << "semantic_serialization=ASCII header thm4137-n10-independent-semantic-v1\\n followed by the ordered key=value fact lines above\n";
    std::cout << "semantic_sha256=" << semantic_hash << '\n';

    if (all_count == EXPECTED_ALL) {
      const std::map<std::string, uint64_t> expected_actual = {
          {"(-2)", 1550812}, {"(-2,-4)", 321}, {"(-4)", 22353},
          {"(0)", 6186633}, {"(0,-2)", 11172}, {"(2)", 1550812},
          {"(2,0)", 11172}, {"(4)", 22353}, {"(4,2)", 321}};
      if (strong_count != EXPECTED_STRONG ||
          all_raw_hash != EXPECTED_ALL_RAW_SHA256 ||
          strong_raw_hash != EXPECTED_STRONG_RAW_SHA256 ||
          profile_hash != EXPECTED_PROFILE_SHA256 ||
          semantic_hash != EXPECTED_SEMANTIC_SHA256 ||
          rational_failures || coset_failures || coset_reorders ||
          minimum_margin != 24 || minimum_margin_profiles.size() != 2 ||
          theta_zero != 8599 ||
          compare(worst_rho, make_rat(34499248, 37325237)) ||
          worst_profiles.size() != 2 ||
          actual_noncentral_only != 3146972 ||
          rational_histogram != std::map<std::string, uint64_t>{{"(0)", EXPECTED_STRONG}} ||
          coset_histogram != std::map<std::string, uint64_t>{{"(0)", EXPECTED_STRONG}} ||
          actual_histogram != expected_actual)
        throw std::runtime_error("full-universe claimed invariant mismatch");
      for (const Profile& p : worst_profiles) {
        if (p.H != 1431 || p.W2 != 2 * 19557 ||
            p.D4x4 != 4 * 111975711 ||
            std::llabs(p.Chdx4) != 4 * 29570784)
          throw std::runtime_error("worst packet mismatch");
      }
      if (worst_profiles[0].Chdx4 != -worst_profiles[1].Chdx4)
        throw std::runtime_error("worst packets are not opposite-tilt paired");
      std::cout << "status=ACCEPT\n";
    } else {
      std::cout << "status=PARTIAL_INPUT_ONLY\n";
    }
  } catch (const std::exception& error) {
    std::cerr << "status=REJECT error=" << error.what() << '\n';
    return 1;
  }
  return 0;
}
