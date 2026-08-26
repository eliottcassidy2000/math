#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;

namespace {

constexpr const char* kWitness =
    "110111111111111111111111011111111110";

struct Tournament {
  std::size_t n;
  std::vector<u64> out;
};

std::size_t order_of(const std::string& bits) {
  std::size_t n = 0;
  while (n * (n - 1) / 2 < bits.size()) ++n;
  if (n * (n - 1) / 2 != bits.size())
    throw std::runtime_error("input length is not a tournament length");
  return n;
}

Tournament parse(const std::string& bits) {
  const std::size_t n = order_of(bits);
  if (n < 3 || n > 12) throw std::runtime_error("supported orders are 3..12");
  Tournament t{n, std::vector<u64>(n, 0)};
  std::size_t p = 0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j, ++p) {
      if (bits[p] == '1')
        t.out[i] |= u64{1} << j;
      else if (bits[p] == '0')
        t.out[j] |= u64{1} << i;
      else
        throw std::runtime_error("input contains a non-binary character");
    }
  }
  return t;
}

bool is_strong(const Tournament& t) {
  const u64 full = (u64{1} << t.n) - 1;
  for (int reverse = 0; reverse <= 1; ++reverse) {
    u64 seen = 1;
    while (true) {
      u64 next = seen;
      for (std::size_t i = 0; i < t.n; ++i) {
        if (!(seen & (u64{1} << i))) continue;
        if (!reverse) {
          next |= t.out[i];
        } else {
          for (std::size_t j = 0; j < t.n; ++j)
            if (t.out[j] & (u64{1} << i)) next |= u64{1} << j;
        }
      }
      if (next == seen) break;
      seen = next;
    }
    if (seen != full) return false;
  }
  return true;
}

struct PathDP {
  const Tournament& t;
  u64 states;
  u64 full;
  std::vector<u64> finish;
  std::vector<u64> paths;

  explicit PathDP(const Tournament& tournament)
      : t(tournament),
        states(u64{1} << t.n),
        full(states - 1),
        finish(static_cast<std::size_t>(states) * t.n, 0),
        paths(states, 0) {
    paths[0] = 1;  // the empty complementary path
    for (std::size_t i = 0; i < t.n; ++i) at(u64{1} << i, i) = 1;
    for (u64 mask = 1; mask < states; ++mask) {
      for (std::size_t last = 0; last < t.n; ++last) {
        if (!(mask & (u64{1} << last))) continue;
        const u64 rest = mask ^ (u64{1} << last);
        if (rest) {
          for (std::size_t prev = 0; prev < t.n; ++prev) {
            if ((rest & (u64{1} << prev)) &&
                (t.out[prev] & (u64{1} << last)))
              at(mask, last) += get(rest, prev);
          }
        }
        paths[mask] += get(mask, last);
      }
    }
  }

  std::size_t index(u64 mask, std::size_t i) const {
    return static_cast<std::size_t>(mask) * t.n + i;
  }
  u64 get(u64 mask, std::size_t i) const { return finish[index(mask, i)]; }
  u64& at(u64 mask, std::size_t i) { return finish[index(mask, i)]; }
};

u64 cycle_count(const Tournament& t) {
  const u64 states = u64{1} << t.n;
  const u64 full = states - 1;
  const auto index = [&t](u64 mask, std::size_t i) {
    return static_cast<std::size_t>(mask) * t.n + i;
  };
  std::vector<u64> from_zero(static_cast<std::size_t>(states) * t.n, 0);
  from_zero[index(1, 0)] = 1;
  for (u64 mask = 1; mask < states; ++mask) {
    if (!(mask & 1)) continue;
    for (std::size_t last = 0; last < t.n; ++last) {
      const u64 ways = from_zero[index(mask, last)];
      if (!ways) continue;
      for (std::size_t next = 0; next < t.n; ++next) {
        if (!(mask & (u64{1} << next)) &&
            (t.out[last] & (u64{1} << next)))
          from_zero[index(mask | (u64{1} << next), next)] += ways;
      }
    }
  }
  u64 answer = 0;
  for (std::size_t last = 1; last < t.n; ++last)
    if (t.out[last] & 1) answer += from_zero[index(full, last)];
  return answer;
}

u64 two_path_cover(const PathDP& dp, std::size_t i, std::size_t j) {
  u64 answer = 0;
  // The component ending at i contains i and not j.  This selects exactly
  // one of the two unordered components, so there is no factor of two.
  for (u64 mask = 1; mask < dp.full; ++mask) {
    if ((mask & (u64{1} << i)) && !(mask & (u64{1} << j)))
      answer += dp.get(mask, i) * dp.get(dp.full ^ mask, j);
  }
  return answer;
}

std::vector<u64> defect_vector(const PathDP& dp) {
  const u64 h = dp.paths[dp.full];
  std::vector<u64> d(dp.t.n, 0);
  for (std::size_t i = 0; i < dp.t.n; ++i) {
    u64 v = 0;
    for (u64 mask = 1; mask <= dp.full; ++mask) {
      if (mask & (u64{1} << i))
        v += dp.get(mask, i) * dp.paths[dp.full ^ mask];
    }
    if (v < h) throw std::runtime_error("split count below H");
    d[i] = v - h;
  }
  return d;
}

template <class T>
void print_vector(const std::string& name, const std::vector<T>& values) {
  std::cout << name << "=(";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) std::cout << ',';
    std::cout << values[i];
  }
  std::cout << ")\n";
}

void direct_witness() {
  const std::string word = kWitness;
  const Tournament t = parse(word);
  PathDP dp(t);
  const u64 h_dp = dp.paths[dp.full];
  const u64 c_dp = cycle_count(t);
  const std::vector<u64> d_dp = defect_vector(dp);
  std::vector<u64> end_dp(t.n, 0);
  for (std::size_t i = 0; i < t.n; ++i) end_dp[i] = dp.get(dp.full, i);

  const u64 states = u64{1} << t.n;
  std::vector<u64> exact_bad(states, 0), cut_cover(states, 0);
  std::vector<u64> end_direct(t.n, 0), d_direct(t.n, 0);
  u64 h_direct = 0;
  std::vector<std::size_t> p(t.n);
  for (std::size_t i = 0; i < t.n; ++i) p[i] = i;
  do {
    std::vector<std::size_t> linear_bad;
    for (std::size_t k = 0; k + 1 < t.n; ++k) {
      if (!(t.out[p[k]] & (u64{1} << p[k + 1])))
        linear_bad.push_back(p[k]);
    }
    if (linear_bad.empty()) {
      ++h_direct;
      ++end_direct[p.back()];
    } else if (linear_bad.size() == 1) {
      ++d_direct[linear_bad[0]];
    }

    // Fixing vertex 0 first chooses one representative modulo rotation.
    if (p[0] == 0) {
      u64 bad = 0;
      for (std::size_t k = 0; k < t.n; ++k) {
        const std::size_t next = (k + 1) % t.n;
        if (!(t.out[p[k]] & (u64{1} << p[next]))) bad |= u64{1} << p[k];
      }
      ++exact_bad[bad];
      for (u64 r = 0; r < states; ++r)
        if ((bad & ~r) == 0) ++cut_cover[r];
    }
  } while (std::next_permutation(p.begin(), p.end()));

  if (h_direct != h_dp || end_direct != end_dp || d_direct != d_dp)
    throw std::runtime_error("literal permutation scan disagrees with path DP");

  std::vector<i64> inverse(cut_cover.begin(), cut_cover.end());
  for (std::size_t bit = 0; bit < t.n; ++bit) {
    for (u64 mask = 0; mask < states; ++mask) {
      if (mask & (u64{1} << bit)) inverse[mask] -= inverse[mask ^ (u64{1} << bit)];
    }
  }
  for (u64 mask = 0; mask < states; ++mask) {
    if (inverse[mask] < 0 || static_cast<u64>(inverse[mask]) != exact_bad[mask])
      throw std::runtime_error("Boolean Mobius inversion failed");
  }

  const u64 c = exact_bad[0];
  if (c != c_dp) throw std::runtime_error("literal cycle scan disagrees with cycle DP");
  std::vector<u64> a(t.n), b_degree(t.n, 0);
  for (std::size_t i = 0; i < t.n; ++i) a[i] = exact_bad[u64{1} << i];
  for (std::size_t i = 0; i < t.n; ++i) {
    for (std::size_t j = i + 1; j < t.n; ++j) {
      const u64 b = exact_bad[(u64{1} << i) | (u64{1} << j)];
      b_degree[i] += b;
      b_degree[j] += b;
    }
  }
  u64 h_formula = t.n * c;
  for (u64 x : a) h_formula += x;
  if (h_formula != h_direct) throw std::runtime_error("H formula failed");
  for (std::size_t i = 0; i < t.n; ++i) {
    if (end_direct[i] != c + a[i]) throw std::runtime_error("End formula failed");
    if (d_direct[i] != (t.n - 1) * a[i] + b_degree[i])
      throw std::runtime_error("D formula failed");
  }

  const auto pair_mask = [](std::size_t i, std::size_t j) {
    return (u64{1} << i) | (u64{1} << j);
  };
  const u64 n67 = cut_cover[pair_mask(6, 7)];
  const u64 n78 = cut_cover[pair_mask(7, 8)];
  const u64 b67 = exact_bad[pair_mask(6, 7)];
  const u64 b78 = exact_bad[pair_mask(7, 8)];
  if (n67 != two_path_cover(dp, 6, 7) || n78 != two_path_cover(dp, 7, 8))
    throw std::runtime_error("literal cuts disagree with two-path-cover DP");

  std::vector<std::size_t> outdegree(t.n, 0);
  for (std::size_t i = 0; i < t.n; ++i)
    for (std::size_t j = 0; j < t.n; ++j)
      if (t.out[i] & (u64{1} << j)) ++outdegree[i];

  std::cout << "mode=direct_witness\n";
  std::cout << "encoding=upper_triangle_row_major;1_means_i_to_j\n";
  std::cout << "word=" << word << "\n";
  std::cout << "order=" << t.n << " strong=" << (is_strong(t) ? "yes" : "no")
            << " cyclic_orders=" << cut_cover[states - 1] << "\n";
  std::cout << "H=" << h_direct << " c=" << c << "\n";
  print_vector("End", end_direct);
  print_vector("A", a);
  print_vector("D", d_direct);
  print_vector("outdegree", outdegree);
  std::cout << "pair=6,7 N=" << n67 << " B=" << b67
            << " End_product=" << end_direct[6] * end_direct[7]
            << " product_plus_min="
            << end_direct[6] * end_direct[7] + std::min(end_direct[6], end_direct[7])
            << "\n";
  std::cout << "pair=7,8 N=" << n78 << " B=" << b78
            << " End_product=" << end_direct[7] * end_direct[8]
            << " product_plus_min="
            << end_direct[7] * end_direct[8] + std::min(end_direct[7], end_direct[8])
            << "\n";
  std::cout << "checks=PASS permutation_vs_DP,all_exact_formulas,"
               "boolean_mobius_512_subsets\n";
}

void census() {
  std::string bits;
  std::size_t rows = 0;
  std::size_t order = 0;
  std::size_t product_failures = 0, repaired_failures = 0;
  std::size_t product_equal_positive = 0, repaired_equal_positive = 0;
  u64 worst_excess = 0;
  while (std::cin >> bits) {
    const Tournament t = parse(bits);
    if (!is_strong(t)) throw std::runtime_error("input contains a non-strong tournament");
    if (!order) order = t.n;
    if (order != t.n) throw std::runtime_error("mixing orders in one census is unsupported");
    PathDP dp(t);
    const u64 h = dp.paths[dp.full];
    const u64 c = cycle_count(t);
    const std::vector<u64> d = defect_vector(dp);
    std::vector<u64> end(t.n), a(t.n), b_degree(t.n, 0);
    for (std::size_t i = 0; i < t.n; ++i) {
      end[i] = dp.get(dp.full, i);
      if (end[i] < c) throw std::runtime_error("negative singleton bad-owner count");
      a[i] = end[i] - c;
    }
    u64 sum_b_edges = 0;
    for (std::size_t i = 0; i < t.n; ++i) {
      for (std::size_t j = i + 1; j < t.n; ++j) {
        const u64 nij = two_path_cover(dp, i, j);
        if (nij + c < end[i] + end[j])
          throw std::runtime_error("negative two-bad-owner count");
        const u64 bij = nij + c - end[i] - end[j];
        b_degree[i] += bij;
        b_degree[j] += bij;
        sum_b_edges += bij;

        const u64 product = end[i] * end[j];
        const u64 repaired = product + std::min(end[i], end[j]);
        if (bij == product && bij) ++product_equal_positive;
        if (bij == repaired && bij) ++repaired_equal_positive;
        if (bij > repaired) ++repaired_failures;
        if (bij > product) {
          ++product_failures;
          const u64 excess = bij - product;
          worst_excess = std::max(worst_excess, excess);
          if (product_failures <= 20)
            std::cout << "product_failure word=" << bits << " pair=" << i << ',' << j
                      << " c=" << c << " End_i=" << end[i] << " End_j=" << end[j]
                      << " N=" << nij << " B=" << bij << " product=" << product
                      << " product_plus_min=" << repaired << "\n";
        }
      }
    }

    u64 h_formula = t.n * c;
    for (u64 x : a) h_formula += x;
    if (h_formula != h) throw std::runtime_error("H formula failed");
    u64 q_left = 0, q_right = 2 * sum_b_edges;
    for (std::size_t i = 0; i < t.n; ++i) {
      if (d[i] != (t.n - 1) * a[i] + b_degree[i])
        throw std::runtime_error("D formula failed");
      q_left += d[i];
      q_right += (t.n - 1) * a[i];
    }
    if (q_left != q_right) throw std::runtime_error("q formula failed");
    ++rows;
  }
  std::cout << "order=" << order << " strong_classes=" << rows
            << " product_failures=" << product_failures
            << " worst_product_excess=" << worst_excess
            << " repaired_failures=" << repaired_failures
            << " product_equal_positive=" << product_equal_positive
            << " repaired_equal_positive=" << repaired_equal_positive << "\n";
  std::cout << "checks=PASS strong_input,H_End_D_q_formulas,nonnegative_B\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--direct-witness") {
      direct_witness();
    } else if (argc == 1) {
      census();
    } else {
      std::cerr << "usage: " << argv[0] << " [--direct-witness]\n";
      return 2;
    }
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
