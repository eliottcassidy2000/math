#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

namespace {

struct Tournament {
  std::size_t n;
  std::vector<u64> out;
};

std::string decimal(i128 x) {
  if (x == 0) return "0";
  bool negative = x < 0;
  if (negative) x = -x;
  std::string answer;
  while (x) {
    answer.push_back(static_cast<char>('0' + x % 10));
    x /= 10;
  }
  if (negative) answer.push_back('-');
  std::reverse(answer.begin(), answer.end());
  return answer;
}

i128 power(i128 base, std::size_t exponent) {
  i128 answer = 1;
  while (exponent--) answer *= base;
  return answer;
}

Tournament family(std::size_t m) {
  if (m < 1 || m > 20) throw std::runtime_error("family scope is 1..20");
  const std::size_t n = m + 5;
  Tournament t{n, std::vector<u64>(n, 0)};
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j) t.out[i] |= u64{1} << j;
  // Reverse 0->3 and 3->z from the transitive order.
  t.out[0] ^= u64{1} << 3;
  t.out[3] |= u64{1} << 0;
  const std::size_t z = n - 1;
  t.out[3] ^= u64{1} << z;
  t.out[z] |= u64{1} << 3;
  return t;
}

std::string encode(const Tournament& t) {
  std::string bits;
  for (std::size_t i = 0; i < t.n; ++i)
    for (std::size_t j = i + 1; j < t.n; ++j)
      bits.push_back((t.out[i] & (u64{1} << j)) ? '1' : '0');
  return bits;
}

std::size_t order_of(const std::string& bits) {
  std::size_t n = 0;
  while (n * (n - 1) / 2 < bits.size()) ++n;
  if (n * (n - 1) / 2 != bits.size())
    throw std::runtime_error("input length is not a tournament length");
  return n;
}

Tournament parse(const std::string& bits) {
  const std::size_t n = order_of(bits);
  if (n < 3 || n > 12) throw std::runtime_error("census scope is 3..12");
  Tournament t{n, std::vector<u64>(n, 0)};
  std::size_t p = 0;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j, ++p) {
      if (bits[p] == '1')
        t.out[i] |= u64{1} << j;
      else if (bits[p] == '0')
        t.out[j] |= u64{1} << i;
      else
        throw std::runtime_error("non-binary tournament word");
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

bool arc(const Tournament& t, std::size_t x, std::size_t y) {
  return t.out[x] & (u64{1} << y);
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
    paths[0] = 1;
    for (std::size_t i = 0; i < t.n; ++i) at(u64{1} << i, i) = 1;
    for (u64 mask = 1; mask < states; ++mask) {
      for (std::size_t last = 0; last < t.n; ++last) {
        if (!(mask & (u64{1} << last))) continue;
        const u64 rest = mask ^ (u64{1} << last);
        if (rest) {
          for (std::size_t previous = 0; previous < t.n; ++previous) {
            if ((rest & (u64{1} << previous)) && arc(t, previous, last))
              at(mask, last) += get(rest, previous);
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
        if (!(mask & (u64{1} << next)) && arc(t, last, next))
          from_zero[index(mask | (u64{1} << next), next)] += ways;
      }
    }
  }
  u64 answer = 0;
  for (std::size_t last = 1; last < t.n; ++last)
    if (arc(t, last, 0)) answer += from_zero[index(full, last)];
  return answer;
}

u64 two_path_cover(const PathDP& dp, std::size_t i, std::size_t j) {
  u64 answer = 0;
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
    for (u64 mask = 1; mask <= dp.full; ++mask)
      if (mask & (u64{1} << i))
        v += dp.get(mask, i) * dp.paths[dp.full ^ mask];
    if (v < h) throw std::runtime_error("v_i below H");
    d[i] = v - h;
  }
  return d;
}

i128 formula_h(std::size_t m) { return 5 * power(2, m) + 5; }
i128 formula_end_r(std::size_t m) { return 5 * power(2, m - 1); }
i128 formula_n_rz(std::size_t m) { return 16 * power(2, m) - 7; }
i128 formula_b_rz(std::size_t m) { return 27 * power(2, m - 1) - 11; }

std::vector<i128> formula_defects(std::size_t m) {
  const std::size_t n = m + 5, z = n - 1, r = n - 2;
  const i128 two_m = power(2, m);
  std::vector<i128> d(n, 0);
  d[1] = two_m + 1;
  d[2] = 5 * two_m + 4;
  d[3] = 3 * two_m + 2;
  for (std::size_t t = 1; t < m; ++t) {
    d[3 + t] = 17 * power(3, t - 1) * power(2, m - t + 1)
                 - 5 * two_m + power(2, t + 4) - 5;
  }
  d[r] = 34 * power(3, m - 1) + 11 * two_m - 12;
  d[z] = 27 * two_m - 12;
  return d;
}

i128 formula_q(std::size_t m) {
  return 34 * power(3, m) + (34 - 5 * static_cast<i128>(m)) * power(2, m)
         - 5 * static_cast<i128>(m) - 44;
}

i128 formula_gap(std::size_t m) {
  const i128 numerator = 4896 * power(6, m) - 8174 * power(4, m)
                         - 476 * power(3, m) + 6504 * power(2, m) + 1298;
  if (numerator % 3) throw std::runtime_error("nonintegral gap formula");
  return numerator / 3;
}

i128 defect_gap(std::size_t n, u64 h, const std::vector<u64>& d) {
  i128 q = 0, squares = 0;
  for (u64 x : d) {
    q += x;
    squares += static_cast<i128>(x) * x;
  }
  return q * q + 2 * static_cast<i128>(n - 4) * h * q
         + static_cast<i128>(n) * (n - 3) * h * h - 5 * squares;
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

void literal_order_ten() {
  const std::size_t m = 5;
  const Tournament t = family(m);
  const std::size_t n = t.n, r = n - 2, z = n - 1;
  std::vector<std::size_t> p(n);
  for (std::size_t i = 0; i < n; ++i) p[i] = i;
  u64 linear_orders = 0, cyclic_orders = 0, h = 0, c = 0, b_rz = 0, n_rz = 0;
  std::vector<u64> end(n, 0), d(n, 0);
  do {
    ++linear_orders;
    std::size_t bad_count = 0, bad_owner = n;
    std::vector<bool> prefix_good(n, true), suffix_good(n, true);
    for (std::size_t k = 0; k + 1 < n; ++k) {
      const bool good = arc(t, p[k], p[k + 1]);
      prefix_good[k + 1] = prefix_good[k] && good;
      if (!good) {
        ++bad_count;
        bad_owner = p[k];
      }
    }
    for (std::size_t k = n - 1; k > 0; --k)
      suffix_good[k - 1] = suffix_good[k] && arc(t, p[k - 1], p[k]);
    if (bad_count == 0) {
      ++h;
      ++end[p.back()];
    } else if (bad_count == 1) {
      ++d[bad_owner];
    }

    // Concatenate the r-ending path first and the z-ending path second.
    // This labels the two otherwise unordered components uniquely.
    if (p.back() == z) {
      for (std::size_t cut = 1; cut < n; ++cut)
        if (p[cut - 1] == r && prefix_good[cut - 1] && suffix_good[cut]) ++n_rz;
    }

    if (p[0] == 0) {
      ++cyclic_orders;
      u64 bad_mask = 0;
      for (std::size_t k = 0; k < n; ++k)
        if (!arc(t, p[k], p[(k + 1) % n])) bad_mask |= u64{1} << p[k];
      if (bad_mask == 0) ++c;
      if (bad_mask == ((u64{1} << r) | (u64{1} << z))) ++b_rz;
    }
  } while (std::next_permutation(p.begin(), p.end()));

  const std::vector<i128> expected_d = formula_defects(m);
  if (linear_orders != 3628800 || cyclic_orders != 362880 || h != formula_h(m)
      || c != 1 || end[r] != formula_end_r(m) || end[z] != 5
      || n_rz != formula_n_rz(m) || b_rz != formula_b_rz(m))
    throw std::runtime_error("literal order-ten scalar mismatch");
  for (std::size_t i = 0; i < n; ++i)
    if (d[i] != expected_d[i]) throw std::runtime_error("literal defect mismatch");

  std::cout << "mode=literal_order_ten\n";
  std::cout << "encoding=upper_triangle_row_major;1_means_i_to_j\n";
  std::cout << "word=" << encode(t) << "\n";
  std::cout << "strong=" << (is_strong(t) ? "yes" : "no")
            << " linear_orders=" << linear_orders << " cyclic_orders=" << cyclic_orders
            << "\n";
  std::cout << "H=" << h << " c=" << c << " N_8,9=" << n_rz
            << " B_8,9=" << b_rz << "\n";
  print_vector("End", end);
  print_vector("D", d);
  std::cout << "product_plus_min=" << end[r] * end[z] + std::min(end[r], end[z])
            << " repaired_excess="
            << decimal(static_cast<i128>(b_rz) - end[r] * end[z]
                       - std::min(end[r], end[z]))
            << "\n";
  std::cout << "five_copy_gap=" << decimal(defect_gap(n, h, d)) << "\n";
  std::cout << "checks=PASS literal_10_factorial,cyclic_9_factorial,formula_profile\n";
}

void family_dp() {
  std::cout << "mode=family_subset_DP range_m=1..12\n";
  for (std::size_t m = 1; m <= 12; ++m) {
    const Tournament t = family(m);
    const std::size_t n = t.n, r = n - 2, z = n - 1;
    PathDP dp(t);
    const u64 h = dp.paths[dp.full];
    const u64 c = cycle_count(t);
    const u64 n_rz = two_path_cover(dp, r, z);
    const u64 b_rz = n_rz + c - dp.get(dp.full, r) - dp.get(dp.full, z);
    const std::vector<u64> d = defect_vector(dp);
    const std::vector<i128> expected_d = formula_defects(m);
    if (!is_strong(t) || h != formula_h(m) || c != 1
        || dp.get(dp.full, r) != formula_end_r(m) || dp.get(dp.full, z) != 5
        || n_rz != formula_n_rz(m) || b_rz != formula_b_rz(m))
      throw std::runtime_error("family scalar formula mismatch");
    i128 q = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (d[i] != expected_d[i]) throw std::runtime_error("family D formula mismatch");
      q += d[i];
    }
    const i128 gap = defect_gap(n, h, d);
    if (q != formula_q(m) || gap != formula_gap(m))
      throw std::runtime_error("family q/gap formula mismatch");
    std::cout << "m=" << m << " n=" << n << " H=" << h << " c=" << c
              << " End_r=" << dp.get(dp.full, r) << " End_z=" << dp.get(dp.full, z)
              << " N_rz=" << n_rz << " B_rz=" << b_rz
              << " repair_excess="
              << decimal(static_cast<i128>(b_rz)
                         - static_cast<i128>(dp.get(dp.full, r)) * dp.get(dp.full, z)
                         - std::min(dp.get(dp.full, r), dp.get(dp.full, z)))
              << " ratio=" << b_rz << '/'
              << dp.get(dp.full, r) * dp.get(dp.full, z)
              << " q=" << decimal(q) << " G=" << decimal(gap) << "\n";
  }
  std::cout << "checks=PASS strong,H,c,End,N,B,D,q,G_all_rows\n";
}

void abstract_firewall() {
  constexpr std::size_t n = 6;
  const u64 full = (u64{1} << n) - 1;
  std::vector<u64> counts(u64{1} << n, 0);
  counts[0] = 1;
  counts[u64{1} << 0] = 5;
  counts[u64{1} << 1] = 5;
  counts[(u64{1} << 0) | (u64{1} << 1)] = 25;
  counts[full] = 1;
  counts[full ^ (u64{1} << 0)] = 5;
  counts[full ^ (u64{1} << 1)] = 5;
  counts[full ^ (u64{1} << 0) ^ (u64{1} << 1)] = 25;
  counts[(u64{1} << 0) | (u64{1} << 1) | (u64{1} << 2)] = 12;
  counts[(u64{1} << 3) | (u64{1} << 4) | (u64{1} << 5)] = 36;

  std::vector<u64> layers(n + 1, 0), marginals(n, 0);
  u64 total = 0, first_moment = 0;
  for (u64 mask = 0; mask <= full; ++mask) {
    const std::size_t size = static_cast<std::size_t>(__builtin_popcountll(mask));
    layers[size] += counts[mask];
    total += counts[mask];
    first_moment += size * counts[mask];
    for (std::size_t i = 0; i < n; ++i)
      if (mask & (u64{1} << i)) marginals[i] += counts[mask];
  }
  const std::vector<u64> expected_layers{1, 10, 25, 48, 25, 10, 1};
  const std::vector<u64> expected_marginals{48, 48, 48, 72, 72, 72};
  const std::vector<u64> indegrees{2, 2, 2, 3, 3, 3};
  bool strict_landau = true;
  u64 degree_sum = 0;
  for (std::size_t k = 1; k < n; ++k) {
    degree_sum += indegrees[k - 1];
    if (degree_sum <= k * (k - 1) / 2) strict_landau = false;
  }
  degree_sum += indegrees.back();

  const u64 c = counts[0], a0 = counts[1], a1 = counts[2], b01 = counts[3];
  const u64 h = n * c + a0 + a1;
  std::vector<u64> d(n, 0);
  d[0] = (n - 1) * a0 + b01;
  d[1] = (n - 1) * a1 + b01;
  const i128 gap = defect_gap(n, h, d);
  if (total != 120 || layers != expected_layers || first_moment != 360
      || marginals != expected_marginals || !strict_landau || degree_sum != 15
      || h != 16 || d[0] != 50 || d[1] != 50 || gap != -3992
      || 25 * b01 > 27 * (c + a0) * (c + a1))
    throw std::runtime_error("abstract firewall mismatch");

  std::cout << "mode=abstract_moment_firewall\n";
  print_vector("layers", layers);
  print_vector("owner_marginals", marginals);
  std::cout << "compatible_indegrees=(2,2,2,3,3,3) strict_Landau=yes\n";
  std::cout << "c=" << c << " A_0=A_1=" << a0 << " B_0,1=" << b01
            << " End_0=End_1=" << c + a0 << " H=" << h
            << " D=(50,50,0,0,0,0) five_copy_gap=" << decimal(gap) << "\n";
  std::cout << "checks=PASS total,reversal_layers,first_moment,owner_marginals,"
               "strong_score_feasibility,27_over_25_local_bound\n";
}

void census() {
  std::string bits, best_word;
  std::size_t rows = 0, order = 0, best_i = 0, best_j = 0, failures_27_25 = 0;
  u64 best_b = 0, best_product = 1;
  while (std::cin >> bits) {
    const Tournament t = parse(bits);
    if (!is_strong(t)) throw std::runtime_error("non-strong census input");
    if (!order) order = t.n;
    if (order != t.n) throw std::runtime_error("mixed census orders");
    PathDP dp(t);
    const u64 c = cycle_count(t);
    std::vector<u64> end(t.n);
    for (std::size_t i = 0; i < t.n; ++i) end[i] = dp.get(dp.full, i);
    for (std::size_t i = 0; i < t.n; ++i) {
      for (std::size_t j = i + 1; j < t.n; ++j) {
        const u64 nij = two_path_cover(dp, i, j);
        if (nij + c < end[i] + end[j]) throw std::runtime_error("negative B");
        const u64 b = nij + c - end[i] - end[j];
        const u64 product = end[i] * end[j];
        if (static_cast<i128>(b) * best_product
            > static_cast<i128>(best_b) * product) {
          best_b = b;
          best_product = product;
          best_word = bits;
          best_i = i;
          best_j = j;
        }
        if (25 * static_cast<i128>(b) > 27 * static_cast<i128>(product))
          ++failures_27_25;
      }
    }
    ++rows;
  }
  const long double ratio = static_cast<long double>(best_b) / best_product;
  std::cout << "order=" << order << " strong_classes=" << rows
            << " max_ratio=" << best_b << '/' << best_product
            << " decimal=" << std::fixed << std::setprecision(12) << ratio
            << " candidate_25B_le_27product_failures=" << failures_27_25 << "\n";
  std::cout << "max_word=" << best_word << " pair=" << best_i << ',' << best_j << "\n";
  std::cout << "checks=PASS strong_input,nonnegative_B,exact_cross_multiplication\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 1) {
      literal_order_ten();
      family_dp();
      abstract_firewall();
    } else if (argc == 2 && std::string(argv[1]) == "--census") {
      census();
    } else {
      std::cerr << "usage: " << argv[0] << " [--census]\n";
      return 2;
    }
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
