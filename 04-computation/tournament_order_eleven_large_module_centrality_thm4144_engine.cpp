#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>

namespace {

constexpr int N = 11;
constexpr int SUBSETS = 1 << N;
constexpr int FULL = SUBSETS - 1;
constexpr int EDGES = N * (N - 1) / 2;

struct Profile {
  bool strong = false;
  int64_t H = 0;
  int64_t twice_W = 0;
  int64_t four_D4 = 0;
  int64_t four_C_hd = 0;
  int64_t coset_margin = std::numeric_limits<int64_t>::max();
};

int64_t absolute(int64_t value) { return value < 0 ? -value : value; }

int64_t ceiling_division(__int128 numerator, __int128 denominator) {
  if (denominator <= 0) {
    std::cerr << "nonpositive ceiling denominator\n";
    std::exit(3);
  }
  __int128 quotient = numerator / denominator;
  __int128 remainder = numerator % denominator;
  if (remainder > 0) ++quotient;
  return static_cast<int64_t>(quotient);
}

std::array<uint16_t, N> adjacency_from_label(const std::string& label) {
  if (static_cast<int>(label.size()) != EDGES) {
    std::cerr << "bad label length\n";
    std::exit(2);
  }
  std::array<uint16_t, N> outgoing{};
  int edge = 0;
  for (int left = 0; left < N; ++left) {
    for (int right = left + 1; right < N; ++right, ++edge) {
      if (label[edge] == '1') {
        outgoing[left] |= uint16_t(1u << right);
      } else if (label[edge] == '0') {
        outgoing[right] |= uint16_t(1u << left);
      } else {
        std::cerr << "nonbinary label\n";
        std::exit(2);
      }
    }
  }
  return outgoing;
}

bool is_strong(const std::array<uint16_t, N>& outgoing) {
  std::array<uint16_t, N> incoming{};
  for (int source = 0; source < N; ++source) {
    for (int target = 0; target < N; ++target) {
      if (outgoing[source] & (1u << target)) {
        incoming[target] |= uint16_t(1u << source);
      }
    }
  }
  for (int reverse = 0; reverse < 2; ++reverse) {
    uint16_t seen = 1;
    uint16_t todo = 1;
    while (todo) {
      int vertex = __builtin_ctz(todo);
      todo &= uint16_t(todo - 1);
      uint16_t next = (reverse ? incoming[vertex] : outgoing[vertex])
                      & uint16_t(~seen);
      seen |= next;
      todo |= next;
    }
    if ((seen & FULL) != FULL) return false;
  }
  return true;
}

Profile evaluate(const std::string& label) {
  const auto outgoing = adjacency_from_label(label);
  Profile profile;
  profile.strong = is_strong(outgoing);
  if (!profile.strong) return profile;

  std::array<uint16_t, N> incoming{};
  for (int source = 0; source < N; ++source) {
    for (int target = 0; target < N; ++target) {
      if (outgoing[source] & (1u << target)) {
        incoming[target] |= uint16_t(1u << source);
      }
    }
  }

  static uint32_t ending[SUBSETS][N];
  static uint32_t starting[SUBSETS][N];
  static uint32_t before[SUBSETS][N];
  static uint32_t after[SUBSETS][N];
  std::fill(&ending[0][0], &ending[0][0] + SUBSETS * N, 0);
  std::fill(&starting[0][0], &starting[0][0] + SUBSETS * N, 0);
  std::fill(&before[0][0], &before[0][0] + SUBSETS * N, 0);
  std::fill(&after[0][0], &after[0][0] + SUBSETS * N, 0);

  for (int vertex = 0; vertex < N; ++vertex) {
    ending[1 << vertex][vertex] = 1;
    starting[1 << vertex][vertex] = 1;
  }
  for (int mask = 1; mask < SUBSETS; ++mask) {
    uint16_t vertices = uint16_t(mask);
    while (vertices) {
      int vertex = __builtin_ctz(vertices);
      vertices &= uint16_t(vertices - 1);
      int rest = mask ^ (1 << vertex);
      if (!rest) continue;
      uint16_t choices = uint16_t(rest) & incoming[vertex];
      while (choices) {
        int previous = __builtin_ctz(choices);
        choices &= uint16_t(choices - 1);
        ending[mask][vertex] += ending[rest][previous];
      }
      choices = uint16_t(rest) & outgoing[vertex];
      while (choices) {
        int next = __builtin_ctz(choices);
        choices &= uint16_t(choices - 1);
        starting[mask][vertex] += starting[rest][next];
      }
    }
  }
  for (int vertex = 0; vertex < N; ++vertex) {
    profile.H += ending[FULL][vertex];
    before[0][vertex] = 1;
    after[0][vertex] = 1;
  }
  for (int mask = 1; mask < SUBSETS; ++mask) {
    for (int vertex = 0; vertex < N; ++vertex) {
      uint16_t choices = uint16_t(mask) & incoming[vertex];
      while (choices) {
        int previous = __builtin_ctz(choices);
        choices &= uint16_t(choices - 1);
        before[mask][vertex] += ending[mask][previous];
      }
      choices = uint16_t(mask) & outgoing[vertex];
      while (choices) {
        int next = __builtin_ctz(choices);
        choices &= uint16_t(choices - 1);
        after[mask][vertex] += starting[mask][next];
      }
    }
  }

  uint32_t capacity[N][N]{};
  int64_t twice_degree[N]{};
  int64_t twice_field[N]{};
  for (int left = 0; left < N; ++left) {
    for (int right = left + 1; right < N; ++right) {
      int source = left;
      int target = right;
      if (!(outgoing[source] & (1u << target))) std::swap(source, target);
      int remainder = FULL ^ (1 << source) ^ (1 << target);
      uint64_t count = 0;
      for (int first = remainder;; first = (first - 1) & remainder) {
        int second = remainder ^ first;
        count += uint64_t(before[first][source]) * after[second][target];
        count += uint64_t(before[first][target]) * after[second][source];
        if (!first) break;
      }
      capacity[source][target] = capacity[target][source] = uint32_t(count);
      profile.twice_W += int64_t(count);
    }
  }
  for (int vertex = 0; vertex < N; ++vertex) {
    for (int other = 0; other < N; ++other) {
      if (vertex == other) continue;
      twice_degree[vertex] += capacity[vertex][other];
      twice_field[vertex] += (outgoing[vertex] & (1u << other))
          ? int64_t(capacity[vertex][other])
          : -int64_t(capacity[vertex][other]);
    }
    profile.four_C_hd += twice_field[vertex] * twice_degree[vertex];
  }
  for (int a = 0; a < N; ++a) {
    for (int b = a + 1; b < N; ++b) {
      for (int c = a + 1; c < N; ++c) {
        if (c == b) continue;
        for (int d = c + 1; d < N; ++d) {
          if (d == b) continue;
          profile.four_D4 += int64_t(capacity[a][b]) * capacity[c][d];
        }
      }
    }
  }

  std::array<int64_t, SUBSETS> values{};
  for (int mask = 0; mask < SUBSETS; ++mask) {
    int64_t value = profile.H;
    for (int source = 0; source < N; ++source) {
      if (!(mask & (1 << source))) continue;
      for (int target = 0; target < N; ++target) {
        if ((mask & (1 << target)) || !(outgoing[source] & (1u << target))) continue;
        value += capacity[source][target];
      }
    }
    values[mask] = value;
  }

  int64_t best_central = std::numeric_limits<int64_t>::min();
  int64_t best_outer = std::numeric_limits<int64_t>::min();
  for (int size = 1; size < N; ++size) {
    int64_t count = 0;
    int64_t sum = 0;
    int64_t square_sum = 0;
    int64_t anchor = -1;
    int64_t lattice = 0;
    for (int mask = 0; mask < SUBSETS; ++mask) {
      if (__builtin_popcount(unsigned(mask)) != size) continue;
      int64_t value = values[mask];
      if (anchor < 0) anchor = value;
      ++count;
      sum += value;
      square_sum += value * value;
      lattice = std::gcd(lattice, absolute(value - anchor));
    }
    int64_t delta = sum - count * profile.H;
    if (delta <= 0) {
      std::cerr << "nonpositive layer mean gap\n";
      std::exit(4);
    }
    int64_t floor = anchor;
    if (lattice) {
      __int128 numerator = __int128(square_sum) - __int128(sum) * profile.H
                           - __int128(anchor) * delta;
      floor += lattice * ceiling_division(numerator, __int128(delta) * lattice);
    }
    if (size == 5 || size == 6) {
      best_central = std::max(best_central, floor);
    } else {
      best_outer = std::max(best_outer, floor);
    }
  }
  profile.coset_margin = best_central - best_outer;
  return profile;
}

bool larger_rho(const Profile& left, const Profile& right) {
  return __int128(absolute(left.four_C_hd)) * right.four_D4
       > __int128(absolute(right.four_C_hd)) * left.four_D4;
}

void print_profile(const char* prefix, const std::string& label,
                   const Profile& profile) {
  if (profile.twice_W % 2 || profile.four_D4 % 4 || profile.four_C_hd % 4) {
    std::cerr << "nonintegral normalized exposure packet\n";
    std::exit(5);
  }
  const int64_t W = profile.twice_W / 2;
  const int64_t D4 = profile.four_D4 / 4;
  const int64_t C_hd = profile.four_C_hd / 4;
  int64_t numerator = 2 * absolute(C_hd);
  int64_t denominator = D4;
  int64_t common = std::gcd(numerator, denominator);
  numerator /= common;
  denominator /= common;
  std::cout << prefix << "_label=" << label
            << ";packet=(" << profile.H << ',' << W << ',' << D4 << ',' << C_hd << ')'
            << ";rho=" << numerator << '/' << denominator
            << ";coset_margin=" << profile.coset_margin << '\n';
}

}  // namespace

int main() {
  std::string label;
  uint64_t rows = 0;
  uint64_t strong_rows = 0;
  uint64_t rational_failures = 0;
  uint64_t coset_failures = 0;
  bool have_maximum = false;
  bool have_minimum = false;
  Profile maximum;
  Profile minimum;
  std::string maximum_label;
  std::string minimum_label;

  while (std::cin >> label) {
    ++rows;
    Profile profile = evaluate(label);
    if (!profile.strong) continue;
    ++strong_rows;
    if (__int128(2) * absolute(profile.four_C_hd) >= profile.four_D4) {
      ++rational_failures;
    }
    if (profile.coset_margin <= 0) ++coset_failures;
    if (!have_maximum || larger_rho(profile, maximum)) {
      have_maximum = true;
      maximum = profile;
      maximum_label = label;
    }
    if (!have_minimum || profile.coset_margin < minimum.coset_margin) {
      have_minimum = true;
      minimum = profile;
      minimum_label = label;
    }
  }
  if (!have_maximum || !have_minimum) {
    std::cerr << "empty strong stream\n";
    return 6;
  }
  std::cout << "rows=" << rows << ";strong=" << strong_rows
            << ";rational_failures=" << rational_failures
            << ";coset_failures=" << coset_failures << '\n';
  print_profile("max_rho", maximum_label, maximum);
  print_profile("min_coset", minimum_label, minimum);
  return 0;
}
