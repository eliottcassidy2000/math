#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using i64 = std::int64_t;
using Matrix = std::vector<std::vector<i64>>;

void need(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

struct Tournament {
  int n{};
  std::vector<std::uint32_t> out;
};

Tournament from_mask(int n, std::uint64_t mask) {
  Tournament t{n, std::vector<std::uint32_t>(static_cast<std::size_t>(n), 0)};
  int bit = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j, ++bit) {
      if ((mask >> bit) & 1ULL) t.out[static_cast<std::size_t>(i)] |= 1U << j;
      else t.out[static_cast<std::size_t>(j)] |= 1U << i;
    }
  }
  return t;
}

Tournament parse(const std::string& bits) {
  int n = 0;
  while (n * (n - 1) / 2 < static_cast<int>(bits.size())) ++n;
  need(n * (n - 1) / 2 == static_cast<int>(bits.size()), "pair-bit label");
  std::uint64_t mask = 0;
  for (std::size_t i = 0; i < bits.size(); ++i) {
    if (bits[i] == '1') mask |= 1ULL << i;
  }
  return from_mask(n, mask);
}

std::string label(const Tournament& t) {
  std::string ans;
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j) {
      ans.push_back((t.out[static_cast<std::size_t>(i)] & (1U << j)) ? '1' : '0');
    }
  }
  return ans;
}

bool arc(const Tournament& t, int u, int v) {
  return (t.out[static_cast<std::size_t>(u)] & (1U << v)) != 0;
}

bool has_sink(const Tournament& t) {
  return std::any_of(t.out.begin(), t.out.end(), [](std::uint32_t row) { return row == 0; });
}

struct Data {
  Tournament t;
  i64 h{};
  Matrix c;
  std::vector<std::array<i64, 2>> starts;
  std::vector<std::array<i64, 2>> ends;
  i64 mass{};
  i64 gate{};
  std::vector<i64> even_component_objects;
};

i64 gate(const Tournament& t, const Matrix& c) {
  std::vector<i64> degree(static_cast<std::size_t>(t.n), 0);
  std::vector<i64> current(static_cast<std::size_t>(t.n), 0);
  i64 mass = 0;
  i64 squares = 0;
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j) {
      const i64 value = c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      mass += value;
      squares += value * value;
      degree[static_cast<std::size_t>(i)] += value;
      degree[static_cast<std::size_t>(j)] += value;
      if (arc(t, i, j)) {
        current[static_cast<std::size_t>(i)] += value;
        current[static_cast<std::size_t>(j)] -= value;
      } else {
        current[static_cast<std::size_t>(i)] -= value;
        current[static_cast<std::size_t>(j)] += value;
      }
    }
  }
  i64 C = 0;
  i64 degree_squares = 0;
  for (int i = 0; i < t.n; ++i) {
    C += degree[static_cast<std::size_t>(i)] * current[static_cast<std::size_t>(i)];
    degree_squares += degree[static_cast<std::size_t>(i)]
                      * degree[static_cast<std::size_t>(i)];
  }
  const i64 D = (mass * mass + squares - degree_squares) / 2;
  return D + 2 * C;
}

Data literal_data(const Tournament& t) {
  Data d{t,
         0,
         Matrix(static_cast<std::size_t>(t.n),
                std::vector<i64>(static_cast<std::size_t>(t.n), 0)),
         std::vector<std::array<i64, 2>>(static_cast<std::size_t>(t.n), {0, 0}),
         std::vector<std::array<i64, 2>>(static_cast<std::size_t>(t.n), {0, 0}),
         0,
         0,
         std::vector<i64>(static_cast<std::size_t>(t.n + 1), 0)};
  std::vector<int> word(static_cast<std::size_t>(t.n));
  std::iota(word.begin(), word.end(), 0);
  do {
    std::vector<bool> good(static_cast<std::size_t>(std::max(0, t.n - 1)), false);
    int bad = 0;
    for (int gap = 0; gap + 1 < t.n; ++gap) {
      good[static_cast<std::size_t>(gap)] = arc(t, word[static_cast<std::size_t>(gap)],
                                                word[static_cast<std::size_t>(gap + 1)]);
      bad += !good[static_cast<std::size_t>(gap)];
    }
    if (bad == 0) {
      ++d.h;
      ++d.starts[static_cast<std::size_t>(word.front())][static_cast<std::size_t>(t.n & 1)];
      ++d.ends[static_cast<std::size_t>(word.back())][static_cast<std::size_t>(t.n & 1)];
    }
    if (bad <= 1) {
      for (int gap = 0; gap + 1 < t.n; ++gap) {
        if (bad == 0 || !good[static_cast<std::size_t>(gap)]) {
          const int u = word[static_cast<std::size_t>(gap)];
          const int v = word[static_cast<std::size_t>(gap + 1)];
          ++d.c[static_cast<std::size_t>(u)][static_cast<std::size_t>(v)];
          ++d.c[static_cast<std::size_t>(v)][static_cast<std::size_t>(u)];
        }
      }
    }
    for (int cut = 1; cut < t.n; ++cut) {
      bool left_good = true;
      bool right_good = true;
      for (int gap = 0; gap < cut - 1; ++gap) {
        left_good = left_good && good[static_cast<std::size_t>(gap)];
      }
      for (int gap = cut; gap + 1 < t.n; ++gap) {
        right_good = right_good && good[static_cast<std::size_t>(gap)];
      }
      if (left_good && right_good) {
        ++d.starts[static_cast<std::size_t>(word.front())][static_cast<std::size_t>(cut & 1)];
        ++d.ends[static_cast<std::size_t>(word[static_cast<std::size_t>(cut - 1)])]
                [static_cast<std::size_t>(cut & 1)];
        if (t.n % 2 == 0 && cut % 2 == 0) {
          ++d.even_component_objects[static_cast<std::size_t>(cut)];
          need(word.front() != word[static_cast<std::size_t>(cut)]
                   && word.front() != word.back()
                   && word[static_cast<std::size_t>(cut - 1)]
                          != word[static_cast<std::size_t>(cut)]
                   && word[static_cast<std::size_t>(cut - 1)] != word.back(),
               "component endpoint disjointness");
        }
      }
    }
  } while (std::next_permutation(word.begin(), word.end()));
  for (int i = 0; i < t.n; ++i) {
    for (int j = i + 1; j < t.n; ++j) {
      d.mass += d.c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    }
  }
  d.gate = gate(t, d.c);
  return d;
}

Tournament ordinal(const Tournament& a, const Tournament& b) {
  Tournament t{a.n + b.n,
               std::vector<std::uint32_t>(static_cast<std::size_t>(a.n + b.n), 0)};
  const std::uint32_t bmask = ((1U << b.n) - 1U) << a.n;
  for (int i = 0; i < a.n; ++i) {
    t.out[static_cast<std::size_t>(i)] = a.out[static_cast<std::size_t>(i)] | bmask;
  }
  for (int i = 0; i < b.n; ++i) {
    t.out[static_cast<std::size_t>(a.n + i)] = b.out[static_cast<std::size_t>(i)] << a.n;
  }
  return t;
}

struct Packet {
  Matrix c;
  i64 total{};
  i64 left{};
  i64 right{};
  i64 remainder{};
};

Packet packet(const Data& a, const Data& b) {
  const int n = a.t.n + b.t.n;
  Matrix c(static_cast<std::size_t>(n), std::vector<i64>(static_cast<std::size_t>(n), 0));
  for (int i = 0; i < a.t.n; ++i) {
    for (int j = i + 1; j < a.t.n; ++j) {
      c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]
          = c[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)]
          = b.h * a.c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    }
  }
  for (int i = 0; i < b.t.n; ++i) {
    for (int j = i + 1; j < b.t.n; ++j) {
      c[static_cast<std::size_t>(a.t.n + i)][static_cast<std::size_t>(a.t.n + j)]
          = c[static_cast<std::size_t>(a.t.n + j)][static_cast<std::size_t>(a.t.n + i)]
          = a.h * b.c[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    }
  }
  for (int i = 0; i < a.t.n; ++i) {
    for (int j = 0; j < b.t.n; ++j) {
      const i64 value = 2 * (a.starts[static_cast<std::size_t>(i)][0]
                                  * b.ends[static_cast<std::size_t>(j)][0]
                              + a.starts[static_cast<std::size_t>(i)][1]
                                    * b.ends[static_cast<std::size_t>(j)][1]);
      c[static_cast<std::size_t>(i)][static_cast<std::size_t>(a.t.n + j)] = value;
      c[static_cast<std::size_t>(a.t.n + j)][static_cast<std::size_t>(i)] = value;
    }
  }
  const Tournament t = ordinal(a.t, b.t);
  const i64 total = gate(t, c);
  const i64 left = b.h * b.h * a.gate;
  const i64 right = a.h * a.h * b.gate;
  return Packet{std::move(c), total, left, right, total - left - right};
}

std::vector<Data> labelled_bank(int n) {
  const int pairs = n * (n - 1) / 2;
  const std::uint64_t count = 1ULL << pairs;
  std::vector<Data> bank;
  bank.reserve(static_cast<std::size_t>(count));
  for (std::uint64_t mask = 0; mask < count; ++mask) bank.push_back(literal_data(from_mask(n, mask)));
  return bank;
}

}  // namespace

int main() {
  try {
    std::vector<std::vector<Data>> banks(7);
    banks[1] = labelled_bank(1);
    banks[2] = labelled_bank(2);
    banks[3] = labelled_bank(3);
    banks[4] = labelled_bank(4);
    banks[5] = labelled_bank(5);
    banks[6] = labelled_bank(6);

    i64 sidecar_checks = 0;
    i64 exchange_checks = 0;
    for (int n = 1; n <= 6; ++n) {
      for (const Data& d : banks[static_cast<std::size_t>(n)]) {
        std::array<i64, 2> s{0, 0};
        std::array<i64, 2> e{0, 0};
        for (int v = 0; v < n; ++v) {
          for (int parity = 0; parity < 2; ++parity) {
            s[static_cast<std::size_t>(parity)]
                += d.starts[static_cast<std::size_t>(v)][static_cast<std::size_t>(parity)];
            e[static_cast<std::size_t>(parity)]
                += d.ends[static_cast<std::size_t>(v)][static_cast<std::size_t>(parity)];
          }
        }
        need(s[0] + s[1] == d.mass + d.h, "literal start-cover identity");
        need(e[0] + e[1] == d.mass + d.h, "literal end-cover identity");
        if (n % 2 == 1) {
          need(s[1] - s[0] == d.h, "literal odd start exchange");
          need(e[1] - e[0] == d.h, "literal odd end exchange");
        } else {
          for (int cut = 2; cut <= n - 2; cut += 2) {
            need(d.even_component_objects[static_cast<std::size_t>(cut)]
                     == d.even_component_objects[static_cast<std::size_t>(n - cut)],
                 "literal even component exchange");
            ++exchange_checks;
          }
        }
        ++sidecar_checks;
      }
    }

    i64 transfer_checks = 0;
    i64 remainder_checks = 0;
    i64 remainder_min = std::numeric_limits<i64>::max();
    for (int na = 1; na <= 4; ++na) {
      for (int nb = 3; nb <= 5; ++nb) {
        for (const Data& a : banks[static_cast<std::size_t>(na)]) {
          for (const Data& b : banks[static_cast<std::size_t>(nb)]) {
            if (has_sink(b.t)) continue;
            const Packet p = packet(a, b);
            need(p.remainder > 0, "labelled ordinal remainder control");
            remainder_min = std::min(remainder_min, p.remainder);
            ++remainder_checks;
            if (na <= 3 && nb <= 4) {
              const Data actual = literal_data(ordinal(a.t, b.t));
              need(actual.c == p.c, "literal direct ordinal transfer");
              ++transfer_checks;
            }
          }
        }
      }
    }

    using Dom = std::tuple<i64, std::string, std::string, int, int, int, i64, i64>;
    using FirstKey = std::tuple<int, int, int, std::string, std::string, int, int, int>;
    using First = std::pair<FirstKey, Dom>;
    std::optional<First> first_left;
    std::optional<First> first_right;
    std::optional<First> first_right_a2;
    std::optional<First> first_right_a3;
    std::map<std::tuple<int, int, int>, Dom> left_sharp;
    std::map<std::tuple<int, int, int>, Dom> right_sharp;
    i64 dominance_presentations = 0;
    for (int na = 1; na <= 4; ++na) {
      for (int nb = 3; nb <= 6; ++nb) {
        if (na + nb > 8) continue;
        for (const Data& a : banks[static_cast<std::size_t>(na)]) {
          for (const Data& b : banks[static_cast<std::size_t>(nb)]) {
            if (has_sink(b.t)) continue;
            const Packet p = packet(a, b);
            const std::string la = label(a.t);
            const std::string lb = label(b.t);
            const auto block = std::make_tuple(na + nb, na, nb);
            for (int u = 0; u < na; ++u) {
              for (int v = 0; v < na; ++v) {
                if (!arc(a.t, u, v)) continue;
                for (int w = 0; w < nb; ++w) {
                  const i64 lhs = p.c[static_cast<std::size_t>(u)][static_cast<std::size_t>(v)];
                  const i64 rhs = p.c[static_cast<std::size_t>(v)]
                                     [static_cast<std::size_t>(na + w)];
                  const Dom row{lhs - rhs, la, lb, u, v, w, lhs, rhs};
                  if (!left_sharp.contains(block) || row < left_sharp[block]) left_sharp[block] = row;
                  if (lhs < rhs) {
                    const First candidate{FirstKey{na + nb, na, nb, la, lb, u, v, w}, row};
                    if (!first_left || candidate < *first_left) first_left = candidate;
                  }
                }
              }
            }
            for (int u = 0; u < nb; ++u) {
              for (int v = 0; v < nb; ++v) {
                if (!arc(b.t, u, v)) continue;
                for (int w = 0; w < na; ++w) {
                  const i64 lhs = p.c[static_cast<std::size_t>(w)]
                                     [static_cast<std::size_t>(na + u)];
                  const i64 rhs = p.c[static_cast<std::size_t>(na + u)]
                                     [static_cast<std::size_t>(na + v)];
                  const Dom row{lhs - rhs, la, lb, w, u, v, lhs, rhs};
                  if (!right_sharp.contains(block) || row < right_sharp[block]) right_sharp[block] = row;
                  if (lhs < rhs) {
                    const First candidate{FirstKey{na + nb, na, nb, la, lb, w, u, v}, row};
                    if (!first_right || candidate < *first_right) first_right = candidate;
                    if (na >= 2 && (!first_right_a2 || candidate < *first_right_a2)) {
                      first_right_a2 = candidate;
                    }
                    if (na >= 3 && (!first_right_a3 || candidate < *first_right_a3)) {
                      first_right_a3 = candidate;
                    }
                  }
                }
              }
            }
            ++dominance_presentations;
          }
        }
      }
    }
    need(first_left && std::get<0>(first_left->first) == 6
                    && std::get<1>(first_left->first) == 3
                    && std::get<2>(first_left->first) == 3,
         "derived first left obstruction order");
    need(first_right && std::get<0>(first_right->first) == 7
                     && std::get<1>(first_right->first) == 1
                     && std::get<2>(first_right->first) == 6,
         "derived first right obstruction order");
    need(first_right_a2 && std::get<0>(first_right_a2->first) == 8
                        && std::get<1>(first_right_a2->first) == 2,
         "derived first right obstruction with A order at least two");
    need(first_right_a3 && std::get<0>(first_right_a3->first) == 8
                        && std::get<1>(first_right_a3->first) == 3,
         "derived first right obstruction with A order at least three");
    const Dom first_left_block_sharp = left_sharp.at(std::make_tuple(6, 3, 3));
    const Dom first_right_block_sharp = right_sharp.at(std::make_tuple(7, 1, 6));
    need(std::get<0>(first_left_block_sharp) == -4, "derived left first-block sharp");
    need(std::get<0>(first_right_block_sharp) == -8, "derived right first-block sharp");

    const Data c3 = literal_data(parse("101"));
    const Data singleton = literal_data(from_mask(1, 0));
    const Data primary_right_hostile = literal_data(parse("111111100111111"));
    const Packet primary_right_packet = packet(singleton, primary_right_hostile);
    need(primary_right_packet.c[0][1] == 22 && primary_right_packet.c[1][6] == 28,
         "primary-canonical first right obstruction");
    const Data primary_right_sharp_hostile = literal_data(parse("111111101111110"));
    const Packet primary_right_sharp_packet = packet(singleton, primary_right_sharp_hostile);
    need(primary_right_sharp_packet.c[0][1] == 18 && primary_right_sharp_packet.c[1][6] == 26,
         "primary-canonical first-block right sharp obstruction");
    const Data right_hostile = literal_data(parse("1111100111"));
    const Packet right_packet = packet(c3, right_hostile);
    need(right_packet.c[0][3] == 20 && right_packet.c[3][7] == 30,
         "right coordinate obstruction");
    const Packet left_packet = packet(c3, c3);
    need(left_packet.c[0][1] == 6 && left_packet.c[1][3] == 10,
         "left coordinate obstruction");
    const Data sharp_a = literal_data(parse("1100111111"));
    const Packet sharp_left = packet(sharp_a, sharp_a);
    need(sharp_left.c[0][1] == 88 && sharp_left.c[1][9] == 730,
         "sharp left stress control");
    const Data sharp_right_a = literal_data(parse("1110101111"));
    const Packet sharp_right = packet(sharp_right_a, right_hostile);
    need(sharp_right.c[4][5] == 90 && sharp_right.c[5][9] == 150,
         "sharp right stress control");

    const auto print_first = [](const char* name, const First& row) {
      const auto& [total, na, nb, la, lb, x, y, z] = row.first;
      const auto& [slack, ignored_a, ignored_b, ignored_x, ignored_y, ignored_z, lhs, rhs] = row.second;
      static_cast<void>(ignored_a);
      static_cast<void>(ignored_b);
      static_cast<void>(ignored_x);
      static_cast<void>(ignored_y);
      static_cast<void>(ignored_z);
      std::cout << name << " total=" << total << " A=" << na << ':' << la
                << " B=" << nb << ':' << lb << " coords=" << x << ',' << y << ',' << z
                << " lhs=" << lhs << " rhs=" << rhs << " slack=" << slack << '\n';
    };
    const auto print_dom = [](const char* name, const Dom& row) {
      const auto& [slack, la, lb, x, y, z, lhs, rhs] = row;
      std::cout << name << " A=" << la << " B=" << lb << " coords=" << x << ',' << y << ',' << z
                << " lhs=" << lhs << " rhs=" << rhs << " slack=" << slack << '\n';
    };

    std::cout << "THM4181_INDEPENDENT_LITERAL_AUDIT\n";
    std::cout << "labelled_counts q1=1 q2=2 q3=8 q4=64 q5=1024 q6=32768\n";
    std::cout << "literal_sidecar_checks " << sidecar_checks << '\n';
    std::cout << "literal_even_component_exchange_checks " << exchange_checks << '\n';
    std::cout << "literal_direct_transfer_checks " << transfer_checks << '\n';
    std::cout << "labelled_remainder_checks " << remainder_checks << '\n';
    std::cout << "labelled_remainder_minimum " << remainder_min << '\n';
    std::cout << "labelled_dominance_presentations " << dominance_presentations << '\n';
    print_first("derived_first_left", *first_left);
    print_first("derived_first_right", *first_right);
    print_first("derived_first_right_A_ge_2", *first_right_a2);
    print_first("derived_first_right_A_ge_3", *first_right_a3);
    print_dom("derived_first_left_block_sharp", first_left_block_sharp);
    print_dom("derived_first_right_block_sharp", first_right_block_sharp);
    std::cout << "primary_canonical_first_right B=111111100111111 lhs=22 rhs=28 slack=-6\n";
    std::cout << "primary_canonical_first_right_block_sharp B=111111101111110 lhs=18 rhs=26 slack=-8\n";
    std::cout << "left_obstruction C3 C3 values 6 10 slack -4\n";
    std::cout << "right_A_ge_3_stress C3 1111100111 values 20 30 slack -10\n";
    std::cout << "sharp_left_5x5 values 88 730 slack -642\n";
    std::cout << "sharp_right_5x5 values 90 150 slack -60\n";
    std::cout << "PASS\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL " << error.what() << '\n';
    return 1;
  }
}
