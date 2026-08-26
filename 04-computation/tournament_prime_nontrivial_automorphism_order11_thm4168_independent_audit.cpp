#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {
constexpr int N = 11;
constexpr int CHILD = 12;
constexpr int CUTS = 1 << N;
constexpr int STATES = 1 << CHILD;
using Parent = std::array<uint16_t, N>;
using Child = std::array<uint16_t, CHILD>;

Parent parse(const std::string& s) {
  if (s.size() != 55) throw std::runtime_error("bad label length");
  Parent out{};
  int k = 0;
  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j, ++k) {
      if (s[k] == '1') out[i] |= uint16_t(1u << j);
      else if (s[k] == '0') out[j] |= uint16_t(1u << i);
      else throw std::runtime_error("bad label character");
    }
  }
  return out;
}

uint64_t child_h(const Parent& parent, int cut) {
  Child out{};
  for (int i = 0; i < N; ++i) out[i] = parent[i];
  for (int v = 0; v < N; ++v) {
    if (cut & (1 << v)) out[N] |= uint16_t(1u << v);
    else out[v] |= uint16_t(1u << N);
  }

  static uint64_t dp[STATES][CHILD];
  std::fill(&dp[0][0], &dp[0][0] + STATES * CHILD, 0);
  for (int v = 0; v < CHILD; ++v) dp[1 << v][v] = 1;
  for (int mask = 1; mask < STATES; ++mask) {
    unsigned vertices = unsigned(mask);
    while (vertices) {
      const int v = __builtin_ctz(vertices);
      vertices &= vertices - 1;
      const int rest = mask ^ (1 << v);
      unsigned predecessors = unsigned(rest);
      while (predecessors) {
        const int u = __builtin_ctz(predecessors);
        predecessors &= predecessors - 1;
        if (out[u] & (1u << v)) dp[mask][v] += dp[rest][u];
      }
    }
  }
  uint64_t h = 0;
  for (int v = 0; v < CHILD; ++v) h += dp[STATES - 1][v];
  return h;
}

int64_t ceil_div(__int128 a, __int128 b) {
  const __int128 q = a / b;
  const __int128 r = a % b;
  return int64_t(q + (r > 0));
}

void audit(const std::string& label) {
  const Parent out = parse(label);
  std::array<int64_t, CUTS> response{};
  for (int cut = 0; cut < CUTS; ++cut) {
    response[cut] = int64_t(child_h(out, cut));
  }

  const int64_t H = response[0];
  int64_t capacity[N][N]{};
  int64_t W2 = 0;
  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      capacity[i][j] = capacity[j][i] =
          response[1 << i] + response[1 << j] - H
          - response[(1 << i) | (1 << j)];
      if (capacity[i][j] <= 0 || (capacity[i][j] & 1)) {
        throw std::runtime_error("bad recovered capacity");
      }
      W2 += capacity[i][j];
    }
  }

  int64_t d2[N]{}, h2[N]{};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (i == j) continue;
      d2[i] += capacity[i][j];
      h2[i] += (out[i] & (1u << j)) ? capacity[i][j] : -capacity[i][j];
    }
  }

  int64_t Chdx4 = 0;
  int64_t D4x4 = 0;
  for (int i = 0; i < N; ++i) Chdx4 += d2[i] * h2[i];
  for (int a = 0; a < N; ++a) {
    for (int b = a + 1; b < N; ++b) {
      for (int c = a + 1; c < N; ++c) {
        if (c == b) continue;
        for (int d = c + 1; d < N; ++d) {
          if (d != b) D4x4 += capacity[a][b] * capacity[c][d];
        }
      }
    }
  }

  int64_t central = std::numeric_limits<int64_t>::min();
  int64_t outer = central;
  uint64_t digest = 14695981039346656037ull;
  for (int cut = 0; cut < CUTS; ++cut) {
    const uint64_t value = uint64_t(response[cut]);
    for (int k = 0; k < 8; ++k) {
      digest ^= (value >> (8 * k)) & 255u;
      digest *= 1099511628211ull;
    }
  }
  for (int m = 1; m < N; ++m) {
    int64_t count = 0, sum = 0, square_sum = 0, anchor = -1, lattice = 0;
    for (int cut = 0; cut < CUTS; ++cut) {
      if (__builtin_popcount(unsigned(cut)) != m) continue;
      const int64_t value = response[cut];
      if (anchor < 0) anchor = value;
      ++count;
      sum += value;
      square_sum += value * value;
      lattice = std::gcd(lattice, std::llabs(value - anchor));
    }
    const int64_t delta = sum - count * H;
    const __int128 numerator = __int128(square_sum) - __int128(sum) * H;
    int64_t floor = anchor;
    if (lattice) {
      floor = anchor + lattice * ceil_div(
          numerator - __int128(anchor) * delta,
          __int128(delta) * lattice);
    }
    if (m == 5 || m == 6) central = std::max(central, floor);
    else outer = std::max(outer, floor);
  }

  const int64_t numerator = 2 * std::llabs(Chdx4);
  const int64_t denominator = D4x4;
  const int64_t divisor = std::gcd(numerator, denominator);
  std::cout << "label=" << label << " H=" << H << " W2=" << W2
            << " D4x4=" << D4x4 << " Chdx4=" << Chdx4
            << " exact_rho=" << numerator / divisor << '/' << denominator / divisor
            << " margin=" << central - outer
            << " response_fnv64=" << std::hex << digest << std::dec << '\n';
}
}  // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
  _setmode(_fileno(stdout), _O_BINARY);
#endif
  try {
    if (argc < 2) {
      std::cerr << "usage: labels...\n";
      return 2;
    }
    for (int i = 1; i < argc; ++i) audit(argv[i]);
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
