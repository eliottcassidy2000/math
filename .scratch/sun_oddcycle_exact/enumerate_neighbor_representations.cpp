#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

namespace {

constexpr std::array<int, 31> PRIMES = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 131, 137};

u64 choose(u64 n, int r) {
  if (n < static_cast<u64>(r)) return 0;
  u128 ans = 1;
  for (int j = 1; j <= r; ++j) {
    ans *= n - static_cast<u64>(r) + static_cast<u64>(j);
    if (ans % static_cast<unsigned>(j) != 0) throw std::logic_error("division");
    ans /= static_cast<unsigned>(j);
  }
  if (ans > std::numeric_limits<u64>::max()) throw std::overflow_error("choose");
  return static_cast<u64>(ans);
}

std::vector<u64> values(int r, u64 first, u64 limit) {
  std::vector<u64> out;
  for (u64 n = first;; ++n) {
    const u64 v = choose(n, r);
    if (v > limit) break;
    out.push_back(v);
  }
  return out;
}

u64 isqrt(u64 n) {
  u64 result = 0;
  u64 bit = u64{1} << 62;
  while (bit > n) bit >>= 2;
  while (bit != 0) {
    if (n >= result + bit) {
      n -= result + bit;
      result = (result >> 1) + bit;
    } else {
      result >>= 1;
    }
    bit >>= 2;
  }
  return result;
}

struct Masks {
  int p;
  std::size_t words;
  std::vector<u64> accept;
};

Masks make_masks(int p, const std::vector<u64>& c4) {
  const std::size_t words = (c4.size() + 63) / 64;
  std::vector<unsigned char> square(static_cast<std::size_t>(p));
  for (int s = 0; s < p; ++s) square[static_cast<std::size_t>((1LL * s * s) % p)] = 1;
  Masks out{p, words, std::vector<u64>(static_cast<std::size_t>(p) * words)};
  for (int b = 0; b < p; ++b) {
    for (std::size_t x = 0; x < c4.size(); ++x) {
      int rem = b - static_cast<int>(c4[x] % static_cast<u64>(p));
      if (rem < 0) rem += p;
      if (square[static_cast<std::size_t>((8 * rem + 1) % p)]) {
        out.accept[static_cast<std::size_t>(b) * words + x / 64] |= u64{1} << (x % 64);
      }
    }
  }
  return out;
}

std::vector<std::array<u64, 4>> enumerate(u64 target) {
  const auto c4 = values(4, 3, target - 1);
  const auto c6 = values(6, 5, target - 1);
  const auto c8 = values(8, 7, target - 1);
  const std::size_t words = (c4.size() + 63) / 64;
  std::vector<Masks> masks;
  for (const int p : PRIMES) masks.push_back(make_masks(p, c4));
  std::vector<u64> candidates(words);
  std::vector<std::array<u64, 4>> reps;
  u64 triples = 0;
  for (std::size_t zi = 0; zi < c8.size(); ++zi) {
    for (std::size_t yi = 0; yi < c6.size(); ++yi) {
      if (c8[zi] > target - 1 - c6[yi]) break;
      const u64 base = target - c8[zi] - c6[yi];
      const std::size_t xmax = static_cast<std::size_t>(
          std::upper_bound(c4.begin(), c4.end(), base - 1) - c4.begin());
      if (xmax == 0) continue;
      triples += static_cast<u64>(xmax);
      const std::size_t active = (xmax + 63) / 64;
      std::fill(candidates.begin(), candidates.begin() + active, ~u64{0});
      candidates[active - 1] =
          (xmax % 64 == 0) ? ~u64{0} : ((u64{1} << (xmax % 64)) - 1);
      bool empty = false;
      for (const auto& pm : masks) {
        const u64* accept = pm.accept.data() + (base % static_cast<u64>(pm.p)) * words;
        u64 any = 0;
        for (std::size_t j = 0; j < active; ++j) {
          candidates[j] &= accept[j];
          any |= candidates[j];
        }
        if (any == 0) {
          empty = true;
          break;
        }
      }
      if (empty) continue;
      for (std::size_t word = 0; word < active; ++word) {
        u64 bits = candidates[word];
        while (bits != 0) {
          const unsigned bit = static_cast<unsigned>(__builtin_ctzll(bits));
          const std::size_t xi = 64 * word + bit;
          bits &= bits - 1;
          if (xi >= xmax) continue;
          const u64 discriminant = 8 * (base - c4[xi]) + 1;
          const u64 root = isqrt(discriminant);
          if (root * root == discriminant) {
            const std::array<u64, 4> rep{
                (root + 1) / 2, static_cast<u64>(xi + 3),
                static_cast<u64>(yi + 5), static_cast<u64>(zi + 7)};
            if (choose(rep[0], 2) + choose(rep[1], 4) + choose(rep[2], 6) +
                    choose(rep[3], 8) != target) {
              throw std::logic_error("square back");
            }
            reps.push_back(rep);
          }
        }
      }
    }
  }
  std::cerr << "target=" << target << " triples=" << triples
            << " representations=" << reps.size() << '\n';
  return reps;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: enumerate_neighbor_representations TARGET OUTPUT.tsv\n";
    return 2;
  }
  const u64 target = std::stoull(argv[1]);
  const auto reps = enumerate(target);
  std::ofstream out(argv[2]);
  if (!out) throw std::runtime_error("open output");
  out << "w\tx\ty\tz\n";
  for (const auto& r : reps) out << r[0] << '\t' << r[1] << '\t' << r[2] << '\t' << r[3] << '\n';
  return 0;
}
