#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

namespace {

// Claimed counterexample from https://gist.github.com/tadamcz/
// 0c578c8b2b3fb92fe8584bc0725187e3.  This verifier was reconstructed
// independently and uses only integer arithmetic for its decisive tests.
constexpr u64 TARGET = 896315812331399ULL;
constexpr std::array<int, 31> PRIMES = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 131, 137};

u64 choose(u64 n, int r) {
  if (n < static_cast<u64>(r)) return 0;
  u128 ans = 1;
  for (int j = 1; j <= r; ++j) {
    ans = ans * (n - static_cast<u64>(r) + static_cast<u64>(j));
    if (ans % static_cast<unsigned>(j) != 0) {
      throw std::logic_error("nonexact binomial division");
    }
    ans /= static_cast<unsigned>(j);
  }
  if (ans > std::numeric_limits<u64>::max()) {
    throw std::overflow_error("binomial does not fit u64");
  }
  return static_cast<u64>(ans);
}

std::vector<u64> binomial_values(int r, u64 first_top, u64 limit) {
  std::vector<u64> out;
  for (u64 top = first_top;; ++top) {
    const u64 value = choose(top, r);
    if (value > limit) break;
    out.push_back(value);
  }
  return out;
}

u64 max_top_with_binomial_at_most(int r, u64 first_top, u64 limit) {
  u64 lo = first_top;
  u64 hi = first_top;
  while (choose(hi, r) <= limit) hi *= 2;
  while (lo + 1 < hi) {
    const u64 mid = lo + (hi - lo) / 2;
    if (choose(mid, r) <= limit) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

u64 integer_sqrt(u64 n) {
  // Restoring square root, using no floating-point operations.
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

struct PrimeMasks {
  int p;
  std::size_t words;
  // accept[b * words + j] has bit x iff
  // 8*(b-C(x+3,4))+1 is a quadratic residue modulo p.
  std::vector<u64> accept;
};

PrimeMasks make_masks(int p, const std::vector<u64>& c4) {
  const std::size_t words = (c4.size() + 63) / 64;
  std::vector<unsigned char> square(static_cast<std::size_t>(p), 0);
  for (int s = 0; s < p; ++s) square[(static_cast<long long>(s) * s) % p] = 1;
  PrimeMasks out{p, words, std::vector<u64>(static_cast<std::size_t>(p) * words, 0)};
  for (int b = 0; b < p; ++b) {
    for (std::size_t x = 0; x < c4.size(); ++x) {
      int remainder = b - static_cast<int>(c4[x] % static_cast<u64>(p));
      if (remainder < 0) remainder += p;
      const int discriminant = (8 * remainder + 1) % p;
      if (square[static_cast<std::size_t>(discriminant)]) {
        out.accept[static_cast<std::size_t>(b) * words + x / 64] |= u64{1} << (x % 64);
      }
    }
  }
  return out;
}

struct Audit {
  u64 n;
  u64 yz_pairs = 0;
  u64 triples = 0;
  u64 modular_survivors = 0;
  u64 exact_squares = 0;
  std::vector<u64> first_empty_frequency;
  std::array<u64, 4> first_witness{};
};

Audit audit(u64 n) {
  if (n == 0) throw std::invalid_argument("positive target required");
  const auto c4 = binomial_values(4, 3, n - 1);
  const auto c6 = binomial_values(6, 5, n - 1);
  const auto c8 = binomial_values(8, 7, n - 1);
  const std::size_t words = (c4.size() + 63) / 64;

  std::vector<PrimeMasks> masks;
  masks.reserve(PRIMES.size());
  for (const int p : PRIMES) masks.push_back(make_masks(p, c4));

  Audit out;
  out.n = n;
  out.first_empty_frequency.assign(PRIMES.size(), 0);
  std::vector<u64> candidates(words, 0);

  std::cout << "target=" << n << '\n';
  std::cout << "top_bounds w<=" << max_top_with_binomial_at_most(2, 2, n)
            << " x<=" << c4.size() + 2
            << " y<=" << c6.size() + 4
            << " z<=" << c8.size() + 6 << '\n';
  std::cout << "index_counts x=" << c4.size() << " y=" << c6.size()
            << " z=" << c8.size() << '\n';

  for (std::size_t zi = 0; zi < c8.size(); ++zi) {
    for (std::size_t yi = 0; yi < c6.size(); ++yi) {
      if (c8[zi] > n - 1 - c6[yi]) break;
      const u64 base = n - c8[zi] - c6[yi];
      const std::size_t xmax = static_cast<std::size_t>(
          std::upper_bound(c4.begin(), c4.end(), base - 1) - c4.begin());
      if (xmax == 0) continue;
      ++out.yz_pairs;
      out.triples += static_cast<u64>(xmax);

      const std::size_t active_words = (xmax + 63) / 64;
      std::fill(candidates.begin(), candidates.begin() + active_words, ~u64{0});
      candidates[active_words - 1] =
          (xmax % 64 == 0) ? ~u64{0} : ((u64{1} << (xmax % 64)) - 1);

      bool empty = false;
      for (std::size_t pi = 0; pi < masks.size(); ++pi) {
        const PrimeMasks& pm = masks[pi];
        const std::size_t residue = static_cast<std::size_t>(base % static_cast<u64>(pm.p));
        const u64* accept = pm.accept.data() + residue * words;
        u64 any = 0;
        for (std::size_t j = 0; j < active_words; ++j) {
          candidates[j] &= accept[j];
          any |= candidates[j];
        }
        if (any == 0) {
          ++out.first_empty_frequency[pi];
          empty = true;
          break;
        }
      }
      if (empty) continue;

      for (std::size_t word = 0; word < active_words; ++word) {
        u64 bits = candidates[word];
        while (bits != 0) {
          const unsigned bit = static_cast<unsigned>(__builtin_ctzll(bits));
          const std::size_t xi = 64 * word + bit;
          bits &= bits - 1;
          if (xi >= xmax) continue;
          ++out.modular_survivors;
          const u64 remainder = base - c4[xi];
          const u64 discriminant = 8 * remainder + 1;
          const u64 root = integer_sqrt(discriminant);
          if (root * root == discriminant) {
            ++out.exact_squares;
            if (out.exact_squares == 1) {
              const u64 w = (root + 1) / 2;
              out.first_witness = {w, static_cast<u64>(xi + 3),
                                   static_cast<u64>(yi + 5),
                                   static_cast<u64>(zi + 7)};
            }
          }
        }
      }
    }
  }

  std::cout << "yz_pairs=" << out.yz_pairs << " triples=" << out.triples
            << " modular_survivors=" << out.modular_survivors
            << " exact_squares=" << out.exact_squares << '\n';
  std::cout << "first_empty_by_prime";
  for (std::size_t i = 0; i < PRIMES.size(); ++i) {
    if (out.first_empty_frequency[i] != 0) {
      std::cout << ' ' << PRIMES[i] << ':' << out.first_empty_frequency[i];
    }
  }
  std::cout << '\n';
  if (out.exact_squares != 0) {
    const auto [w, x, y, z] = out.first_witness;
    const u64 sum = choose(w, 2) + choose(x, 4) + choose(y, 6) + choose(z, 8);
    if (sum != n) throw std::logic_error("witness square-back failed");
    std::cout << "first_witness=" << w << ',' << x << ',' << y << ',' << z << '\n';
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc == 2) {
    audit(std::stoull(argv[1]));
    return 0;
  }
  if (argc != 1) {
    std::cerr << "usage: verifier [target]\n";
    return 2;
  }
  const Audit left = audit(TARGET - 1);
  std::cout << '\n';
  const Audit center = audit(TARGET);
  std::cout << '\n';
  const Audit right = audit(TARGET + 1);
  std::cout << '\n';
  if (center.yz_pairs != 248160 || center.triples != 2755643831ULL ||
      center.modular_survivors != 0 || center.exact_squares != 0) {
    throw std::logic_error("central certificate regression");
  }
  if (left.exact_squares != 89 ||
      left.first_witness != std::array<u64, 4>{33663667, 9433, 16, 9}) {
    throw std::logic_error("left control regression");
  }
  if (right.exact_squares != 67 ||
      right.first_witness != std::array<u64, 4>{40920205, 6138, 22, 13}) {
    throw std::logic_error("right control regression");
  }
  std::cout << "certificate=PASS\n";
  return 0;
}
