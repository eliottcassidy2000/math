// Independent exact verifier for the proposed THM-4026 counterexample.
//
// Sun's domain is w,x,y,z >= 2 with C(n,k)=0 for 0<=n<k.  For existence,
// the repeated zero values may be represented uniquely by x=3,y=5,z=7.
// For fixed x,y,z put R=N-C(x,4)-C(y,6)-C(z,8).  Then w>=2 exists iff
// R>=1 and 8R+1 is an odd square at least 9, in which case
// w=(1+sqrt(8R+1))/2.
//
// Two identities give a useful norm formulation:
//
//   8*C(w,2)+1 = (2w-1)^2,
//   24*C(x,4)+1 = (x^2-3x+1)^2.
//
// Writing A=x^2-3x+1, B=2w-1 and S=C(y,6)+C(z,8), multiplication of
// N-S=C(w,2)+C(x,4) by 24 gives
//
//   A^2+3B^2 = 24(N-S)+4.
//
// Thus the residual problem is a Q(sqrt(-3)) field-norm equation with the
// extra thin constraint that A belongs to the quadratic image x^2-3x+1.
//
// Reproduction commands are frozen in
// 05-knowledge/results/sun_2468_counterexample_thm4026_independent.out.

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 choose_small_k(u64 n, unsigned k) {
    if (n < k) return 0;
    u128 ans = 1;
    for (unsigned j = 1; j <= k; ++j) {
        ans = ans * (n - k + j) / j;
    }
    if (ans > std::numeric_limits<u64>::max()) {
        std::cerr << "binomial overflow\n";
        std::exit(2);
    }
    return static_cast<u64>(ans);
}

static u64 max_upper(u64 target, unsigned k, u64 lo) {
    u64 hi = std::max<u64>(lo + 1, 2);
    while (choose_small_k(hi, k) <= target) {
        if (hi > 1000000000ULL) {
            std::cerr << "range search unexpectedly large\n";
            std::exit(2);
        }
        hi *= 2;
    }
    u64 left = lo;
    while (left + 1 < hi) {
        u64 mid = left + (hi - left) / 2;
        if (choose_small_k(mid, k) <= target) left = mid;
        else hi = mid;
    }
    return left;
}

static u64 isqrt_u64(u64 n) {
    if (n == 0) return 0;
    u64 x = u64(1) << ((std::bit_width(n) + 1U) / 2U);
    for (;;) {
        const u64 y = static_cast<u64>((static_cast<u128>(x) + n / x) / 2);
        if (y >= x) return x;
        x = y;
    }
}

struct Filter {
    unsigned p;
    std::vector<u64> table;
};

struct Result {
    u64 yz_pairs = 0;
    u64 raw_triples = 0;
    u64 sieve_survivors = 0;
    u64 nonnegative_survivors = 0;
    u64 square_hits = 0;
    std::vector<std::tuple<u64,u64,u64,u64>> witnesses;
};

static Result verify(u64 target, const std::vector<unsigned>& primes) {
    const u64 xmin = 3, ymin = 5, zmin = 7;
    const u64 wmax = max_upper(target, 2, 2);
    const u64 xmax = max_upper(target, 4, xmin);
    const u64 ymax = max_upper(target, 6, ymin);
    const u64 zmax = max_upper(target, 8, zmin);

    std::vector<u64> c4, c6, c8;
    c4.reserve(xmax - xmin + 1);
    c6.reserve(ymax - ymin + 1);
    c8.reserve(zmax - zmin + 1);
    for (u64 x = xmin; x <= xmax; ++x) c4.push_back(choose_small_k(x, 4));
    for (u64 y = ymin; y <= ymax; ++y) c6.push_back(choose_small_k(y, 6));
    for (u64 z = zmin; z <= zmax; ++z) c8.push_back(choose_small_k(z, 8));

    const std::size_t nx = c4.size();
    const std::size_t words = (nx + 63) / 64;
    std::vector<Filter> filters;
    filters.reserve(primes.size());

    for (unsigned p : primes) {
        std::vector<unsigned char> qr(p, 0);
        for (unsigned a = 0; a < p; ++a) qr[(u64(a) * a) % p] = 1;
        Filter f{p, std::vector<u64>(std::size_t(p) * words, 0)};
        const u64 tp = target % p;
        for (unsigned r = 0; r < p; ++r) {
            u64* row = f.table.data() + std::size_t(r) * words;
            for (std::size_t ix = 0; ix < nx; ++ix) {
                u64 q = (tp + p - r) % p;
                q = (q + p - (c4[ix] % p)) % p;
                const unsigned d = unsigned((8 * q + 1) % p);
                if (qr[d]) row[ix >> 6] |= u64(1) << (ix & 63);
            }
        }
        filters.push_back(std::move(f));
    }

    Result total;

#pragma omp parallel
    {
        Result local;
        std::vector<u64> bits(words);
#pragma omp for schedule(dynamic, 1)
        for (std::int64_t iy_signed = 0; iy_signed < static_cast<std::int64_t>(c6.size()); ++iy_signed) {
            const std::size_t iy = static_cast<std::size_t>(iy_signed);
            const u64 vy = c6[iy];
            for (std::size_t iz = 0; iz < c8.size(); ++iz) {
                const u64 vz = c8[iz];
                if (vy > target || vz > target - vy) break;
                const u64 s = vy + vz;
                ++local.yz_pairs;
                local.raw_triples += static_cast<u64>(
                    std::upper_bound(c4.begin(), c4.end(), target - s) - c4.begin());

                const Filter& f0 = filters[0];
                const u64* row0 = f0.table.data() + std::size_t(s % f0.p) * words;
                std::copy(row0, row0 + words, bits.begin());
                for (std::size_t j = 1; j < filters.size(); ++j) {
                    const Filter& f = filters[j];
                    const u64* row = f.table.data() + std::size_t(s % f.p) * words;
                    for (std::size_t q = 0; q < words; ++q) bits[q] &= row[q];
                }

                for (std::size_t q = 0; q < words; ++q) {
                    u64 word = bits[q];
                    local.sieve_survivors += static_cast<u64>(__builtin_popcountll(word));
                    while (word) {
                        const unsigned b = static_cast<unsigned>(__builtin_ctzll(word));
                        const std::size_t ix = (q << 6) + b;
                        word &= word - 1;
                        if (ix >= nx) continue;
                        if (c4[ix] > target - s) continue;
                        ++local.nonnegative_survivors;
                        const u64 rem = target - s - c4[ix];
                        const u64 disc = 8 * rem + 1;
                        const u64 root = isqrt_u64(disc);
                        if (static_cast<u128>(root) * root == disc && root >= 3 && (root & 1)) {
                            ++local.square_hits;
                            if (local.witnesses.size() < 32) {
                                local.witnesses.emplace_back((root + 1) / 2,
                                    xmin + ix, ymin + iy, zmin + iz);
                            }
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            total.yz_pairs += local.yz_pairs;
            total.raw_triples += local.raw_triples;
            total.sieve_survivors += local.sieve_survivors;
            total.nonnegative_survivors += local.nonnegative_survivors;
            total.square_hits += local.square_hits;
            total.witnesses.insert(total.witnesses.end(), local.witnesses.begin(), local.witnesses.end());
        }
    }

    std::sort(total.witnesses.begin(), total.witnesses.end());
    std::cout << "target=" << target << "\n";
    std::cout << "canonical_ranges=w=2.." << wmax << " x=" << xmin << ".." << xmax
              << " y=" << ymin << ".." << ymax << " z=" << zmin << ".." << zmax << "\n";
    std::cout << "boundary_next=C(" << (wmax+1) << ",2)=" << choose_small_k(wmax+1,2)
              << " C(" << (xmax+1) << ",4)=" << choose_small_k(xmax+1,4)
              << " C(" << (ymax+1) << ",6)=" << choose_small_k(ymax+1,6)
              << " C(" << (zmax+1) << ",8)=" << choose_small_k(zmax+1,8) << "\n";
    std::cout << "primes=";
    for (std::size_t i = 0; i < primes.size(); ++i) std::cout << (i ? "," : "") << primes[i];
    std::cout << "\n";
    std::cout << "yz_pairs=" << total.yz_pairs << " raw_xyz_triples=" << total.raw_triples
              << " sieve_survivors=" << total.sieve_survivors
              << " nonnegative_survivors=" << total.nonnegative_survivors
              << " representations=" << total.square_hits << "\n";
    for (const auto& [w,x,y,z] : total.witnesses) {
        const u64 sum = choose_small_k(w,2) + choose_small_k(x,4)
                      + choose_small_k(y,6) + choose_small_k(z,8);
        std::cout << "witness=" << w << "," << x << "," << y << "," << z
                  << " sum=" << sum << "\n";
    }
    return total;
}

int main(int argc, char** argv) {
    const u64 target = argc > 1 ? std::stoull(argv[1]) : 896315812331399ULL;
    const int bank = argc > 2 ? std::stoi(argv[2]) : 1;
    const std::vector<unsigned> bank1{3,5,7,11,13,17,19,23,29,31,37,41,43,47,53};
    const std::vector<unsigned> bank2{59,61,67,71,73,79,83,89,97,101,103,107,109,113,127};
    std::vector<unsigned> bank3 = bank1;
    bank3.insert(bank3.end(), bank2.begin(), bank2.end());
    const std::vector<unsigned> tail{131,137,139,149,151,157,163};
    std::vector<unsigned> bank4 = bank3;
    bank4.insert(bank4.end(), tail.begin(), tail.end());
    std::vector<unsigned> selected;
    if (bank >= 41 && bank <= 47) {
        selected = bank3;
        selected.insert(selected.end(), tail.begin(), tail.begin() + (bank - 40));
    } else {
        selected = bank == 4 ? bank4 : (bank == 3 ? bank3 : (bank == 2 ? bank2 : bank1));
    }
    if (argc > 3) {
        selected.clear();
        std::string spec = argv[3];
        std::size_t start = 0;
        while (start < spec.size()) {
            const std::size_t comma = spec.find(',', start);
            selected.push_back(static_cast<unsigned>(std::stoul(spec.substr(start, comma - start))));
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
    }
    const Result result = verify(target, selected);
    return result.square_hits == 0 ? 0 : 1;
}
