#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using u64 = std::uint64_t;
using u128 = unsigned __int128;

namespace {

// Explicit proof gates remain active even in builds that define NDEBUG.
void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "CHECK FAILED: " << message << '\n';
        std::abort();
    }
}

constexpr u64 N = 896315812331399ULL;
constexpr std::uint32_t MOD1 = 3U * 5U * 7U * 11U * 13U;       // 15015
constexpr std::uint32_t MOD2 = 17U * 19U * 23U * 29U;          // 215441
constexpr u64 FNV_OFFSET = 14695981039346656037ULL;
constexpr u64 FNV_PRIME = 1099511628211ULL;

u64 choose_u64(u64 n, unsigned r) {
    if (n < r) {
        return 0;
    }
    r = static_cast<unsigned>(std::min<u64>(r, n - r));
    u128 value = 1;
    for (unsigned i = 1; i <= r; ++i) {
        value = value * (n - r + i) / i;
    }
    require(value <= std::numeric_limits<u64>::max(), "binomial fits u64");
    return static_cast<u64>(value);
}

u64 largest_k(unsigned r, u64 first, u64 limit) {
    require(choose_u64(first, r) <= limit, "first binomial is feasible");
    u64 lo = first;
    u64 hi = std::max<u64>(first + 1, 2);
    while (choose_u64(hi, r) <= limit) {
        lo = hi;
        hi *= 2;
    }
    while (lo + 1 < hi) {
        const u64 mid = lo + (hi - lo) / 2;
        if (choose_u64(mid, r) <= limit) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

std::vector<u64> binomial_values(unsigned r, u64 first, u64 limit) {
    const u64 last = largest_k(r, first, limit);
    std::vector<u64> values;
    values.reserve(static_cast<std::size_t>(last - first + 1));
    for (u64 k = first; k <= last; ++k) {
        values.push_back(choose_u64(k, r));
    }
    return values;
}

// Restoring binary square root.  This returns floor(sqrt(n)) using integer
// operations only; it is independent of libm and floating-point rounding.
u64 isqrt_u64(u64 n) {
    u64 remainder = n;
    u64 root = 0;
    u64 bit = u64{1} << 62;  // greatest power of four representable in u64
    while (bit > remainder) {
        bit >>= 2;
    }
    while (bit != 0) {
        if (remainder >= root + bit) {
            remainder -= root + bit;
            root = (root >> 1) + bit;
        } else {
            root >>= 1;
        }
        bit >>= 2;
    }
    return root;
}

std::vector<std::uint8_t> triangular_residues(std::uint32_t modulus) {
    require(modulus % 2 == 1, "triangular residue modulus is odd");
    std::vector<std::uint8_t> table(modulus, 0);
    // For odd M, C(k+M,2)-C(k,2) is divisible by M, so k=0,...,M-1
    // covers every possible triangular residue modulo M.
    for (u64 k = 0; k < modulus; ++k) {
        const u64 t = k == 0 ? 0 : k * (k - 1) / 2;
        table[static_cast<std::size_t>(t % modulus)] = 1;
    }
    return table;
}

u64 mix64(u64 x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

u64 fnv_word(u64 hash, u64 word) {
    for (unsigned shift = 0; shift < 64; shift += 8) {
        hash ^= (word >> shift) & 0xffU;
        hash *= FNV_PRIME;
    }
    return hash;
}

struct Representation {
    u64 w = 0;
    u64 x = 0;
    u64 y = 0;
    u64 z = 0;
};

struct ZStats {
    u64 z = 0;
    u64 c8 = 0;
    u64 y_rows = 0;
    u64 feasible_xyz = 0;
    u64 zero_m = 0;
    u64 sieve1 = 0;
    u64 sieve2 = 0;
    u64 square_hits = 0;
    u64 survivor_hash = FNV_OFFSET;
    std::optional<Representation> first;
};

struct Tables {
    std::vector<u64> c4;
    std::vector<u64> c6;
    std::vector<u64> c8;
    std::vector<std::uint32_t> c4_mod1;
    std::vector<std::uint32_t> c4_mod2;
    std::vector<std::uint8_t> tri1;
    std::vector<std::uint8_t> tri2;
};

Tables make_tables(u64 n) {
    Tables t;
    t.c4 = binomial_values(4, 3, n);
    t.c6 = binomial_values(6, 5, n);
    t.c8 = binomial_values(8, 7, n);
    t.c4_mod1.reserve(t.c4.size());
    t.c4_mod2.reserve(t.c4.size());
    for (const u64 value : t.c4) {
        t.c4_mod1.push_back(static_cast<std::uint32_t>(value % MOD1));
        t.c4_mod2.push_back(static_cast<std::uint32_t>(value % MOD2));
    }
    t.tri1 = triangular_residues(MOD1);
    t.tri2 = triangular_residues(MOD2);
    return t;
}

std::uint32_t subtract_residue(std::uint32_t a, std::uint32_t b,
                               std::uint32_t modulus) {
    return a >= b ? a - b : a + modulus - b;
}

std::vector<ZStats> exhaustive_scan(u64 n, const Tables& t) {
    std::vector<ZStats> all(t.c8.size());

#pragma omp parallel for schedule(dynamic, 1)
    for (std::int64_t zi_signed = 0;
         zi_signed < static_cast<std::int64_t>(t.c8.size()); ++zi_signed) {
        const std::size_t zi = static_cast<std::size_t>(zi_signed);
        ZStats st;
        st.z = 7 + zi;
        st.c8 = t.c8[zi];
        const u64 after_z = n - st.c8;

        for (std::size_t yi = 0; yi < t.c6.size() && t.c6[yi] <= after_z; ++yi) {
            ++st.y_rows;
            const u64 residual = after_z - t.c6[yi];
            const auto x_end_it = std::upper_bound(t.c4.begin(), t.c4.end(), residual);
            const std::size_t x_end = static_cast<std::size_t>(x_end_it - t.c4.begin());
            st.feasible_xyz += x_end;
            const std::uint32_t residual1 = static_cast<std::uint32_t>(residual % MOD1);
            const std::uint32_t residual2 = static_cast<std::uint32_t>(residual % MOD2);

            for (std::size_t xi = 0; xi < x_end; ++xi) {
                const u64 m = residual - t.c4[xi];
                if (m == 0) {
                    ++st.zero_m;  // w>=2 forces C(w,2)>=1
                    continue;
                }

                const std::uint32_t m1 = subtract_residue(
                    residual1, t.c4_mod1[xi], MOD1);
#ifdef DEEP_RESIDUE_AUDIT
                require(m1 == static_cast<std::uint32_t>(m % MOD1),
                        "optimized MOD1 residue agrees with division");
#endif
                if (!t.tri1[m1]) {
                    continue;
                }
                ++st.sieve1;

                const std::uint32_t m2 = subtract_residue(
                    residual2, t.c4_mod2[xi], MOD2);
#ifdef DEEP_RESIDUE_AUDIT
                require(m2 == static_cast<std::uint32_t>(m % MOD2),
                        "optimized MOD2 residue agrees with division");
#endif
                if (!t.tri2[m2]) {
                    continue;
                }
                ++st.sieve2;

                const u64 d = 8 * m + 1;
                const u64 s = isqrt_u64(d);
                require(static_cast<u128>(s) * s <= d, "isqrt lower bound");
                require(static_cast<u128>(s + 1) * (s + 1) > d,
                        "isqrt upper bound");
                const u64 token = mix64(d) ^ mix64((3 + xi) << 20)
                                  ^ mix64((5 + yi) << 10) ^ mix64(7 + zi);
                st.survivor_hash = fnv_word(st.survivor_hash, token);
                if (s * s != d) {
                    continue;
                }

                ++st.square_hits;
                require((s & 1U) == 1U, "square discriminant has odd root");
                const u64 w = (s + 1) / 2;
                const Representation rep{w, 3 + xi, 5 + yi, 7 + zi};
                require(w >= 2, "recovered triangular index is admissible");
                require(choose_u64(rep.w, 2) + choose_u64(rep.x, 4)
                            + choose_u64(rep.y, 6) + choose_u64(rep.z, 8)
                        == n,
                        "recovered witness resums to target");
                if (!st.first.has_value()) {
                    st.first = rep;
                }
            }
        }
        all[zi] = st;
    }
    return all;
}

std::optional<Representation> find_control(u64 n) {
    const Tables t = make_tables(n);
    std::atomic<bool> found{false};
    std::vector<std::optional<Representation>> by_z(t.c8.size());

#pragma omp parallel for schedule(dynamic, 1)
    for (std::int64_t zi_signed = 0;
         zi_signed < static_cast<std::int64_t>(t.c8.size()); ++zi_signed) {
        if (found.load(std::memory_order_relaxed)) {
            continue;
        }
        const std::size_t zi = static_cast<std::size_t>(zi_signed);
        const u64 after_z = n - t.c8[zi];
        for (std::size_t yi = 0;
             yi < t.c6.size() && t.c6[yi] <= after_z
             && !found.load(std::memory_order_relaxed); ++yi) {
            const u64 residual = after_z - t.c6[yi];
            const std::uint32_t residual1 = static_cast<std::uint32_t>(residual % MOD1);
            const std::uint32_t residual2 = static_cast<std::uint32_t>(residual % MOD2);
            const std::size_t x_end = static_cast<std::size_t>(
                std::upper_bound(t.c4.begin(), t.c4.end(), residual) - t.c4.begin());
            for (std::size_t xi = 0;
                 xi < x_end && !found.load(std::memory_order_relaxed); ++xi) {
                const u64 m = residual - t.c4[xi];
                if (m == 0) {
                    continue;
                }
                const auto m1 = subtract_residue(residual1, t.c4_mod1[xi], MOD1);
                const auto m2 = subtract_residue(residual2, t.c4_mod2[xi], MOD2);
                if (!t.tri1[m1] || !t.tri2[m2]) {
                    continue;
                }
                const u64 d = 8 * m + 1;
                const u64 s = isqrt_u64(d);
                require(static_cast<u128>(s) * s <= d, "control isqrt lower bound");
                require(static_cast<u128>(s + 1) * (s + 1) > d,
                        "control isqrt upper bound");
                if (s * s == d) {
                    const Representation rep{(s + 1) / 2, 3 + xi, 5 + yi, 7 + zi};
                    require(choose_u64(rep.w, 2) + choose_u64(rep.x, 4)
                                + choose_u64(rep.y, 6) + choose_u64(rep.z, 8)
                            == n,
                            "control witness resums to target");
                    by_z[zi] = rep;
                    found.store(true, std::memory_order_relaxed);
                    break;
                }
            }
        }
    }

    for (const auto& rep : by_z) {
        if (rep.has_value()) {
            return rep;
        }
    }
    return std::nullopt;
}

void print_bound(unsigned r, u64 first, u64 n) {
    const u64 last = largest_k(r, first, n);
    std::cout << "BOUND r=" << r << " first=" << first << " last=" << last
              << " count=" << (last - first + 1)
              << " value_last=" << choose_u64(last, r)
              << " value_next=" << choose_u64(last + 1, r) << '\n';
}

void print_rep(u64 n, const std::optional<Representation>& rep) {
    if (!rep.has_value()) {
        std::cout << "CONTROL n=" << n << " result=NONE\n";
        return;
    }
    const auto& r = *rep;
    std::cout << "CONTROL n=" << n << " result=REP"
              << " w=" << r.w << " x=" << r.x << " y=" << r.y << " z=" << r.z
              << " terms=" << choose_u64(r.w, 2) << ',' << choose_u64(r.x, 4)
              << ',' << choose_u64(r.y, 6) << ',' << choose_u64(r.z, 8) << '\n';
}

}  // namespace

int main() {
    static_assert(MOD1 == 15015);
    static_assert(MOD2 == 215441);
    require(choose_u64(2, 2) == 1, "triangular lower-domain convention");
    require(choose_u64(3, 4) == 0, "quartic zero convention");
    require(choose_u64(5, 6) == 0, "sextic zero convention");
    require(choose_u64(7, 8) == 0, "octic zero convention");
    require(isqrt_u64(0) == 0 && isqrt_u64(1) == 1, "small isqrt controls I");
    require(isqrt_u64(2) == 1 && isqrt_u64(15) == 3 && isqrt_u64(16) == 4,
            "small isqrt controls II");

    std::cout << "CERTIFICATE sun_2468_exact_integer_sieve_v1\n";
    std::cout << "N " << N << '\n';
#ifdef _OPENMP
    std::cout << "OPENMP threads=" << omp_get_max_threads() << '\n';
#else
    std::cout << "OPENMP disabled\n";
#endif
    std::cout << "ZERO_CONVENTIONS C(2,2)=1 C(3,4)=0 C(5,6)=0 C(7,8)=0\n";
    print_bound(2, 2, N);
    print_bound(4, 3, N);
    print_bound(6, 5, N);
    print_bound(8, 7, N);
    std::cout << "SIEVES modulus1=" << MOD1 << " modulus2=" << MOD2 << '\n';
    std::cout << "ISQRT_CHECK floor_inequalities_asserted_for_every_stage2_candidate\n";

    const auto boundary_zero = find_control(0);
    const auto boundary_one = find_control(1);
    require(!boundary_zero.has_value(), "zero is excluded by w>=2");
    require(boundary_one.has_value(), "one has the boundary representation");
    require(boundary_one->w == 2 && boundary_one->x == 3
                && boundary_one->y == 5 && boundary_one->z == 7,
            "boundary-one witness uses the intended zero conventions");
    print_rep(0, boundary_zero);
    print_rep(1, boundary_one);

    const auto start = std::chrono::steady_clock::now();
    const Tables tables = make_tables(N);
    const auto stats = exhaustive_scan(N, tables);
    const auto after_scan = std::chrono::steady_clock::now();

    u64 y_rows = 0;
    u64 feasible = 0;
    u64 zero_m = 0;
    u64 sieve1 = 0;
    u64 sieve2 = 0;
    u64 hits = 0;
    std::optional<Representation> first;
    std::cout << "Z_ROWS_BEGIN\n";
    for (const ZStats& st : stats) {
        y_rows += st.y_rows;
        feasible += st.feasible_xyz;
        zero_m += st.zero_m;
        sieve1 += st.sieve1;
        sieve2 += st.sieve2;
        hits += st.square_hits;
        if (!first.has_value() && st.first.has_value()) {
            first = st.first;
        }
        std::cout << "Z z=" << st.z << " c8=" << st.c8
                  << " y_rows=" << st.y_rows
                  << " feasible=" << st.feasible_xyz
                  << " zero_m=" << st.zero_m
                  << " sieve1=" << st.sieve1
                  << " sieve2=" << st.sieve2
                  << " hits=" << st.square_hits
                  << " fnv64=" << std::hex << std::setw(16) << std::setfill('0')
                  << st.survivor_hash << std::dec << std::setfill(' ') << '\n';
    }
    std::cout << "Z_ROWS_END\n";
    std::cout << "TOTAL y_rows=" << y_rows << " feasible_xyz=" << feasible
              << " positive_m=" << (feasible - zero_m) << " zero_m=" << zero_m
              << " sieve1=" << sieve1 << " sieve2=" << sieve2
              << " square_hits=" << hits << '\n';
    if (first.has_value()) {
        print_rep(N, first);
    } else {
        std::cout << "RESULT n=" << N << " representations=0\n";
    }

    const auto control_minus = find_control(N - 1);
    const auto control_plus = find_control(N + 1);
    print_rep(N - 1, control_minus);
    print_rep(N + 1, control_plus);
    require(control_minus.has_value() && control_plus.has_value(),
            "both neighboring positive controls are represented");
    const auto finish = std::chrono::steady_clock::now();

    const double scan_seconds = std::chrono::duration<double>(after_scan - start).count();
    const double total_seconds = std::chrono::duration<double>(finish - start).count();
    std::cout << std::fixed << std::setprecision(6)
              << "TIMING scan_seconds=" << scan_seconds
              << " total_seconds=" << total_seconds << '\n';
    return hits == 0 ? 0 : 2;
}
