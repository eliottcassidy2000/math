// Clean-room FINITE-EXACT all-parity audit of the raw LRC network projections.
//
// Universe: primitive 1 <= a < b < c <= H with 3 not dividing abc.
// Each live integer carrier is enumerated directly from
//   a*C0+b*C1+c*C2=0,
//   14*|Ci| < 3*(sum of the other two speeds),
//   3 not dividing Ci.
// No short-relation classifier, slice formula, or project implementation is
// used to compute a projection.  A Diophantine row parametrization merely
// avoids scanning an empty quadratic box.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using I = std::int64_t;
using W = __int128_t;
using Triple = std::array<I, 3>;

static void need(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

static I floor_div(I n, I d) {
    need(d > 0, "floor_div denominator");
    I q = n / d, r = n % d;
    if (r < 0) --q;
    return q;
}

static I ceil_div(I n, I d) {
    return -floor_div(-n, d);
}

static I positive_mod(I x, I m) {
    I r = x % m;
    return r < 0 ? r + m : r;
}

static I inverse_mod(I a, I m) {
    if (m == 1) return 0;
    I old_r = a, r = m, old_s = 1, s = 0;
    while (r) {
        const I q = old_r / r;
        const I next_r = old_r - q * r;
        old_r = r; r = next_r;
        const I next_s = old_s - q * s;
        old_s = s; s = next_s;
    }
    need(old_r == 1, "inverse gcd");
    return positive_mod(old_s, m);
}

static std::string triple_string(const Triple& x) {
    return std::to_string(x[0]) + "," + std::to_string(x[1]) + "," + std::to_string(x[2]);
}

static std::string fraction(I p, I q) {
    const I g = std::gcd(p, q);
    return std::to_string(p / g) + "/" + std::to_string(q / g);
}

struct Row {
    I denominator = 1;
    Triple numerator{0, 0, 0};
    I positive_carriers = 0; // one from each +/- pair
};

static Row raw_row(const Triple& w) {
    const I a = w[0], b = w[1], c = w[2];
    Row row;
    row.denominator = 14 * a * b * c;
    const I m0 = (3 * (b + c) - 1) / 14;
    const I m1 = (3 * (a + c) - 1) / 14;
    const I m2 = (3 * (a + b) - 1) / 14;
    const I d = std::gcd(b, c);
    need(std::gcd(a, d) == 1, "raw_row called on imprimitive speeds");
    const I b0 = b / d, c0 = c / d;
    const I inv = inverse_mod(positive_mod(b0, c0), c0);

    // C0 is nonzero, so global sign reversal lets us take C0>0 and double.
    for (I C0 = d; C0 <= m0; C0 += d) {
        if (C0 % 3 == 0) continue;
        const I rhs = -a * C0;
        need(rhs % d == 0, "divisibility parametrization");
        const I n = rhs / d;
        const I C1_base = c0 == 1 ? 0 : positive_mod(positive_mod(n, c0) * inv, c0);
        const I C2_base = (rhs - b * C1_base) / c;
        need(b * C1_base + c * C2_base == rhs, "particular solution");

        I lo = std::max(ceil_div(-m1 - C1_base, c0),
                        ceil_div(C2_base - m2, b0));
        I hi = std::min(floor_div(m1 - C1_base, c0),
                        floor_div(C2_base + m2, b0));
        if (lo > hi) continue;

        int allowed_residue = -1;
        for (int r = 0; r < 3; ++r) {
            const I C1_mod = positive_mod(C1_base + c0 * r, 3);
            const I C2_mod = positive_mod(C2_base - b0 * r, 3);
            if (C1_mod != 0 && C2_mod != 0) {
                need(allowed_residue == -1, "multiple live t residues");
                allowed_residue = r;
            }
        }
        need(allowed_residue >= 0, "missing live t residue");
        I t = lo + positive_mod(allowed_residue - lo, 3);
        for (; t <= hi; t += 3) {
            const I C1 = C1_base + c0 * t;
            const I C2 = C2_base - b0 * t;
            need(C1 != 0 && C2 != 0 && C1 % 3 != 0 && C2 % 3 != 0,
                 "non-live parametrized carrier");
            need(a * C0 + b * C1 + c * C2 == 0, "carrier equation");
            need(std::abs(C1) <= m1 && std::abs(C2) <= m2, "carrier roof");
            const I p0 = 3 * (b + c) - 14 * C0;
            const I p1 = 3 * (a + c) - 14 * std::abs(C1);
            const I p2 = 3 * (a + b) - 14 * std::abs(C2);
            need(p0 > 0 && p1 > 0 && p2 > 0, "positive margin");
            // Common denominator 14abc.  The cap q=3/(7c) has numerator
            // 6ab in that denominator.
            row.numerator[0] += 2 * a * std::min(p0, 6 * b);
            row.numerator[1] += 2 * b * std::min(p1, 6 * a);
            row.numerator[2] += 2 * std::min(c * p2, 6 * a * b);
            ++row.positive_carriers;
        }
    }
    return row;
}

static Row brute_row(const Triple& w) {
    const I a = w[0], b = w[1], c = w[2];
    Row row;
    row.denominator = 14 * a * b * c;
    const I m0 = (3 * (b + c) - 1) / 14;
    const I m1 = (3 * (a + c) - 1) / 14;
    const I m2 = (3 * (a + b) - 1) / 14;
    for (I C0 = 1; C0 <= m0; ++C0) {
        if (C0 % 3 == 0) continue;
        for (I C1 = -m1; C1 <= m1; ++C1) {
            if (C1 % 3 == 0) continue;
            const I residual = -a * C0 - b * C1;
            if (residual % c) continue;
            const I C2 = residual / c;
            if (std::abs(C2) > m2 || C2 % 3 == 0) continue;
            const I p0 = 3 * (b + c) - 14 * C0;
            const I p1 = 3 * (a + c) - 14 * std::abs(C1);
            const I p2 = 3 * (a + b) - 14 * std::abs(C2);
            row.numerator[0] += 2 * a * std::min(p0, 6 * b);
            row.numerator[1] += 2 * b * std::min(p1, 6 * a);
            row.numerator[2] += 2 * std::min(c * p2, 6 * a * b);
            ++row.positive_carriers;
        }
    }
    return row;
}

// For sorted positive distinct speeds these are exactly the feasible signed
// relations with coefficient magnitudes (1,1,1), (1,1,2), or (1,2,2).
static unsigned low_circuit_mask(const Triple& w) {
    const I a = w[0], b = w[1], c = w[2];
    unsigned mask = 0;
    if (c == a + b) mask |= 1U;
    if (c == 2 * a + b || c == a + 2 * b || 2 * b == a + c) mask |= 2U;
    if (a + 2 * b == 2 * c || 2 * a + b == 2 * c ||
        c == 2 * a + 2 * b || 2 * b == 2 * a + c) mask |= 4U;
    return mask;
}

static bool brute_low_circuit(const Triple& w) {
    const Triple patterns[] = {{1, 1, 1}, {1, 1, 2}, {1, 2, 2}};
    for (Triple p : patterns) {
        std::sort(p.begin(), p.end());
        do {
            for (unsigned signs = 0; signs < 8; ++signs) {
                I dot = 0;
                for (int i = 0; i < 3; ++i)
                    dot += ((signs >> i) & 1U ? -p[i] : p[i]) * w[i];
                if (dot == 0) return true;
            }
        } while (std::next_permutation(p.begin(), p.end()));
    }
    return false;
}

struct Leader {
    I numerator = -1, denominator = 1;
    Triple w{0, 0, 0}, projections{0, 0, 0};
    void update(I p, const Row& row, const Triple& speeds) {
        if (numerator < 0 || W(p) * denominator > W(numerator) * row.denominator) {
            numerator = p; denominator = row.denominator;
            w = speeds; projections = row.numerator;
        }
    }
    void print(const std::string& label) const {
        std::cout << label << " " << fraction(numerator, denominator)
                  << " w=" << triple_string(w) << " E="
                  << fraction(projections[0], denominator) << ","
                  << fraction(projections[1], denominator) << ","
                  << fraction(projections[2], denominator) << "\n";
    }
};

int main(int argc, char** argv) {
    try {
        const I H = argc > 1 ? std::stoll(argv[1]) : 199;
        need(H >= 7 && H <= 611, "height must satisfy 7<=H<=611");

        I control_rows = 0;
        for (I c = 3; c <= std::min<I>(H, 31); ++c) {
            if (c % 3 == 0) continue;
            for (I b = 2; b < c; ++b) {
                if (b % 3 == 0) continue;
                const I d = std::gcd(b, c);
                for (I a = 1; a < b; ++a) {
                    if (a % 3 == 0 || std::gcd(a, d) != 1) continue;
                    const Triple w{a, b, c};
                    const Row fast = raw_row(w), slow = brute_row(w);
                    need(fast.denominator == slow.denominator &&
                         fast.numerator == slow.numerator &&
                         fast.positive_carriers == slow.positive_carriers,
                         "brute carrier mismatch " + triple_string(w));
                    need(bool(low_circuit_mask(w)) == brute_low_circuit(w),
                         "low-circuit classification mismatch " + triple_string(w));
                    ++control_rows;
                }
            }
        }
        std::cout << "CLEANROOM ALL-PARITY RAW HEAD H=" << H
                  << " BRUTE_CONTROL_ROWS=" << control_rows << "\n";

        I rows = 0, low_rows = 0, low_overlap_rows = 0, total_positive_carriers = 0;
        I all_min_failures = 0, nonlow_any_gt = 0, nonlow_any_eq = 0, nonlow_min_ge = 0;
        std::vector<std::pair<Triple, int>> nonlow_projection_equalities;
        Triple first_failure{0, 0, 0};
        std::array<I, 8> parity_rows{}, parity_failures{};
        std::array<I, 8> circuit_profile{}, circuit_min_gt{}, circuit_min_eq{};
        Leader all_min, nonlow_min, nonlow_individual;
        std::uint64_t hash_sum = 0, hash_xor = 0;
        auto mix = [](std::uint64_t x) {
            x += 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        };

        for (I c = 3; c <= H; ++c) {
            if (c % 3 == 0) continue;
            for (I b = 2; b < c; ++b) {
                if (b % 3 == 0) continue;
                const I d = std::gcd(b, c);
                for (I a = 1; a < b; ++a) {
                    if (a % 3 == 0 || std::gcd(a, d) != 1) continue;
                    const Triple w{a, b, c};
                    const Row row = raw_row(w);
                    const unsigned low = low_circuit_mask(w);
                    const unsigned parity = unsigned((a & 1) | ((b & 1) << 1) | ((c & 1) << 2));
                    const I minimum = *std::min_element(row.numerator.begin(), row.numerator.end());
                    ++rows; ++parity_rows[parity]; ++circuit_profile[low];
                    total_positive_carriers += row.positive_carriers;
                    all_min.update(minimum, row, w);
                    const W min_cmp = W(77) * minimum - W(6) * row.denominator;
                    if (min_cmp > 0) {
                        ++all_min_failures; ++parity_failures[parity];
                        ++circuit_min_gt[low];
                        if (first_failure[2] == 0) first_failure = w;
                    }
                    if (min_cmp == 0) ++circuit_min_eq[low];
                    if (low) {
                        ++low_rows;
                        if ((low & (low - 1)) != 0) ++low_overlap_rows;
                    } else {
                        nonlow_min.update(minimum, row, w);
                        if (min_cmp >= 0) ++nonlow_min_ge;
                        for (int coordinate = 0; coordinate < 3; ++coordinate) {
                            const I p = row.numerator[coordinate];
                            nonlow_individual.update(p, row, w);
                            const W comparison = W(77) * p - W(6) * row.denominator;
                            if (comparison > 0) ++nonlow_any_gt;
                            if (comparison == 0) {
                                ++nonlow_any_eq;
                                nonlow_projection_equalities.push_back({w, coordinate});
                            }
                        }
                    }
                    std::uint64_t h = mix(std::uint64_t(a));
                    h = mix(h ^ std::uint64_t(b));
                    h = mix(h ^ std::uint64_t(c));
                    h = mix(h ^ std::uint64_t(row.denominator));
                    for (I p : row.numerator) h = mix(h ^ std::uint64_t(p));
                    h = mix(h ^ std::uint64_t(low));
                    hash_sum += h; hash_xor ^= h;
                }
            }
            if (c % 100 == 1)
                std::cerr << "progress c=" << c << " rows=" << rows << "\n";
        }

        const Triple hostile{2, 5, 7};
        const Row hostile_row = raw_row(hostile);
        std::cout << "UNIVERSE rows=" << rows << " low_rows=" << low_rows
                  << " low_overlap_rows=" << low_overlap_rows
                  << " positive_carrier_pairs=" << total_positive_carriers << "\n";
        std::cout << "TARGET all_min_gt=" << all_min_failures
                  << " nonlow_min_ge=" << nonlow_min_ge
                  << " nonlow_any_projection_gt=" << nonlow_any_gt
                  << " nonlow_any_projection_eq=" << nonlow_any_eq
                  << " first_failure=" << triple_string(first_failure) << "\n";
        for (const auto& equality : nonlow_projection_equalities)
            std::cout << "NONLOW_PROJECTION_EQUALITY w=" << triple_string(equality.first)
                      << " coordinate=" << equality.second << "\n";
        std::cout << "HOSTILE_2_5_7 mask=" << low_circuit_mask(hostile) << " E="
                  << fraction(hostile_row.numerator[0], hostile_row.denominator) << ","
                  << fraction(hostile_row.numerator[1], hostile_row.denominator) << ","
                  << fraction(hostile_row.numerator[2], hostile_row.denominator)
                  << " positive_carrier_pairs=" << hostile_row.positive_carriers << "\n";
        all_min.print("MAX_ALL_MIN");
        nonlow_min.print("MAX_NONLOW_MIN");
        nonlow_individual.print("MAX_NONLOW_INDIVIDUAL");
        for (unsigned mask = 0; mask < 8; ++mask)
            if (parity_rows[mask])
                std::cout << "PARITY mask=" << mask << " rows=" << parity_rows[mask]
                          << " min_failures=" << parity_failures[mask] << "\n";
        for (unsigned mask = 0; mask < 8; ++mask)
            if (circuit_profile[mask])
                std::cout << "CIRCUITS mask=" << mask << " rows=" << circuit_profile[mask]
                          << " min_gt=" << circuit_min_gt[mask]
                          << " min_eq=" << circuit_min_eq[mask] << "\n";
        std::cout << "COMMUTATIVE_HASH sum=" << hash_sum << " xor=" << hash_xor << "\n";

        need(nonlow_min_ge == 0 && nonlow_any_gt == 0,
             "all-parity nonlow weak target failure");
        need(first_failure == hostile, "wrong first failure");
        if (H == 199) need(rows == 333960, "H199 row count");
        if (H == 611) {
            need(rows == 9720930, "H611 row count");
            need(all_min_failures == 14220, "H611 failure count");
        }
        std::cout << "PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << "\n";
        return 1;
    }
    return 0;
}
