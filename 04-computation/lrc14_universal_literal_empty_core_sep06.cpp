// FINITE-EXACT independent native interval verification of the network head.
// Universe: primitive 1<=a<b<c<=H, odd and nonzero modulo three, H<=601.
// Primary engine: all six literal sheet assignments, exact endpoints on
// denominator 42abc, and a simultaneous three-pointer contact scan.
// No carrier congruence, lattice-direction classifier, or raw roof formula
// is used to compute the projections.
//
// Definitions independently transcribed from the native sheet construction in
// 04-computation/lrc14_one_ray_overnight_hexagon_sep05.py.
// Reproduce: c++ -std=c++17 -O3 -DNDEBUG <this.cpp> -o /tmp/lrc14_literal_sep06
//           /tmp/lrc14_literal_sep06 601
// Optional second argument is an exact per-row TSV output path.

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using I = std::int64_t;
using V = std::array<I, 3>;
using Wide = __int128_t;

static void need(bool value, const std::string& label) {
    if (!value) throw std::runtime_error(label);
}

static std::string fraction(I p, I q) {
    const I g = std::gcd(p, q);
    return std::to_string(p / g) + "/" + std::to_string(q / g);
}

static std::string triple(const V& w) {
    return std::to_string(w[0]) + "," + std::to_string(w[1]) + "," + std::to_string(w[2]);
}

struct Cursor {
    I left, right, width;
    I step, radius, denominator;
    int index, speed, color;
    bool valid;

    Cursor(int speed_, int color_, I multiplier, I denominator_)
        : step(42 * multiplier), radius(3 * multiplier), denominator(denominator_),
          index(0), speed(speed_), color(color_), valid(true) {
        // Sorted residues (3*owner-speed*color) mod(3*speed) are precisely
        // r0, r0+3, ..., r0+3*(speed-1), r0=(-speed*color) mod3.
        if (color == 0) {
            left = 0;
            right = radius;
        } else {
            const I residue = (3 - (speed * color) % 3) % 3;
            const I center = 14 * residue * multiplier;
            left = center - radius;
            right = center + radius;
        }
        width = right - left;
    }

    void advance() {
        ++index;
        if (color == 0) {
            if (index > speed) {
                valid = false;
                return;
            }
            if (index == speed) {
                left = denominator - radius;
                right = denominator;
            } else {
                left = index * step - radius;
                right = index * step + radius;
            }
        } else {
            if (index == speed) {
                valid = false;
                return;
            }
            left += step;
            right += step;
        }
        width = right - left;
    }
};

struct Native {
    I denominator;
    V projections{};
    I mass = 0;
    I contacts = 0;
};

static Native literal(const V& w) {
    Native result;
    result.denominator = 42 * w[0] * w[1] * w[2];
    constexpr int permutations[6][3] = {
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2},
        {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
    };
    for (const auto& assignment : permutations) {
        Cursor x(static_cast<int>(w[0]), assignment[0], result.denominator / (42 * w[0]), result.denominator);
        Cursor y(static_cast<int>(w[1]), assignment[1], result.denominator / (42 * w[1]), result.denominator);
        Cursor z(static_cast<int>(w[2]), assignment[2], result.denominator / (42 * w[2]), result.denominator);
        while (x.valid && y.valid && z.valid) {
            const I left = std::max(x.left, std::max(y.left, z.left));
            const I right = std::min(x.right, std::min(y.right, z.right));
            if (left < right) {
                // All three intervals meet. For each omitted sheet, use the
                // complete intersection length of the other two intervals,
                // capped by the complete omitted interval length.
                const I yz = std::min(y.right, z.right) - std::max(y.left, z.left);
                const I xz = std::min(x.right, z.right) - std::max(x.left, z.left);
                const I xy = std::min(x.right, y.right) - std::max(x.left, y.left);
                result.projections[0] += std::min(yz, x.width);
                result.projections[1] += std::min(xz, y.width);
                result.projections[2] += std::min(xy, z.width);
                result.mass += right - left;
                ++result.contacts;
            }
            // At least one interval ends. Ties advance simultaneously and
            // zero-length contacts never contribute.
            if (x.right == right) x.advance();
            if (y.right == right) y.advance();
            if (z.right == right) z.advance();
        }
    }
    return result;
}

struct Leader {
    I numerator = 0, denominator = 1;
    V speeds{};
    V projections{};
    I contacts = 0;

    void update(I numerator_, const Native& row, const V& w) {
        if (Wide(numerator_) * denominator > Wide(numerator) * row.denominator) {
            numerator = numerator_;
            denominator = row.denominator;
            speeds = w;
            projections = row.projections;
            contacts = row.contacts;
        }
    }

    void print(const std::string& label) const {
        std::cout << label << " " << fraction(numerator, denominator)
                  << " w=" << triple(speeds) << " E="
                  << fraction(projections[0], denominator) << ","
                  << fraction(projections[1], denominator) << ","
                  << fraction(projections[2], denominator)
                  << " native_contacts=" << contacts << "\n";
    }
};

static void controls() {
    struct Control { V w; V p; V q; I mp, mq, contacts; };
    const Control cases[] = {
        {{1, 5, 7}, {8, 6, 4}, {245, 49, 35}, 8, 245, 2},
        {{1, 5, 11}, {6, 6, 6}, {77, 77, 77}, 6, 77, 2},
        {{1, 19, 23}, {12, 12, 12}, {437, 161, 161}, 12, 437, 4},
        {{17, 23, 25}, {106, 12, 2546}, {4025, 425, 68425}, 744, 68425, 4},
        {{19, 23, 29}, {156, 192, 3840}, {4669, 3857, 88711}, 154, 12673, 6},
        {{1, 1201, 1205}, {310116, 516, 516}, {10130435, 8435, 8435}, 310116, 10130435, 172},
        {{1, 599, 607}, {38940, 132, 132}, {2545151, 4249, 4249}, 38940, 2545151, 44},
        {{5, 1001, 1003}, {156680, 122, 848864}, {7028021, 5015, 35140105}, 767384, 35140105, 58}
    };
    for (const auto& c : cases) {
        const Native row = literal(c.w);
        for (int i = 0; i < 3; ++i)
            need(Wide(row.projections[i]) * c.q[i] == Wide(c.p[i]) * row.denominator,
                 "known literal/raw control mismatch " + triple(c.w));
        need(Wide(row.mass) * c.mq == Wide(c.mp) * row.denominator,
             "known physical mass mismatch " + triple(c.w));
        need(row.contacts == 3 * c.contacts, "native/raw contact multiplicity mismatch " + triple(c.w) + " native=" + std::to_string(row.contacts));
        std::cout << "CONTROL w=" << triple(c.w) << " E="
                  << fraction(row.projections[0], row.denominator) << ","
                  << fraction(row.projections[1], row.denominator) << ","
                  << fraction(row.projections[2], row.denominator)
                  << " mass=" << fraction(row.mass, row.denominator)
                  << " native_contacts=" << row.contacts << "\n";
    }
}

int main(int argc, char** argv) {
    try {
        const int H = argc > 1 ? std::stoi(argv[1]) : 601;
        need(H >= 5 && H <= 601, "height must satisfy 5<=H<=601");
        std::ofstream dump;
        if (argc > 2) {
            dump.open(argv[2]);
            need(bool(dump), "cannot open optional exact-row dump");
            dump << "a\tb\tc\tdenominator\tE0_numerator\tE1_numerator\tE2_numerator\tmass_numerator\tcontacts\n";
        }
        std::cout << "FINITE-EXACT independent native six-sheet interval engine H=" << H << "\n";
        controls();
        std::cout.flush();

        I rows = 0, norm4 = 0, empty = 0, strict_rows = 0;
        I all_min_failures = 0, non_norm4_every_failures = 0;
        I contacts_total = 0;
        std::vector<V> equalities;
        Leader all_min, non_norm4_max, non_norm4_min;
        std::uint64_t checksum = 14695981039346656037ULL;
        auto hash_integer = [&](I x) {
            auto value = static_cast<std::uint64_t>(x);
            for (int i = 0; i < 8; ++i) {
                checksum ^= value & 255U;
                checksum *= 1099511628211ULL;
                value >>= 8;
            }
        };
        for (I c = 5; c <= H; c += 2) {
            if (c % 3 == 0) continue;
            for (I b = 3; b < c; b += 2) {
                if (b % 3 == 0) continue;
                const I bc_gcd = std::gcd(b, c);
                for (I a = 1; a < b; a += 2) {
                    if (a % 3 == 0 || std::gcd(a, bc_gcd) != 1) continue;
                    const V w{a, b, c};
                    const Native row = literal(w);
                    ++rows;
                    contacts_total += row.contacts;
                    empty += row.contacts == 0;
                    const I minimum = *std::min_element(row.projections.begin(), row.projections.end());
                    const I maximum = *std::max_element(row.projections.begin(), row.projections.end());
                    all_min.update(minimum, row, w);
                    const Wide min_compare = Wide(77) * minimum - Wide(6) * row.denominator;
                    if (min_compare > 0) ++all_min_failures;
                    if (min_compare == 0) equalities.push_back(w);
                    // This is only a speed-identity classification, not an
                    // inferred carrier support or congruence filter. Every row
                    // above was evaluated by the full native interval engine.
                    const bool has_norm4 = c == 2 * a + b || c == a + 2 * b || 2 * b == a + c;
                    if (has_norm4) {
                        ++norm4;
                    } else {
                        ++strict_rows;
                        non_norm4_max.update(maximum, row, w);
                        non_norm4_min.update(minimum, row, w);
                        if (Wide(77) * maximum >= Wide(6) * row.denominator)
                            ++non_norm4_every_failures;
                    }
                    for (I value : w) hash_integer(value);
                    hash_integer(row.denominator);
                    for (I value : row.projections) hash_integer(value);
                    hash_integer(row.mass);
                    hash_integer(row.contacts);
                    if (dump.is_open())
                        dump << a << '\t' << b << '\t' << c << '\t' << row.denominator
                             << '\t' << row.projections[0] << '\t' << row.projections[1]
                             << '\t' << row.projections[2] << '\t' << row.mass << '\t' << row.contacts << '\n';
                }
            }
            if (c % 100 == 1) {
                std::cerr << "native interval progress c=" << c << " rows=" << rows << "\n";
            }
        }
        std::cout << "UNIVERSE rows=" << rows << " norm4_speed_identity=" << norm4
                  << " without_norm4_identity=" << strict_rows << " empty=" << empty
                  << " native_contacts_total=" << contacts_total << "\n";
        std::cout << "TARGET all_min_gt_6/77=" << all_min_failures
                  << " non_norm4_any_E_ge_6/77=" << non_norm4_every_failures << "\n";
        all_min.print("MAX_ALL_MIN_E");
        non_norm4_max.print("MAX_NON_NORM4_E");
        non_norm4_min.print("MAX_NON_NORM4_MIN_E");
        std::cout << "MIN_EQUALITIES count=" << equalities.size();
        for (const auto& w : equalities) std::cout << " w=" << triple(w);
        std::cout << "\nSEMANTIC_FNV1A64 " << checksum << " (noncryptographic exact-row checksum)\n";
        need(all_min_failures == 0 && non_norm4_every_failures == 0, "finite-head target failed");
        const std::vector<V> expected = H >= 11 ? std::vector<V>{{1, 5, 11}} : std::vector<V>{};
        need(equalities == expected, "unexpected equality locus");
        if (H == 79) need(rows == 2910, "H79 universe count");
        if (H == 601) need(rows == 1317935, "H601 universe count");
        std::cout << "PASS: every eligible row checked from literal intervals; strict every-projection bound"
                     " without norm4 identity; universal minimum bound and stated equality locus.\n";
        std::cout << "SCOPE: finite network head only; no all-height inference or LRC14 conclusion.\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << "\n";
        return 1;
    }
    return 0;
}
