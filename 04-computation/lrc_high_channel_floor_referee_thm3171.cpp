#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

using i128 = __int128_t;
using i64 = long long;

struct Moments { i128 s0, s1, s2; };

static Moments floor_moments(i64 n, i64 m, i64 a, i64 b) {
    if (n == 0) return {0, 0, 0};
    const i128 sn = (i128)n * (n - 1) / 2;
    const i128 sn2 = (i128)n * (n - 1) * (2 * (i128)n - 1) / 6;
    const i64 qa = a / m, ar = a % m;
    const i64 qb = b / m, br = b % m;
    i128 b0 = (i128)qa * sn + (i128)qb * n;
    i128 b1 = (i128)qa * sn2 + (i128)qb * sn;
    i128 b2 = (i128)qa * qa * sn2 + 2 * (i128)qa * qb * sn
              + (i128)qb * qb * n;
    if (ar == 0) return {b0, b1, b2};
    const i64 y = (ar * (n - 1) + br) / m;
    if (y == 0) return {b0, b1, b2};
    const Moments u = floor_moments(y, ar, m, m - br + ar - 1);
    const i128 r0 = (i128)n * y - u.s0;
    const i128 r1 = (i128)y * sn - (u.s2 - u.s0) / 2;
    const i128 r2 = (i128)n * y * y - 2 * u.s1 - u.s0;
    return {b0 + r0, b1 + r1,
            b2 + 2 * (i128)qa * r1 + 2 * (i128)qb * r0 + r2};
}

struct Prefix { i128 count, sum; };

static Prefix residue_prefix(i64 n, i64 m, i64 a, i64 b, i64 t,
                             const Moments &base) {
    // Count and sum of residues (a*k+b) mod m that are < t.
    const Moments shifted = floor_moments(n, m, a, b + m - t);
    const i128 d0 = shifted.s0 - base.s0;
    const i128 d1 = shifted.s1 - base.s1;
    const i128 y0d = (shifted.s2 - base.s2 - d0) / 2;
    const i128 high_sum = (i128)a * d1 + (i128)b * d0 - (i128)m * y0d;
    const i128 total = (i128)a * n * (n - 1) / 2 + (i128)b * n
                       - (i128)m * base.s0;
    return {(i128)n - d0, total - high_sum};
}

static i128 triangle_sum(i64 n, i64 m, i64 a, i64 b, i64 C,
                         const Moments &base, i128 total) {
    if (C <= 0) return 0;
    const i64 T = (C - 1) / 168;
    if (2 * T >= m) throw std::runtime_error("overlapping triangle tails");
    const Prefix low = residue_prefix(n, m, a, b, T + 1, base);
    const Prefix before_high = residue_prefix(n, m, a, b, m - T, base);
    const i128 high_count = (i128)n - before_high.count;
    const i128 high_sum = total - before_high.sum;
    return (i128)C * low.count - 168 * low.sum
           + ((i128)C - 168 * (i128)m) * high_count + 168 * high_sum;
}

struct Fraction { i128 num, den; };

static Fraction pair_mass(int e, i64 p, int f, i64 q) {
    if (p > q) return pair_mass(f, q, e, p);
    const i64 z = 168 * p - e, w = 168 * q - f;
    const i64 r = (90 * e) % 168, s = (90 * f) % 168;
    const i64 determinant = r * w - s * z;
    if (determinant % 168 != 0) throw std::runtime_error("nonintegral phase");
    i64 b = (determinant / 168) % z;
    if (b < 0) b += z;
    const i64 a = w % z;
    const Moments base = floor_moments(p, z, a, b);
    const i128 total = (i128)a * p * (p - 1) / 2 + (i128)b * p
                       - (i128)z * base.s0;
    const i128 outer = triangle_sum(p, z, a, b, 12 * (z + w), base, total);
    const i128 inner = triangle_sum(p, z, a, b, 12 * (w - z), base, total);
    return {outer - inner, (i128)z * w};
}

static std::string text128(i128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string s;
    while (x) { s.push_back(char('0' + x % 10)); x /= 10; }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

struct Witness {
    bool set = false;
    Fraction value{0, 1};
    i64 p = 0, q = 0;
    int edge = 0, orientation = 0;
};

static bool less_fraction(const Fraction &x, const Fraction &y) {
    return x.num * y.den < y.num * x.den;
}

int main(int argc, char **argv) {
    const i64 maximum = argc > 1 ? std::stoll(argv[1]) : 3000;
    const int labels[6] = {1, 2, 3, 4, 6, 12};
    const int edge_i[9] = {0,0,0,0,0,1,1,1,1};
    const int edge_j[9] = {1,2,3,4,5,2,3,4,5};
    unsigned threads_n = std::max(1u, std::thread::hardware_concurrency());
    if (argc > 2) threads_n = std::stoul(argv[2]);
    std::atomic<i64> next_p{6};
    std::atomic<bool> failed{false};
    std::mutex lock;
    Witness failure;
    std::vector<Witness> minima(18);
    std::atomic<unsigned long long> pairs{0}, checks{0}, exceptional{0};

    auto worker = [&]() {
        std::vector<Witness> local(18);
        unsigned long long local_pairs = 0, local_checks = 0, local_exceptional = 0;
        while (!failed.load(std::memory_order_relaxed)) {
            const i64 p = next_p.fetch_add(1);
            if (p > maximum) break;
            for (i64 q = p + 1; q <= maximum && q <= 2 * p; ++q) {
                const i64 g = std::gcd(p, q), P = p / g, Q = q / g;
                if (P + Q < 8) continue;
                ++local_pairs;
                for (int edge = 0; edge < 9; ++edge) {
                    const int u = edge_i[edge], v = edge_j[edge];
                    for (int orientation = 0; orientation < 2; ++orientation) {
                        Fraction value = orientation == 0
                            ? pair_mass(labels[u], p, labels[v], q)
                            : pair_mass(labels[u], q, labels[v], p);
                        const int index = 2 * edge + orientation;
                        if (!local[index].set || less_fraction(value, local[index].value)) {
                            local[index] = {true, value, p, q, edge, orientation};
                        }
                        ++local_checks;
                        const bool exception = edge == 3 && P == 3 && Q == 5;
                        if (exception) {
                            ++local_exceptional;
                            if (value.num * 280393 < (i128)2030 * value.den) {
                                std::lock_guard<std::mutex> guard(lock);
                                failed = true;
                                failure = {true, value, p, q, edge, orientation};
                                break;
                            }
                        } else if (value.num * 105 < value.den) {
                            std::lock_guard<std::mutex> guard(lock);
                            failed = true;
                            failure = {true, value, p, q, edge, orientation};
                            break;
                        }
                    }
                    if (failed) break;
                }
                if (failed) break;
            }
        }
        std::lock_guard<std::mutex> guard(lock);
        pairs += local_pairs; checks += local_checks; exceptional += local_exceptional;
        for (int i = 0; i < 18; ++i) {
            if (local[i].set && (!minima[i].set || less_fraction(local[i].value, minima[i].value)))
                minima[i] = local[i];
        }
    };

    std::vector<std::thread> pool;
    for (unsigned i = 0; i < threads_n; ++i) pool.emplace_back(worker);
    for (auto &thread : pool) thread.join();
    std::cout << "max_level=" << maximum << ";threads=" << threads_n
              << ";ordered_level_pairs=" << pairs << ";edge_orientation_checks=" << checks
              << ";exception_checks=" << exceptional << "\n";
    if (failed) {
        std::cout << "COUNTEREXAMPLE edge=" << failure.edge
                  << ";orientation=" << failure.orientation << ";p=" << failure.p
                  << ";q=" << failure.q << ";num=" << text128(failure.value.num)
                  << ";den=" << text128(failure.value.den) << "\n";
        return 1;
    }
    std::cout << "PASS\n";
    for (int i = 0; i < 18; ++i) {
        const Witness &w = minima[i];
        std::cout << "edge_orientation=" << i << ";p=" << w.p << ";q=" << w.q
                  << ";num=" << text128(w.value.num) << ";den=" << text128(w.value.den)
                  << "\n";
    }
    return 0;
}
