// Exact finite bank for THM-1231.
//
// This program contains no floating-point arithmetic.  It proves the finite
// implication used by THM-1231: every projective five-tuple carried by a
// spanning tree of heavy ratios has pair mass at least 13/81.  A pair is
// heavy when rho < 1/49 - 71/63504.  The analytic part of THM-1231 proves
// that every putative counterexample has a connected heavy graph and that
// every reduced heavy ratio a:b has b <= 127.
//
// Compile and run, for example, with
//   clang++ -O3 -std=c++17 this_file.cpp -o /tmp/thm1202_bank
//   /tmp/thm1202_bank

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

using std::array;
using std::cerr;
using std::cout;
using std::gcd;
using std::hash;
using std::lock_guard;
using std::max;
using std::mutex;
using std::pair;
using std::reverse;
using std::sort;
using std::string;
using std::thread;
using std::uint64_t;
using std::uint8_t;
using std::unordered_set;
using std::vector;

using u64 = uint64_t;
using u128 = unsigned __int128;

static constexpr u64 Q = 1000000000ULL;

struct Key {
    array<u64, 5> a{};
    uint8_t n = 0;

    bool operator==(const Key &other) const {
        return n == other.n && a == other.a;
    }
};

struct KeyHash {
    std::size_t operator()(const Key &key) const {
        std::size_t h = key.n;
        for (int i = 0; i < key.n; ++i) {
            std::size_t z = hash<u64>{}(key.a[i]);
            h ^= z + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

bool key_less(const Key &left, const Key &right) {
    return std::lexicographical_compare(
        left.a.begin(), left.a.begin() + left.n,
        right.a.begin(), right.a.begin() + right.n);
}

struct PairData {
    u64 numerator;
    u64 denominator;
};

PairData pair_data(u64 a, u64 b) {
    u64 g = gcd(a, b);
    a /= g;
    b /= g;
    if (a > b) {
        std::swap(a, b);
    }
    auto folded = [](u64 z) {
        u64 r = z % 14;
        return r * (14 - r);
    };
    int correction = static_cast<int>(folded(a + b))
        - static_cast<int>(folded(b - a));
    u64 product = a * b;
    long long signed_numerator = static_cast<long long>(4 * product)
        + correction;
    return {
        static_cast<u64>(signed_numerator),
        196 * product,
    };
}

u64 floor_q_rho(u64 a, u64 b) {
    PairData data = pair_data(a, b);
    return static_cast<u64>(
        static_cast<u128>(Q) * data.numerator / data.denominator);
}

bool defect_gt_theta(u64 a, u64 b) {
    // (1/49-rho) > 71/63504, compared by exact cross multiplication.
    PairData data = pair_data(a, b);
    if (data.denominator <= 49 * data.numerator) {
        return false;
    }
    return static_cast<u128>(data.denominator - 49 * data.numerator)
            * 63504
        > static_cast<u128>(71) * 49 * data.denominator;
}

u64 floor_sum(const Key &key) {
    u64 result = 0;
    for (int i = 0; i < key.n; ++i) {
        for (int j = i + 1; j < key.n; ++j) {
            result += floor_q_rho(key.a[i], key.a[j]);
        }
    }
    return result;
}

bool certified_at_least_target(u64 lower_numerator, int remaining_pairs) {
    // lower_numerator/Q + remaining_pairs/91 >= 13/81.
    return (static_cast<u128>(lower_numerator) * 91
                + static_cast<u128>(remaining_pairs) * Q)
            * 81
        >= static_cast<u128>(13) * Q * 91;
}

Key add_ratio(const Key &key, int parent, u64 p, u64 q) {
    // Add y = key[parent] * p/q, then take the primitive sorted integer row.
    u128 raw_numerator = static_cast<u128>(key.a[parent]) * p;
    if (raw_numerator > std::numeric_limits<u64>::max()) {
        cerr << "numerator overflow\n";
        std::exit(2);
    }
    u64 numerator = static_cast<u64>(raw_numerator);
    u64 denominator = q;
    u64 common = gcd(numerator, denominator);
    numerator /= common;
    denominator /= common;

    Key result;
    result.n = key.n + 1;
    for (int i = 0; i < key.n; ++i) {
        u128 value = static_cast<u128>(key.a[i]) * denominator;
        if (value > std::numeric_limits<u64>::max()) {
            cerr << "coordinate overflow\n";
            std::exit(2);
        }
        result.a[i] = static_cast<u64>(value);
    }
    result.a[key.n] = numerator;

    u64 row_gcd = 0;
    for (int i = 0; i < result.n; ++i) {
        row_gcd = gcd(row_gcd, result.a[i]);
    }
    for (int i = 0; i < result.n; ++i) {
        result.a[i] /= row_gcd;
    }
    sort(result.a.begin(), result.a.begin() + result.n);
    return result;
}

bool distinct(const Key &key) {
    for (int i = 1; i < key.n; ++i) {
        if (key.a[i] == key.a[i - 1]) {
            return false;
        }
    }
    return true;
}

string decimal(u128 value) {
    if (value == 0) {
        return "0";
    }
    string result;
    while (value != 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    reverse(result.begin(), result.end());
    return result;
}

bool pair_is(u64 a, u64 b, u64 numerator, u64 denominator) {
    PairData data = pair_data(a, b);
    return static_cast<u128>(data.numerator) * denominator
        == static_cast<u128>(numerator) * data.denominator;
}

int main() {
    vector<pair<u64, u64>> heavy_types;
    for (u64 b = 2; b <= 127; ++b) {
        for (u64 a = 1; a < b; ++a) {
            if (gcd(a, b) == 1 && defect_gt_theta(a, b)) {
                heavy_types.push_back({b, a});
            }
        }
    }
    cout << "heavy_threshold=71/63504 denominator_cutoff=127 types="
         << heavy_types.size() << "\n";
    if (heavy_types.size() != 61
        || heavy_types.back().first != 97) {
        return 3;
    }

    Key singleton;
    singleton.n = 1;
    singleton.a[0] = 1;
    unordered_set<Key, KeyHash> states;
    states.insert(singleton);

    for (int n = 2; n <= 4; ++n) {
        unordered_set<Key, KeyHash> next;
        next.reserve(states.size() * 10);
        u64 attempts = 0;
        u64 pruned = 0;
        for (const Key &state : states) {
            for (int parent = 0; parent < state.n; ++parent) {
                for (auto [p, q] : heavy_types) {
                    for (int reverse_ratio = 0; reverse_ratio < 2;
                         ++reverse_ratio) {
                        Key row = add_ratio(
                            state,
                            parent,
                            reverse_ratio ? q : p,
                            reverse_ratio ? p : q);
                        if (!distinct(row)) {
                            continue;
                        }
                        ++attempts;
                        int remaining_pairs = 10 - n * (n - 1) / 2;
                        if (certified_at_least_target(
                                floor_sum(row), remaining_pairs)) {
                            ++pruned;
                            continue;
                        }
                        next.insert(row);
                    }
                }
            }
        }
        states.swap(next);
        cout << "level=" << n
             << " states=" << states.size()
             << " attempts=" << attempts
             << " pruned=" << pruned << "\n";
    }

    vector<Key> four_rows(states.begin(), states.end());
    unsigned thread_count = max(1u, thread::hardware_concurrency());
    vector<thread> workers;
    mutex result_mutex;
    u64 checked = 0;
    u64 passed = 0;
    u64 global_lower = 10 * Q;
    Key global_key;
    unordered_set<Key, KeyHash> fallback_rows;

    for (unsigned worker = 0; worker < thread_count; ++worker) {
        workers.emplace_back([&, worker]() {
            u64 local_checked = 0;
            u64 local_passed = 0;
            u64 local_lower = 10 * Q;
            Key local_key;
            vector<Key> local_fallback;
            for (std::size_t row_index = worker;
                 row_index < four_rows.size();
                 row_index += thread_count) {
                const Key &state = four_rows[row_index];
                for (int parent = 0; parent < 4; ++parent) {
                    for (auto [p, q] : heavy_types) {
                        for (int reverse_ratio = 0; reverse_ratio < 2;
                             ++reverse_ratio) {
                            Key row = add_ratio(
                                state,
                                parent,
                                reverse_ratio ? q : p,
                                reverse_ratio ? p : q);
                            if (!distinct(row)) {
                                continue;
                            }
                            ++local_checked;
                            u64 lower = floor_sum(row);
                            if (lower < local_lower
                                || (lower == local_lower
                                    && key_less(row, local_key))) {
                                local_lower = lower;
                                local_key = row;
                            }
                            if (certified_at_least_target(lower, 0)) {
                                ++local_passed;
                            } else {
                                local_fallback.push_back(row);
                            }
                        }
                    }
                }
            }
            lock_guard<mutex> guard(result_mutex);
            checked += local_checked;
            passed += local_passed;
            for (const Key &row : local_fallback) {
                fallback_rows.insert(row);
            }
            if (local_lower < global_lower
                || (local_lower == global_lower
                    && key_less(local_key, global_key))) {
                global_lower = local_lower;
                global_key = local_key;
            }
        });
    }
    for (thread &worker : workers) {
        worker.join();
    }

    cout << "five_extensions=" << checked
         << " certified=" << passed
         << " failures=" << checked - passed << "\n";
    cout << "minimum_fixed_lower_numerator=" << global_lower
         << " / " << Q << " witness=";
    for (int i = 0; i < 5; ++i) {
        cout << (i ? "," : "") << global_key.a[i];
    }
    cout << "\n";

    if (checked != 347319154 || passed != checked || !fallback_rows.empty()) {
        return 4;
    }

    cout << "witness_pair_word="
         << "2*(1/63+17/819+1/84+1/91+23/1092)"
         << " exact_R=44/273 target_gap=5/7371\n";
    u128 target_margin = static_cast<u128>(global_lower) * 81
        - static_cast<u128>(13) * Q;
    cout << "fixed_grid_target_margin_numerator=" << decimal(target_margin)
         << " / " << 81 * Q << "\n";
    return 0;
}
