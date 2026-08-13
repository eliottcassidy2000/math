// Exact finite-head compiler for the disconnected-low LRC14 pair floor.
//
// Scope: the 2,530 body-safe contexts on the 29 rulers L<4592 and every
// primitive non-(3,5) channel P<Q<8P, P+Q>=8, at common dilation g=1,2,3,
// with raw p=gP<264.  This file does not implement or use the affine-tail
// argument.  It uses the already-proved THM-3350 midpoint lower bound as an
// exact screen and evaluates every unresolved row with an independent
// __int128 port of the THM-3352 Euclidean physical-mass engine.

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

namespace {

constexpr int SCALE_BITS = 56;
constexpr u64 SCALE = u64{1} << SCALE_BITS;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

std::string text_i128(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    u128 magnitude = negative ? u128(-(value + 1)) + 1 : u128(value);
    std::string answer;
    while (magnitude) {
        answer.push_back(char('0' + magnitude % 10));
        magnitude /= 10;
    }
    if (negative) answer.push_back('-');
    std::reverse(answer.begin(), answer.end());
    return answer;
}

i64 checked_i64(i128 value, const char* label) {
    if (value < std::numeric_limits<i64>::min() || value > std::numeric_limits<i64>::max()) {
        fail(std::string("int64 overflow in ") + label + ": " + text_i128(value));
    }
    return i64(value);
}

i64 positive_mod(i64 value, i64 modulus) {
    i64 answer = value % modulus;
    return answer < 0 ? answer + modulus : answer;
}

struct Fraction {
    i64 numerator = 0;
    i64 denominator = 1;

    Fraction() = default;
    Fraction(i64 n, i64 d) : numerator(n), denominator(d) {
        require(d != 0, "zero fraction denominator");
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
        i64 divisor = std::gcd(numerator < 0 ? -numerator : numerator, denominator);
        numerator /= divisor;
        denominator /= divisor;
    }
};

bool fraction_less(const Fraction& left, const Fraction& right) {
    return i128(left.numerator) * right.denominator < i128(right.numerator) * left.denominator;
}

std::string fraction_text(const Fraction& value) {
    if (value.denominator == 1) return std::to_string(value.numerator);
    return std::to_string(value.numerator) + "/" + std::to_string(value.denominator);
}

struct Moments {
    i64 zero;
    i64 one;
    i64 square;
};

Moments floor_moments(i64 n, i64 m, i64 a, i64 b) {
    if (n == 0) return {0, 0, 0};
    require(n >= 0 && m > 0, "invalid floor_moments domain");
    i128 s1 = i128(n) * (n - 1) / 2;
    i128 s2 = i128(n) * (n - 1) * (2 * n - 1) / 6;
    i64 qa = a / m, a0 = a % m;
    i64 qb = b / m, b0 = b % m;
    if (a0 < 0) { --qa; a0 += m; }
    if (b0 < 0) { --qb; b0 += m; }
    i128 base0 = i128(qa) * s1 + i128(qb) * n;
    i128 base1 = i128(qa) * s2 + i128(qb) * s1;
    i128 base2 = i128(qa) * qa * s2 + i128(2) * qa * qb * s1 + i128(qb) * qb * n;
    if (a0 == 0) {
        return {checked_i64(base0, "floor base0"), checked_i64(base1, "floor base1"),
                checked_i64(base2, "floor base2")};
    }
    i64 height = checked_i64((i128(a0) * (n - 1) + b0) / m, "floor height");
    if (height == 0) {
        return {checked_i64(base0, "floor zero base0"), checked_i64(base1, "floor zero base1"),
                checked_i64(base2, "floor zero base2")};
    }
    Moments upper = floor_moments(height, a0, m, m - b0 + a0 - 1);
    i128 r0 = i128(n) * height - upper.zero;
    i128 r1 = i128(height) * s1 - (i128(upper.square) - upper.zero) / 2;
    i128 r2 = i128(n) * height * height - i128(2) * upper.one - upper.zero;
    return {
        checked_i64(base0 + r0, "floor recurse0"),
        checked_i64(base1 + r1, "floor recurse1"),
        checked_i64(base2 + i128(2) * qa * r1 + i128(2) * qb * r0 + r2, "floor recurse2")
    };
}

struct Prefix {
    i64 count;
    i64 sum;
};

Prefix residue_prefix(i64 n, i64 m, i64 a, i64 b, i64 threshold,
                      const Moments& base, i64 total) {
    if (threshold <= 0) return {0, 0};
    if (threshold >= m) return {n, total};
    Moments shifted = floor_moments(n, m, a, b + m - threshold);
    i64 d0 = shifted.zero - base.zero;
    i64 d1 = shifted.one - base.one;
    i64 y0d = (shifted.square - base.square - d0) / 2;
    i64 high_sum = checked_i64(i128(a) * d1 + i128(b) * d0 - i128(m) * y0d,
                               "residue high sum");
    return {n - d0, total - high_sum};
}

i64 triangle_sum(i64 n, i64 m, i64 a, i64 b, i64 peak, i64 L,
                 const Moments& base, i64 total) {
    if (peak <= 0) return 0;
    i64 radius = (peak - 1) / L;
    i64 turns = radius / m;
    i64 tail = radius % m;
    i128 answer = i128(n) * (i128(2) * turns * peak - i128(L) * m * turns * turns);
    Prefix low = residue_prefix(n, m, a, b, tail + 1, base, total);
    answer += i128(peak - L * turns * m) * low.count - i128(L) * low.sum;
    Prefix before = residue_prefix(n, m, a, b, m - tail, base, total);
    i64 high_count = n - before.count;
    i64 high_sum = total - before.sum;
    answer += i128(peak - L * (turns + 1) * m) * high_count + i128(L) * high_sum;
    return checked_i64(answer, "triangle sum");
}

struct MassFraction {
    i64 numerator;
    i64 denominator;
};

MassFraction physical_mass(i64 L, i64 cell, i64 e, i64 p, i64 f, i64 q) {
    require(L >= 168 && L % 14 == 0 && e >= 1 && e <= 14 && f >= 1 && f <= 14,
            "physical mass outside ruler domain");
    require(p >= 1 && q >= 1, "nonpositive level");
    i64 z = L * p - e, w = L * q - f;
    if (z > w) return physical_mass(L, cell, f, q, e, p);
    i64 r = positive_mod(e * cell, L), s = positive_mod(f * cell, L);
    i64 determinant = checked_i64(i128(r) * w - i128(s) * z, "phase determinant");
    require(determinant % L == 0, "nonintegral reflected phase");
    i64 b = positive_mod(determinant / L, z), a = w % z;
    Moments base = floor_moments(p, z, a, b);
    i64 total = checked_i64(i128(a) * p * (p - 1) / 2 + i128(b) * p - i128(z) * base.zero,
                            "residue total");
    i64 unit = L / 14;
    i64 outer = checked_i64(i128(unit) * (z + w), "outer peak");
    i64 inner = checked_i64(i128(unit) * (w - z), "inner peak");
    i64 numerator = triangle_sum(p, z, a, b, outer, L, base, total)
                  - triangle_sum(p, z, a, b, inner, L, base, total);
    i64 denominator = checked_i64(i128(z) * w, "mass denominator");
    require(numerator >= 0 && denominator > 0, "negative physical mass");
    i64 divisor = std::gcd(numerator, denominator);
    return {numerator / divisor, denominator / divisor};
}

u64 floor_scaled(u64 numerator, u64 denominator) {
    return u64((u128(numerator) * SCALE) / denominator);
}

u64 ceil_scaled(u64 numerator, u64 denominator) {
    u128 value = u128(numerator) * SCALE;
    return u64((value + denominator - 1) / denominator);
}

struct Context {
    int L;
    int cell;
    int e;
    int f;
};

struct Group {
    int L;
    int e;
    int f;
    std::vector<int> cells;
};

struct Channel {
    int g;
    int P;
    int Q;
};

struct Witness {
    bool valid = false;
    Fraction value;
    int g = 0, P = 0, Q = 0, L = 0, cell = 0, e = 0, f = 0;
};

auto witness_key(const Witness& w) {
    return std::tuple(w.g, w.P, w.Q, w.L, w.cell, w.e, w.f);
}

bool witness_better(const Witness& candidate, const Witness& incumbent) {
    if (!incumbent.valid) return true;
    if (fraction_less(candidate.value, incumbent.value)) return true;
    if (fraction_less(incumbent.value, candidate.value)) return false;
    return witness_key(candidate) < witness_key(incumbent);
}

u64 splitmix64(u64 x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

void hash_field(u64& state, u64 value) {
    state ^= splitmix64(value + state);
    state *= 0x100000001b3ULL;
}

struct Result {
    u64 target_bound_rows = 0;
    u64 phase_a_exact_rows = 0;
    u64 phase_b_exact_rows = 0;
    u64 failures = 0;
    u64 phase_a_hash = 0xcbf29ce484222325ULL;
    u64 phase_b_hash = 0xcbf29ce484222325ULL;
    Witness minimum;
    Witness first_failure;
    std::vector<std::uint16_t> candidate_band;
};

void hash_exact(u64& state, const Channel& channel, const Group& group, int cell,
                const MassFraction& mass) {
    for (i64 value : {i64(channel.g), i64(channel.P), i64(channel.Q), i64(group.L),
                      i64(cell), i64(group.e), i64(group.f), mass.numerator, mass.denominator}) {
        hash_field(state, u64(value));
    }
}

struct Tables {
    std::vector<Group> groups;
    std::vector<Channel> channels;
    std::vector<u64> first_error;
    std::vector<std::vector<u64>> second_error;
    std::vector<std::size_t> second_offset;
    int max_P[4] = {0, 263, 131, 87};
    int max_Q[4] = {0, 2103, 1047, 695};

    std::size_t first_index(int g, int P, std::size_t group) const {
        return (std::size_t(g - 1) * 264 + P) * groups.size() + group;
    }
};

u64 gamma_term_scaled(int k, int g, int L, int endpoint) {
    u64 numerator, denominator;
    if (k % 2 == 0) {
        numerator = u64(endpoint) * k;
        denominator = u64(2) * (u64(g) * L * k - endpoint);
    } else {
        numerator = u64(endpoint) * (u64(k) * k + 1);
        denominator = u64(2) * k * (u64(g) * L * k - endpoint);
    }
    return ceil_scaled(numerator, denominator);
}

Tables make_tables(std::vector<Group> groups, std::vector<Channel> channels) {
    Tables tables;
    tables.groups = std::move(groups);
    tables.channels = std::move(channels);
    const std::size_t G = tables.groups.size();
    tables.first_error.resize(std::size_t(3) * 264 * G);
    for (int g = 1; g <= 3; ++g) {
        for (int P = 1; P <= tables.max_P[g]; ++P) {
            for (std::size_t index = 0; index < G; ++index) {
                const Group& group = tables.groups[index];
                tables.first_error[tables.first_index(g, P, index)] =
                    gamma_term_scaled(P, g, group.L, group.e);
            }
        }
    }
    tables.second_error.resize(4);
    for (int g = 1; g <= 3; ++g) {
        tables.second_error[g].resize(std::size_t(tables.max_Q[g] + 1) * G);
        for (int Q = 1; Q <= tables.max_Q[g]; ++Q) {
            for (std::size_t index = 0; index < G; ++index) {
                const Group& group = tables.groups[index];
                tables.second_error[g][std::size_t(Q) * G + index] =
                    gamma_term_scaled(Q, g, group.L, group.f);
            }
        }
    }
    return tables;
}

i64 midpoint_lower_units(const Tables& tables, const Channel& channel, std::size_t group_index) {
    const Group& group = tables.groups[group_index];
    u64 product = u64(channel.P) * channel.Q;
    u64 phase;
    if (product <= 22) {
        phase = floor_scaled(1, 105);
    } else {
        phase = floor_scaled(product - 12, 49 * product);
    }
    u64 first = tables.first_error[tables.first_index(channel.g, channel.P, group_index)];
    u64 second = tables.second_error[channel.g][std::size_t(channel.Q) * tables.groups.size() + group_index];
    u64 determinant = u64(std::abs(channel.Q * group.e - channel.P * group.f));
    u64 third_numerator = determinant * (determinant / group.L + 1);
    u64 third_denominator = u64(2) * channel.g * channel.g * group.L * channel.P * channel.Q;
    u64 third = ceil_scaled(third_numerator, third_denominator);
    i128 answer = i128(phase) - first - second - third;
    return checked_i64(answer, "fixed-point midpoint lower bound");
}

bool lower_at_least(i64 lower_units, const Fraction& threshold) {
    if (lower_units < 0 || threshold.numerator < 0) return false;
    return i128(lower_units) * threshold.denominator >= i128(threshold.numerator) * SCALE;
}

void evaluate_group(const Channel& channel, const Group& group, Result& result, bool phase_b) {
    u64& exact_rows = phase_b ? result.phase_b_exact_rows : result.phase_a_exact_rows;
    u64& exact_hash = phase_b ? result.phase_b_hash : result.phase_a_hash;
    for (int cell : group.cells) {
        MassFraction raw = physical_mass(group.L, cell, group.e, channel.g * channel.P,
                                         group.f, channel.g * channel.Q);
        ++exact_rows;
        hash_exact(exact_hash, channel, group, cell, raw);
        Witness witness;
        witness.valid = true;
        witness.value = Fraction(raw.numerator, raw.denominator);
        witness.g = channel.g; witness.P = channel.P; witness.Q = channel.Q;
        witness.L = group.L; witness.cell = cell; witness.e = group.e; witness.f = group.f;
        if (witness_better(witness, result.minimum)) result.minimum = witness;
        if (i128(raw.numerator) * 294 < raw.denominator) {
            ++result.failures;
            if (!result.first_failure.valid) result.first_failure = witness;
        }
    }
}

template <class Function>
void parallel_channels(std::size_t count, int threads, Function function) {
    std::atomic<std::size_t> next{0};
    std::vector<std::thread> workers;
    for (int worker = 0; worker < threads; ++worker) {
        workers.emplace_back([&]() {
            while (true) {
                std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
                if (index >= count) return;
                function(index);
            }
        });
    }
    for (auto& worker : workers) worker.join();
}

std::string hex64(u64 value) {
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

std::string witness_text(const Witness& witness) {
    if (!witness.valid) return "none";
    std::ostringstream out;
    out << fraction_text(witness.value) << "@(g,P,Q;L,cell,e,f)=("
        << witness.g << ',' << witness.P << ',' << witness.Q << ';'
        << witness.L << ',' << witness.cell << ',' << witness.e << ',' << witness.f << ')';
    return out.str();
}

std::vector<Context> read_contexts(const std::string& path) {
    std::ifstream input(path);
    require(bool(input), "cannot open context atlas");
    std::vector<Context> contexts;
    Context row;
    while (input >> row.L >> row.cell >> row.e >> row.f) contexts.push_back(row);
    require(input.eof(), "malformed context atlas");
    require(contexts.size() == 2530, "unexpected context count");
    require(std::is_sorted(contexts.begin(), contexts.end(), [](const Context& a, const Context& b) {
        return std::tie(a.L, a.cell, a.e, a.f) < std::tie(b.L, b.cell, b.e, b.f);
    }), "context atlas is not sorted");
    require(std::adjacent_find(contexts.begin(), contexts.end(), [](const Context& a, const Context& b) {
        return std::tie(a.L, a.cell, a.e, a.f) == std::tie(b.L, b.cell, b.e, b.f);
    }) == contexts.end(), "duplicate context");
    for (const Context& context : contexts) {
        require(context.L < 4592 && context.L >= 168 && context.L % 14 == 0,
                "context outside small-ruler universe");
        int unit = context.L / 14;
        int r = int(positive_mod(context.e * context.cell, context.L));
        int s = int(positive_mod(context.f * context.cell, context.L));
        require(unit <= r && r + context.e <= context.L - unit, "unsafe e endpoint");
        require(unit <= s && s + context.f <= context.L - unit, "unsafe f endpoint");
    }
    return contexts;
}

std::vector<Group> make_groups(const std::vector<Context>& contexts) {
    std::map<std::tuple<int, int, int>, std::vector<int>> grouped;
    for (const Context& context : contexts) {
        grouped[{context.L, context.e, context.f}].push_back(context.cell);
    }
    std::vector<Group> groups;
    for (auto& [key, cells] : grouped) {
        std::sort(cells.begin(), cells.end());
        cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
        auto [L, e, f] = key;
        groups.push_back({L, e, f, std::move(cells)});
    }
    require(groups.size() == 1304, "unexpected grouped-context count");
    return groups;
}

std::vector<Channel> make_channels() {
    std::vector<Channel> channels;
    std::size_t counts[4] = {};
    for (int g = 1; g <= 3; ++g) {
        for (int P = 1; g * P < 264; ++P) {
            for (int Q = P + 1; Q < 8 * P; ++Q) {
                if (std::gcd(P, Q) != 1 || P + Q < 8 || (P == 3 && Q == 5)) continue;
                channels.push_back({g, P, Q});
                ++counts[g];
            }
        }
    }
    require(counts[1] == 148110 && counts[2] == 36978 && counts[3] == 16286,
            "unexpected channel counts");
    require(channels.size() == 201374, "unexpected total channel count");
    return channels;
}

int run_probe(const std::string& path) {
    std::ifstream input(path);
    require(bool(input), "cannot open probe input");
    i64 L, cell, e, p, f, q;
    while (input >> L >> cell >> e >> p >> f >> q) {
        MassFraction value = physical_mass(L, cell, e, p, f, q);
        std::cout << L << ' ' << cell << ' ' << e << ' ' << p << ' ' << f << ' ' << q
                  << ' ' << value.numerator << ' ' << value.denominator << '\n';
    }
    require(input.eof(), "malformed probe input");
    return 0;
}

int run_scan(const std::string& context_path, const std::string& summary_path, int thread_count) {
    std::vector<Context> contexts = read_contexts(context_path);
    std::set<int> rulers;
    for (const Context& context : contexts) rulers.insert(context.L);
    require(rulers.size() == 29, "unexpected ruler count");
    std::vector<Group> groups = make_groups(contexts);
    std::vector<Channel> channels = make_channels();
    Tables tables = make_tables(std::move(groups), std::move(channels));

    // This is an actual included row, not an external lower bound.  It caps
    // the only band that can contain the global minimum.
    const Fraction seed(92, 7645);
    MassFraction seed_check = physical_mass(168, 90, 12, 4, 6, 5);
    require(seed_check.numerator == seed.numerator && seed_check.denominator == seed.denominator,
            "candidate seed control mismatch");
    const Fraction target(1, 294);

    std::vector<Result> results(tables.channels.size());
    parallel_channels(tables.channels.size(), thread_count, [&](std::size_t channel_index) {
        const Channel& channel = tables.channels[channel_index];
        Result& result = results[channel_index];
        result.candidate_band.reserve(64);
        for (std::size_t group_index = 0; group_index < tables.groups.size(); ++group_index) {
            i64 lower = midpoint_lower_units(tables, channel, group_index);
            const Group& group = tables.groups[group_index];
            if (lower_at_least(lower, target)) {
                result.target_bound_rows += group.cells.size();
                if (!lower_at_least(lower, seed)) {
                    result.candidate_band.push_back(std::uint16_t(group_index));
                }
            } else {
                evaluate_group(channel, group, result, false);
            }
        }
    });

    Witness phase_a_minimum;
    u64 phase_a_failures = 0;
    for (const Result& result : results) {
        phase_a_failures += result.failures;
        if (result.minimum.valid && witness_better(result.minimum, phase_a_minimum)) {
            phase_a_minimum = result.minimum;
        }
    }
    require(phase_a_minimum.valid, "phase A evaluated no physical row");

    // Every row outside candidate_band has midpoint lower bound >=seed, and
    // phase_a_minimum<=seed.  Recheck only the target-safe part of that band
    // against the exact phase-A candidate.  If it produces a smaller value,
    // previously skipped rows remain safe because the candidate only fell.
    if (fraction_less(target, phase_a_minimum.value)) {
        parallel_channels(tables.channels.size(), thread_count, [&](std::size_t channel_index) {
            const Channel& channel = tables.channels[channel_index];
            Result& result = results[channel_index];
            for (std::uint16_t packed_index : result.candidate_band) {
                std::size_t group_index = packed_index;
                i64 lower = midpoint_lower_units(tables, channel, group_index);
                require(lower_at_least(lower, target), "candidate band lost target certificate");
                if (!lower_at_least(lower, phase_a_minimum.value)) {
                    evaluate_group(channel, tables.groups[group_index], result, true);
                }
            }
        });
    }

    Witness global_minimum;
    Witness first_failure;
    u64 target_bound_rows = 0, phase_a_exact_rows = 0, phase_b_exact_rows = 0;
    u64 candidate_band_groups = 0, failures = 0;
    for (const Result& result : results) {
        target_bound_rows += result.target_bound_rows;
        phase_a_exact_rows += result.phase_a_exact_rows;
        phase_b_exact_rows += result.phase_b_exact_rows;
        candidate_band_groups += result.candidate_band.size();
        failures += result.failures;
        if (result.minimum.valid && witness_better(result.minimum, global_minimum)) global_minimum = result.minimum;
        if (result.first_failure.valid && (!first_failure.valid || witness_key(result.first_failure) < witness_key(first_failure))) {
            first_failure = result.first_failure;
        }
    }
    const u64 total_rows = 509476220ULL;
    require(target_bound_rows + phase_a_exact_rows == total_rows, "phase-A row partition mismatch");
    require(global_minimum.valid, "missing global minimum");
    require(!fraction_less(phase_a_minimum.value, global_minimum.value), "phase-B minimum ordering failure");

    std::ofstream summary(summary_path, std::ios::binary);
    require(bool(summary), "cannot create semantic summary");
    for (std::size_t index = 0; index < tables.channels.size(); ++index) {
        const Channel& channel = tables.channels[index];
        const Result& result = results[index];
        summary << channel.g << ' ' << channel.P << ' ' << channel.Q << ' '
                << result.target_bound_rows << ' ' << result.phase_a_exact_rows << ' '
                << result.candidate_band.size() << ' ' << result.phase_b_exact_rows << ' '
                << result.failures << ' ' << hex64(result.phase_a_hash) << ' '
                << hex64(result.phase_b_hash) << ' ' << witness_text(result.minimum) << '\n';
    }
    require(bool(summary), "failed writing semantic summary");

    i128 gap_numerator = i128(global_minimum.value.numerator) * 294 - global_minimum.value.denominator;
    i128 gap_denominator = i128(global_minimum.value.denominator) * 294;
    i64 gap_gcd = std::gcd(checked_i64(gap_numerator, "gap numerator"),
                           checked_i64(gap_denominator, "gap denominator"));
    Fraction gap(checked_i64(gap_numerator / gap_gcd, "reduced gap numerator"),
                 checked_i64(gap_denominator / gap_gcd, "reduced gap denominator"));

    std::cout << "LRC14 DISCONNECTED-LOW FINITE HEAD EXACT CERTIFICATE\n";
    std::cout << "status=FINITE-EXACT physical pair-floor census; no affine-tail claim\n";
    std::cout << "contexts=" << contexts.size() << ";groups=" << tables.groups.size()
              << ";rulers=" << rulers.size() << "\n";
    std::cout << "channels_by_g=((1,148110),(2,36978),(3,16286));channels_total="
              << tables.channels.size() << "\n";
    std::cout << "physical_rows=" << total_rows << ";dyadic_screen_bits=" << SCALE_BITS << "\n";
    std::cout << "target_bound_rows=" << target_bound_rows
              << ";phase_a_exact_rows=" << phase_a_exact_rows
              << ";candidate_band_groups=" << candidate_band_groups
              << ";phase_b_exact_rows=" << phase_b_exact_rows << "\n";
    std::cout << "phase_a_minimum=" << witness_text(phase_a_minimum) << "\n";
    std::cout << "global_minimum=" << witness_text(global_minimum) << "\n";
    std::cout << "global_minimum_minus_1/294=" << fraction_text(gap) << "\n";
    std::cout << "failures_below_1/294=" << failures
              << ";first_failure=" << witness_text(first_failure) << "\n";
    std::cout << "seed_control=92/7645@(1,4,5;168,90,12,6)"
              << ";phase_a_failures=" << phase_a_failures << "\n";
    std::cout << "summary_rows=" << tables.channels.size() << "\n";
    std::cout << "conclusion=every finite-head physical mass is at least 1/294 iff failures_below_1/294=0\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--probe") return run_probe(argv[2]);
        if (argc < 3 || argc > 4) {
            std::cerr << "usage: " << argv[0] << " CONTEXTS SUMMARY [THREADS]\n"
                      << "   or: " << argv[0] << " --probe PROBES\n";
            return 2;
        }
        int threads = argc == 4 ? std::stoi(argv[3]) : int(std::thread::hardware_concurrency());
        if (threads <= 0) threads = 1;
        return run_scan(argv[1], argv[2], threads);
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
