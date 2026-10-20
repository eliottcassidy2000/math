// Independent exhaustive strong-induction verification. No producer imports.
// Compile: g++ -std=c++17 -O2 -Wall -Wextra -Wconversion -pedantic audit_signed_descent.cpp -o /tmp/audit_signed_descent
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using U = std::uint64_t;
constexpr U N = U{1} << 26;
#ifndef AUDIT_STEP_LIMIT
#define AUDIT_STEP_LIMIT 1000000
#endif
constexpr U STEP_LIMIT = AUDIT_STEP_LIMIT;

U step(U n, bool positive) {
    if (n == 0) throw std::runtime_error("zero magnitude entered");
    if ((n & 1U) == 0) return n / 2;
    // Check the unhalved numerator before performing multiplication/addition.
    // The negative check is slightly conservative and can only cause failure.
    constexpr U top = std::numeric_limits<U>::max();
    if (positive) {
        if (n > (top - 1) / 3) throw std::overflow_error("3*n+1 overflow");
        return (3*n + 1)/2;
    }
    if (n > top/3) throw std::overflow_error("3*n-1 overflow");
    return (3*n - 1)/2;
}

void verify_cycle(const std::vector<U>& cycle, bool positive) {
    if (cycle.empty()) throw std::runtime_error("empty cycle");
    auto copy = cycle;
    std::sort(copy.begin(), copy.end());
    if (std::adjacent_find(copy.begin(), copy.end()) != copy.end())
        throw std::runtime_error("nonprimitive known cycle");
    if (cycle[0] != copy[0]) throw std::runtime_error("cycle base is not minimum");
    for (std::size_t i=0; i<cycle.size(); ++i)
        if (step(cycle[i], positive) != cycle[(i+1)%cycle.size()])
            throw std::runtime_error("incorrect known cycle transition");
}

void run(bool positive) {
    const std::vector<std::vector<U>> cycles = positive
        ? std::vector<std::vector<U>>{{1,2}}
        : std::vector<std::vector<U>>{{1},{5,7,10},{17,25,37,55,82,41,61,91,136,68,34}};
    std::vector<U> bases;
    for (const auto& cycle : cycles) { verify_cycle(cycle, positive); bases.push_back(cycle[0]); }
    U base_count=0, descent_count=0, total_steps=0, max_steps=0, max_steps_start=0;
    U max_value=0, max_value_start=0;
    // Each non-base start is accepted iff a finite exact trajectory reaches x<start.
    // Thus previously accepted starts suffice by strong induction, not conjecture.
    for (U start=1; start<=N; ++start) {
        if (std::find(bases.begin(), bases.end(), start) != bases.end()) {
            ++base_count;
            continue;
        }
        U x=start, steps=0;
        do {
            if (steps >= STEP_LIMIT) throw std::runtime_error("step limit at start " + std::to_string(start));
            x = step(x, positive);
            ++steps;
            if (x > max_value) { max_value=x; max_value_start=start; }
        } while (x >= start);
        if (!(x>=1 && x<start)) throw std::runtime_error("invalid descent witness");
        if (total_steps > std::numeric_limits<U>::max() - steps) throw std::overflow_error("step sum overflow");
        total_steps += steps;
        ++descent_count;
        if (steps>max_steps) { max_steps=steps; max_steps_start=start; }
    }
    if (base_count + descent_count != N || base_count != bases.size())
        throw std::runtime_error("coverage count mismatch");
    std::cout << (positive ? "positive" : "negative_magnitude")
              << " N=" << N << " known_cycle_bases=" << base_count
              << " certified_descents=" << descent_count << " total_steps=" << total_steps
              << " max_steps=" << max_steps << " max_steps_start=" << max_steps_start
              << " max_value=" << max_value << " max_value_start=" << max_value_start
              << " step_limit=" << STEP_LIMIT << " PASS\n";
}
void self_tests() {
    for (bool sign : {false, true}) {
        for (U bad : {U{0}, std::numeric_limits<U>::max()}) {
            bool failed=false;
            try { (void)step(bad, sign); }
            catch (const std::exception&) { failed=true; }
            if (!failed) throw std::runtime_error("guard hostile was incorrectly accepted");
        }
    }
    if (step(U{6148914691236517203}, true) != U{9223372036854775805}
        || step(U{6148914691236517205}, false) != U{9223372036854775807})
        throw std::runtime_error("safe uint64 boundary control failed");
}
int main() {
    try {
        self_tests();
        run(true);
        run(false);
        std::cout << "BOTH_SIGNED_UNIVERSES_VERIFIED\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
