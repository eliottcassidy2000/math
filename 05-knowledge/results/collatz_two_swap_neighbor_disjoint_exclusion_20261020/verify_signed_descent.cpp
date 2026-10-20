#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Mathematical trajectories use checked unsigned 128-bit arithmetic.
// The optional trace checksum alone intentionally wraps modulo 2^64.
using U128 = unsigned __int128;
constexpr std::uint64_t LIMIT = UINT64_C(1) << 26;
constexpr std::uint64_t STEP_LIMIT = 1000000;
constexpr U128 MAX128 = ~U128(0);

std::string decimal(U128 value) {
    if (value == 0) return "0";
    std::string result;
    while (value) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

bool known(U128 x, int sign) {
    return x == 1 || (sign == -1 && (x == 5 || x == 17));
}

U128 odd_numerator(U128 x, int sign) {
    if (x == 0 || (x & 1) == 0)
        throw std::runtime_error("Odd-step input is not positive and odd");
    if (sign == 1) {
        if (x > (MAX128 - 1) / 3)
            throw std::runtime_error("Positive trajectory arithmetic overflow");
        return 3*x + 1;
    }
    if (sign != -1 || x > MAX128 / 3)
        throw std::runtime_error("Negative trajectory arithmetic overflow");
    return 3*x - 1;
}

U128 half_step(U128 x, int sign) {
    return (x & 1) ? odd_numerator(x, sign)/2 : x/2;
}

void verify_cycle(const std::vector<std::uint64_t>& cycle, int sign) {
    if (cycle.empty()) throw std::runtime_error("Empty claimed cycle");
    for (std::size_t i = 0; i < cycle.size(); ++i) {
        if (half_step(cycle[i], sign) != cycle[(i+1) % cycle.size()])
            throw std::runtime_error("Claimed cycle edge failed");
    }
}

struct Receipt {
    int sign;
    std::uint64_t odd_starts = 0;
    std::uint64_t exempt_known_minima = 0;
    std::uint64_t total_odd_steps = 0;
    std::uint64_t maximum_descent_odd_steps = 0;
    std::uint64_t maximum_step_start = 0;
    U128 maximum_prehalving_numerator = 0;
    std::uint64_t trace = UINT64_C(1469598103934665603);
};

Receipt verify_sign(int sign) {
    Receipt result;
    result.sign = sign;
    for (std::uint64_t start = 1; start <= LIMIT; start += 2) {
        if (known(start, sign)) {
            ++result.exempt_known_minima;
            continue;
        }
        ++result.odd_starts;
        U128 x = start;
        std::uint64_t steps = 0;
        while (x >= start && !known(x, sign)) {
            if (++steps > STEP_LIMIT)
                throw std::runtime_error("Descent step cap reached");
            U128 numerator = odd_numerator(x, sign);
            result.maximum_prehalving_numerator =
                std::max(result.maximum_prehalving_numerator, numerator);
            x = numerator;
            do { x /= 2; } while ((x & 1) == 0);
        }
        if (!(x > 0 && (x < start || known(x, sign))))
            throw std::runtime_error("No valid descent or known-cycle terminal");
        if (steps > std::numeric_limits<std::uint64_t>::max()
                        - result.total_odd_steps)
            throw std::runtime_error("Step aggregate overflow");
        result.total_odd_steps += steps;
        if (steps > result.maximum_descent_odd_steps) {
            result.maximum_descent_odd_steps = steps;
            result.maximum_step_start = start;
        }
        // Deterministic trace, not used for any mathematical decision.
        for (std::uint64_t value : {
                start, steps, static_cast<std::uint64_t>(x),
                static_cast<std::uint64_t>(x >> 64)}) {
            result.trace ^= value;
            result.trace *= UINT64_C(1099511628211);
        }
    }
    const std::uint64_t expected_exempt = sign == 1 ? 1 : 3;
    if (result.exempt_known_minima != expected_exempt
        || result.odd_starts + result.exempt_known_minima != LIMIT/2)
        throw std::runtime_error("Incomplete odd-start coverage");
    return result;
}

void print_receipt(const Receipt& r) {
    std::cout
        << "    {\"sign\": " << r.sign
        << ", \"tested_odd_starts\": " << r.odd_starts
        << ", \"known_minimum_starts\": " << r.exempt_known_minima
        << ", \"even_starts_covered_by_immediate_halving\": " << LIMIT/2
        << ", \"total_odd_descent_steps\": " << r.total_odd_steps
        << ", \"maximum_odd_descent_steps\": " << r.maximum_descent_odd_steps
        << ", \"first_maximum_step_start\": " << r.maximum_step_start
        << ", \"maximum_prehalving_numerator\": \""
        << decimal(r.maximum_prehalving_numerator)
        << "\", \"trace_mod_2_to_64\": \"" << r.trace << "\"}";
}

int main() {
    try {
        verify_cycle({1, 2}, 1);
        verify_cycle({1}, -1);
        verify_cycle({5, 7, 10}, -1);
        verify_cycle({17, 25, 37, 55, 82, 41, 61, 91, 136, 68, 34}, -1);
        const Receipt positive = verify_sign(1);
        const Receipt negative = verify_sign(-1);
        std::cout
            << "{\n  \"status\": \"PASS_EXACT_BOUNDED_SIGNED_DESCENT\",\n"
            << "  \"magnitude_limit_inclusive\": " << LIMIT << ",\n"
            << "  \"arithmetic\": \"checked unsigned 128-bit trajectories\",\n"
            << "  \"failure_policy\": \"nonzero exit on overflow, step cap, "
               "cycle error, or incomplete coverage\",\n"
            << "  \"per_start_odd_step_cap\": " << STEP_LIMIT << ",\n"
            << "  \"proof_method\": \"induction on the start; even starts halve "
               "immediately, every other odd start reaches a smaller value "
               "or an explicitly verified known cycle\",\n"
            << "  \"maps\": [\n";
        print_receipt(positive);
        std::cout << ",\n";
        print_receipt(negative);
        std::cout
            << "\n  ],\n"
            << "  \"known_positive_cycle\": [1, 2],\n"
            << "  \"known_negative_magnitude_cycles\": "
               "[[1], [5, 7, 10], "
               "[17, 25, 37, 55, 82, 41, 61, 91, 136, 68, 34]],\n"
            << "  \"scope\": \"All integer magnitudes at most the stated bound "
               "only; no claim of global Collatz convergence.\"\n}\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "CERTIFICATE FAILURE: " << e.what() << "\n";
        return 1;
    }
}
