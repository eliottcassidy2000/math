// Independent exact audit for THM-4233.
//
// This program deliberately does not use the primary common-grid/event sweep.
// It constructs each runner's safe teeth as reduced rational intervals,
// intersects the two ordered interval families directly, and integrates the
// centered indicator primitive with exact rational arithmetic.

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i128 = __int128_t;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

i64 checked_i64(i128 value, const char* where) {
    if (value < static_cast<i128>(std::numeric_limits<i64>::min()) ||
        value > static_cast<i128>(std::numeric_limits<i64>::max())) {
        fail(std::string("signed 64-bit overflow in ") + where);
    }
    return static_cast<i64>(value);
}

i64 magnitude(i64 value) {
    require(value != std::numeric_limits<i64>::min(), "unrepresentable magnitude");
    return value < 0 ? -value : value;
}

struct Rational {
    i64 num = 0;
    i64 den = 1;

    Rational() = default;

    Rational(i64 numerator, i64 denominator = 1) : num(numerator), den(denominator) {
        require(den != 0, "zero rational denominator");
        if (den < 0) {
            den = -den;
            num = -num;
        }
        const i64 divisor = std::gcd(magnitude(num), den);
        num /= divisor;
        den /= divisor;
    }

    std::string str() const {
        if (den == 1) return std::to_string(num);
        return std::to_string(num) + "/" + std::to_string(den);
    }
};

Rational operator-(const Rational& value) {
    return Rational(-value.num, value.den);
}

Rational operator+(const Rational& left, const Rational& right) {
    const i64 common = std::gcd(left.den, right.den);
    const i64 left_scale = right.den / common;
    const i64 right_scale = left.den / common;
    const i128 numerator = static_cast<i128>(left.num) * left_scale +
                           static_cast<i128>(right.num) * right_scale;
    const i128 denominator = static_cast<i128>(left.den) * left_scale;
    return Rational(checked_i64(numerator, "rational addition numerator"),
                    checked_i64(denominator, "rational addition denominator"));
}

Rational operator-(const Rational& left, const Rational& right) {
    return left + (-right);
}

Rational operator*(const Rational& left, const Rational& right) {
    const i64 cross_left = std::gcd(magnitude(left.num), right.den);
    const i64 cross_right = std::gcd(magnitude(right.num), left.den);
    const i128 numerator = static_cast<i128>(left.num / cross_left) *
                           (right.num / cross_right);
    const i128 denominator = static_cast<i128>(left.den / cross_right) *
                             (right.den / cross_left);
    return Rational(checked_i64(numerator, "rational multiplication numerator"),
                    checked_i64(denominator, "rational multiplication denominator"));
}

bool operator==(const Rational& left, const Rational& right) {
    return left.num == right.num && left.den == right.den;
}

bool operator<(const Rational& left, const Rational& right) {
    return static_cast<i128>(left.num) * right.den <
           static_cast<i128>(right.num) * left.den;
}

bool operator<=(const Rational& left, const Rational& right) {
    return !(right < left);
}

Rational rat_min(const Rational& left, const Rational& right) {
    return right < left ? right : left;
}

Rational rat_max(const Rational& left, const Rational& right) {
    return left < right ? right : left;
}

struct Interval {
    Rational lo;
    Rational hi;
};

bool operator==(const Interval& left, const Interval& right) {
    return left.lo == right.lo && left.hi == right.hi;
}

// The safe set G_m = {t in [0,1] : ||m t|| >= 1/14} consists, up to
// null endpoints, of these m disjoint teeth-complements.
std::vector<Interval> safe_teeth(i64 speed) {
    require(speed > 0, "speed must be positive");
    std::vector<Interval> intervals;
    intervals.reserve(static_cast<std::size_t>(speed));
    const i64 denominator = checked_i64(static_cast<i128>(14) * speed,
                                        "safe tooth denominator");
    for (i64 index = 0; index < speed; ++index) {
        intervals.push_back({Rational(14 * index + 1, denominator),
                             Rational(14 * index + 13, denominator)});
    }
    return intervals;
}

// Ordered two-pointer interval intersection.  This is independent of any
// common denominator or midpoint classification.
std::vector<Interval> intersect_ordered(const std::vector<Interval>& left,
                                        const std::vector<Interval>& right) {
    std::vector<Interval> answer;
    answer.reserve(left.size() + right.size());
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < left.size() && j < right.size()) {
        const Rational lo = rat_max(left[i].lo, right[j].lo);
        const Rational hi = rat_min(left[i].hi, right[j].hi);
        if (lo < hi) {
            if (!answer.empty() && !(answer.back().hi < lo)) {
                answer.back().hi = rat_max(answer.back().hi, hi);
            } else {
                answer.push_back({lo, hi});
            }
        }

        if (left[i].hi < right[j].hi) {
            ++i;
        } else if (right[j].hi < left[i].hi) {
            ++j;
        } else {
            ++i;
            ++j;
        }
    }
    return answer;
}

std::vector<Rational> boundary_points(const std::vector<Interval>& intervals) {
    std::vector<Rational> answer;
    answer.reserve(2 * intervals.size());
    for (const Interval& interval : intervals) {
        answer.push_back(interval.lo);
        answer.push_back(interval.hi);
    }
    return answer;
}

struct Audit {
    i64 first = 0;
    i64 second = 0;
    std::vector<Interval> components;
    Rational beta;
    Rational minimum;
    Rational minimum_at;
    Rational maximum;
    Rational maximum_at;
    Rational omega;
    std::vector<Rational> common_boundaries;
};

Audit audit_pair(i64 first, i64 second) {
    const std::vector<Interval> first_teeth = safe_teeth(first);
    const std::vector<Interval> second_teeth = safe_teeth(second);
    std::vector<Interval> components = intersect_ordered(first_teeth, second_teeth);
    require(!components.empty(), "empty safe intersection");

    Rational beta;
    for (std::size_t index = 0; index < components.size(); ++index) {
        const Interval& interval = components[index];
        require(interval.lo < interval.hi, "nonpositive component length");
        if (index != 0) {
            require(components[index - 1].hi < interval.lo,
                    "intersection components were not maximally merged");
        }
        beta = beta + (interval.hi - interval.lo);
    }

    // H(t) = integral_0^t (1_{G_first intersect G_second} - beta) ds.
    // H is affine between component endpoints, so all extrema occur at 0, 1,
    // or at one of the rational endpoints below.
    Rational covered;
    Rational minimum;
    Rational maximum;
    Rational minimum_at;
    Rational maximum_at;
    for (const Interval& interval : components) {
        const Rational at_left = covered - beta * interval.lo;
        if (at_left < minimum) {
            minimum = at_left;
            minimum_at = interval.lo;
        }
        if (maximum < at_left) {
            maximum = at_left;
            maximum_at = interval.lo;
        }

        covered = covered + (interval.hi - interval.lo);
        const Rational at_right = covered - beta * interval.hi;
        if (at_right < minimum) {
            minimum = at_right;
            minimum_at = interval.hi;
        }
        if (maximum < at_right) {
            maximum = at_right;
            maximum_at = interval.hi;
        }
    }
    require(covered == beta, "component total disagrees with beta");
    require(covered - beta == Rational(0), "centered primitive does not close");

    const std::vector<Rational> first_boundaries = boundary_points(first_teeth);
    const std::vector<Rational> second_boundaries = boundary_points(second_teeth);
    std::vector<Rational> common_boundaries;
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < first_boundaries.size() && j < second_boundaries.size()) {
        if (first_boundaries[i] < second_boundaries[j]) {
            ++i;
        } else if (second_boundaries[j] < first_boundaries[i]) {
            ++j;
        } else {
            common_boundaries.push_back(first_boundaries[i]);
            ++i;
            ++j;
        }
    }

    return {first,
            second,
            std::move(components),
            beta,
            minimum,
            minimum_at,
            maximum,
            maximum_at,
            maximum - minimum,
            std::move(common_boundaries)};
}

struct Fnv1a64 {
    u64 state = 14695981039346656037ULL;

    void word(u64 value) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (value >> (8 * byte)) & 0xffULL;
            state *= 1099511628211ULL;
        }
    }

    void rational(const Rational& value) {
        word(static_cast<u64>(value.num));
        word(static_cast<u64>(value.den));
    }
};

u64 semantic_ledger(const Audit& target,
                    const Audit& hostile,
                    const Audit& small,
                    const Audit& fibonacci) {
    Fnv1a64 hash;
    hash.word(static_cast<u64>(target.first));
    hash.word(static_cast<u64>(target.second));
    hash.word(static_cast<u64>(target.components.size()));
    for (const Interval& interval : target.components) {
        hash.rational(interval.lo);
        hash.rational(interval.hi);
    }
    hash.rational(target.beta);
    hash.rational(target.minimum);
    hash.rational(target.minimum_at);
    hash.rational(target.maximum);
    hash.rational(target.maximum_at);
    hash.rational(target.omega);
    for (const Rational& collision : target.common_boundaries) {
        hash.rational(collision);
    }
    hash.rational(hostile.beta);
    hash.rational(hostile.omega);
    hash.rational(small.beta);
    hash.rational(small.omega);
    hash.rational(fibonacci.beta);
    hash.rational(fibonacci.omega);
    return hash.state;
}

std::string hex64(u64 value) {
    std::ostringstream stream;
    stream << std::hex << std::setfill('0') << std::setw(16) << value;
    return stream.str();
}

i64 tick_on_grid(const Rational& point, i64 grid) {
    require(grid % point.den == 0, "rational endpoint is off the asserted grid");
    return checked_i64(static_cast<i128>(point.num) * (grid / point.den),
                       "endpoint grid tick");
}

}  // namespace

int main() {
    try {
        constexpr i64 u = 3713;
        constexpr i64 v = 5149;
        constexpr i64 pool_maximum = 290;
        const i64 grid = checked_i64(static_cast<i128>(14) * u * v, "target grid");
        const i64 raw_primitive_denominator =
            checked_i64(static_cast<i128>(grid) * grid, "raw primitive denominator");

        const Audit target = audit_pair(u, v);
        const Audit reversed = audit_pair(v, u);
        const Audit hostile = audit_pair(1, 13);
        const Audit small = audit_pair(1, 2);
        const Audit fibonacci = audit_pair(2584, 4181);
        const Audit dilation_two = audit_pair(2 * u, 2 * v);

        const Rational baseline(66, 91);
        const Rational threshold(1650, 8281LL * 3467LL);
        const Rational beta_margin = target.beta - baseline;
        const Rational threshold_margin = threshold - target.omega;

        require(std::gcd(u, v) == 1, "target pair is not primitive");
        require(u > pool_maximum && v > pool_maximum, "target pair meets the old pool");
        require(grid == 267655318, "target common denominator");
        require(2 + 2 * u + 2 * v == 17726, "raw endpoint-entry count");
        require(target.common_boundaries.size() == 2, "boundary-collision count");
        require(2 + 2 * u + 2 * v -
                    static_cast<i64>(target.common_boundaries.size()) ==
                    17724,
                "unique endpoint count");
        require(tick_on_grid(target.common_boundaries[0], grid) == 95591185,
                "first collision tick");
        require(tick_on_grid(target.common_boundaries[1], grid) == 172064133,
                "second collision tick");
        require(target.components.size() == 7595, "safe-component count");
        require(17724 - 1 - static_cast<i64>(target.components.size()) == 10128,
                "unsafe-cell count");
        require(target.beta == Rational(98322360, 133827659), "target beta");
        require(target.beta * Rational(grid) == Rational(196644720),
                "target safe tick count");
        require(beta_margin == Rational(16387038, 1739759567),
                "target beta margin above 66/91");
        require(target.minimum == Rational(-138535899, 4823550313337),
                "target primitive minimum");
        require(target.minimum_at == Rational(55959, 72086),
                "target primitive minimum location");
        require(tick_on_grid(target.minimum_at, grid) == 207775767,
                "target primitive minimum grid tick");
        require(target.maximum == Rational(138535899, 4823550313337),
                "target primitive maximum");
        require(target.maximum_at == Rational(16127, 72086),
                "target primitive maximum location");
        require(tick_on_grid(target.maximum_at, grid) == 59879551,
                "target primitive maximum grid tick");
        require(target.omega == Rational(277071798, 4823550313337),
                "target primitive oscillation");
        require(target.minimum * Rational(raw_primitive_denominator) ==
                    Rational(-2057535171948LL),
                "target raw primitive minimum numerator");
        require(target.maximum * Rational(raw_primitive_denominator) ==
                    Rational(2057535171948LL),
                "target raw primitive maximum numerator");
        require(threshold == Rational(1650, 28710227), "THM-4228 threshold");
        require(threshold_margin ==
                    Rational(82934716896LL, 2826229070241355051LL),
                "target oscillation threshold margin");

        require(reversed.components == target.components,
                "swapping the speeds changed the safe intervals");
        require(reversed.beta == target.beta && reversed.omega == target.omega,
                "swapping the speeds changed the observable");
        require(hostile.beta == Rational(66, 91), "(1,13) beta control");
        require(hostile.omega == Rational(990, 8281), "(1,13) omega control");
        require(small.beta == Rational(11, 14), "(1,2) beta control");
        require(small.omega == Rational(11, 98), "(1,2) omega control");
        require(fibonacci.omega == Rational(42967355, 553336008694LL),
                "Fibonacci predecessor omega control");
        require(threshold < fibonacci.omega,
                "Fibonacci predecessor unexpectedly passes the transfer threshold");
        require(dilation_two.beta == target.beta, "dilation changed beta");
        require(dilation_two.omega == target.omega * Rational(1, 2),
                "dilation-two primitive scaling");

        const u64 ledger = semantic_ledger(target, hostile, small, fibonacci);
        require(ledger == 0x635148be446ddc28ULL,
                "independent semantic interval ledger");

        std::cout << "THM-4233 independent reduced-rational safe-teeth audit\n";
        std::cout << "pair=" << u << ',' << v << " gcd=" << std::gcd(u, v)
                  << " pool_max=" << pool_maximum
                  << " outside_pool=" << (u > pool_maximum && v > pool_maximum) << '\n';
        std::cout << "method=ordered_reduced_rational_interval_intersection"
                     " primitive=exact_rational_endpoint_integration"
                     " imports_primary=0\n";
        std::cout << "grid=" << grid
                  << " raw_endpoint_entries=" << 2 + 2 * u + 2 * v
                  << " unique_endpoints="
                  << 2 + 2 * u + 2 * v - target.common_boundaries.size()
                  << " collision_ticks="
                  << tick_on_grid(target.common_boundaries[0], grid) << ','
                  << tick_on_grid(target.common_boundaries[1], grid) << '\n';
        std::cout << "safe_components=" << target.components.size()
                  << " unsafe_cells="
                  << 17724 - 1 - static_cast<i64>(target.components.size())
                  << " safe_ticks=" << (target.beta * Rational(grid)).str()
                  << " beta=" << target.beta.str()
                  << " beta_minus_66_over_91=" << beta_margin.str() << '\n';
        std::cout << "primitive_min_raw="
                  << (target.minimum * Rational(raw_primitive_denominator)).str()
                  << " at_tick=" << tick_on_grid(target.minimum_at, grid)
                  << " primitive_max_raw="
                  << (target.maximum * Rational(raw_primitive_denominator)).str()
                  << " at_tick=" << tick_on_grid(target.maximum_at, grid)
                  << " raw_den=" << raw_primitive_denominator << '\n';
        std::cout << "primitive_min=" << target.minimum.str()
                  << " at=" << target.minimum_at.str()
                  << " primitive_max=" << target.maximum.str()
                  << " at=" << target.maximum_at.str()
                  << " omega=" << target.omega.str() << '\n';
        std::cout << "thm4228_threshold=" << threshold.str()
                  << " threshold_margin=" << threshold_margin.str() << '\n';
        std::cout << "semantic_interval_ledger_fnv1a64_le=" << hex64(ledger) << '\n';
        std::cout << "hostile_1_13_beta=" << hostile.beta.str()
                  << " hostile_1_13_omega=" << hostile.omega.str()
                  << " small_1_2_beta=" << small.beta.str()
                  << " small_1_2_omega=" << small.omega.str() << '\n';
        std::cout << "fibonacci_2584_4181_beta=" << fibonacci.beta.str()
                  << " fibonacci_2584_4181_omega=" << fibonacci.omega.str()
                  << " fibonacci_predecessor_pass=" << (fibonacci.omega <= threshold)
                  << '\n';
        std::cout << "dilation_g2_beta=" << dilation_two.beta.str()
                  << " dilation_g2_omega=" << dilation_two.omega.str()
                  << " expected_half=" << (target.omega * Rational(1, 2)).str() << '\n';
        std::cout << "checks=PASS reduced_teeth,ordered_intersection,null_collisions,"
                     "rational_primitive,swap_symmetry,controls,dilation_scaling\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
