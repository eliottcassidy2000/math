#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

namespace {

using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr int kMinM = 2;
constexpr int kMaxM = 1000;
constexpr int kDirectPairMaxM = 100;
constexpr int kDirectTailMax = 101;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

struct Rational {
    i64 numerator;
    i64 denominator;

    Rational(i64 n = 0, i64 d = 1) : numerator(n), denominator(d) {
        require(d != 0, "zero rational denominator");
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
        const i64 divisor = std::gcd(numerator < 0 ? -numerator : numerator,
                                     denominator);
        numerator /= divisor;
        denominator /= divisor;
    }
};

bool operator==(const Rational& left, const Rational& right) {
    return left.numerator == right.numerator
        && left.denominator == right.denominator;
}

bool operator<(const Rational& left, const Rational& right) {
    return left.numerator * right.denominator
         < right.numerator * left.denominator;
}

bool operator<=(const Rational& left, const Rational& right) {
    return !(right < left);
}

Rational operator+(const Rational& left, const Rational& right) {
    return Rational(left.numerator * right.denominator
                      + right.numerator * left.denominator,
                    left.denominator * right.denominator);
}

Rational operator-(const Rational& left, const Rational& right) {
    return Rational(left.numerator * right.denominator
                      - right.numerator * left.denominator,
                    left.denominator * right.denominator);
}

Rational operator*(const Rational& value, i64 scalar) {
    return Rational(value.numerator * scalar, value.denominator);
}

Rational operator/(const Rational& value, i64 scalar) {
    return Rational(value.numerator, value.denominator * scalar);
}

i64 floor_nonnegative(const Rational& value) {
    require(value.numerator >= 0, "negative floor input");
    return value.numerator / value.denominator;
}

bool circle_gap_safe(i64 speed, const Rational& phase) {
    const Rational scaled = phase * speed;
    i64 residue = scaled.numerator % scaled.denominator;
    if (residue < 0) {
        residue += scaled.denominator;
    }
    const i64 gap_numerator = std::min(residue,
                                       scaled.denominator - residue);
    return 14 * gap_numerator >= scaled.denominator;
}

bool speed_safe_on_interval(i64 speed,
                            const Rational& lower,
                            const Rational& upper) {
    const Rational scaled_lower = lower * speed;
    const Rational scaled_upper = upper * speed;
    const i64 first = std::max<i64>(0, floor_nonnegative(scaled_lower) - 1);
    const i64 last = floor_nonnegative(scaled_upper) + 1;
    for (i64 band = first; band <= last; ++band) {
        const Rational safe_lower(14 * band + 1, 14);
        const Rational safe_upper(14 * band + 13, 14);
        if (safe_lower <= scaled_lower && scaled_upper <= safe_upper) {
            return true;
        }
    }
    return false;
}

int block_lower(int m, int band) {
    return (14 * band + 1) * m;
}

int block_upper(int m, int band) {
    const int s = 12 * m + 1;
    return ((14 * band + 13) * s) / 16;
}

std::vector<int> proposed_pool(int m) {
    std::vector<int> pool;
    for (int band = 0; band <= 2; ++band) {
        const int lower = block_lower(m, band);
        const int upper = block_upper(m, band);
        require(lower <= upper, "empty proposed block");
        for (int speed = lower; speed <= upper; ++speed) {
            pool.push_back(speed);
        }
    }
    require(std::adjacent_find(pool.begin(), pool.end()) == pool.end(),
            "duplicate pool speed");
    require(std::is_sorted(pool.begin(), pool.end()),
            "pool is not sorted");
    return pool;
}

std::vector<int> proposed_pool_m1() {
    std::vector<int> pool;
    for (int band = 0; band <= 3; ++band) {
        const int lower = block_lower(1, band);
        const int upper = block_upper(1, band);
        require(lower <= upper, "empty m=1 proposed block");
        for (int speed = lower; speed <= upper; ++speed) {
            pool.push_back(speed);
        }
    }
    return pool;
}

bool in_proposed_pool_m1(int speed) {
    for (int band = 0; band <= 3; ++band) {
        if (block_lower(1, band) <= speed
            && speed <= block_upper(1, band)) {
            return true;
        }
    }
    return false;
}

bool in_proposed_pool(int speed, int m) {
    for (int band = 0; band <= 2; ++band) {
        if (block_lower(m, band) <= speed
            && speed <= block_upper(m, band)) {
            return true;
        }
    }
    return false;
}

u64 choose_u64(unsigned n, unsigned k) {
    if (k > n) {
        return 0;
    }
    k = std::min(k, n - k);
    u64 result = 1;
    for (unsigned index = 1; index <= k; ++index) {
        u64 numerator = static_cast<u64>(n - k + index);
        u64 denominator = static_cast<u64>(index);
        const u64 first_gcd = std::gcd(numerator, denominator);
        numerator /= first_gcd;
        denominator /= first_gcd;
        const u64 second_gcd = std::gcd(result, denominator);
        result /= second_gcd;
        denominator /= second_gcd;
        require(denominator == 1, "binomial cancellation failed");
        require(result <= std::numeric_limits<u64>::max() / numerator,
                "binomial overflow");
        result *= numerator;
    }
    return result;
}

class BigUInt {
public:
    explicit BigUInt(u64 value = 0) {
        do {
            digits_.push_back(static_cast<std::uint32_t>(value % kBase));
            value /= kBase;
        } while (value != 0);
    }

    bool is_zero() const {
        return digits_.size() == 1 && digits_[0] == 0;
    }

    int compare(const BigUInt& other) const {
        if (digits_.size() != other.digits_.size()) {
            return digits_.size() < other.digits_.size() ? -1 : 1;
        }
        for (std::size_t index = digits_.size(); index-- > 0;) {
            if (digits_[index] != other.digits_[index]) {
                return digits_[index] < other.digits_[index] ? -1 : 1;
            }
        }
        return 0;
    }

    void add(const BigUInt& other) {
        const std::size_t size = std::max(digits_.size(), other.digits_.size());
        digits_.resize(size, 0);
        u64 carry = 0;
        for (std::size_t index = 0; index < size; ++index) {
            const u64 right = index < other.digits_.size()
                            ? other.digits_[index] : 0;
            const u64 total = static_cast<u64>(digits_[index]) + right + carry;
            digits_[index] = static_cast<std::uint32_t>(total % kBase);
            carry = total / kBase;
        }
        if (carry != 0) {
            digits_.push_back(static_cast<std::uint32_t>(carry));
        }
    }

    void subtract(const BigUInt& other) {
        require(compare(other) >= 0, "negative BigUInt subtraction");
        i64 borrow = 0;
        for (std::size_t index = 0; index < digits_.size(); ++index) {
            const i64 right = index < other.digits_.size()
                            ? static_cast<i64>(other.digits_[index]) : 0;
            i64 value = static_cast<i64>(digits_[index]) - right - borrow;
            if (value < 0) {
                value += static_cast<i64>(kBase);
                borrow = 1;
            } else {
                borrow = 0;
            }
            digits_[index] = static_cast<std::uint32_t>(value);
        }
        require(borrow == 0, "BigUInt subtraction borrow drift");
        trim();
    }

    void multiply_small(std::uint32_t multiplier) {
        if (multiplier == 0 || is_zero()) {
            digits_.assign(1, 0);
            return;
        }
        u64 carry = 0;
        for (std::uint32_t& digit : digits_) {
            const u64 product = static_cast<u64>(digit) * multiplier + carry;
            digit = static_cast<std::uint32_t>(product % kBase);
            carry = product / kBase;
        }
        while (carry != 0) {
            digits_.push_back(static_cast<std::uint32_t>(carry % kBase));
            carry /= kBase;
        }
    }

    std::uint32_t divide_small(std::uint32_t divisor) {
        require(divisor != 0, "BigUInt division by zero");
        u64 remainder = 0;
        for (std::size_t index = digits_.size(); index-- > 0;) {
            const u64 current = remainder * kBase + digits_[index];
            digits_[index] = static_cast<std::uint32_t>(current / divisor);
            remainder = current % divisor;
        }
        trim();
        return static_cast<std::uint32_t>(remainder);
    }

    std::uint32_t modulo_small(std::uint32_t divisor) const {
        require(divisor != 0, "BigUInt modulus by zero");
        u64 remainder = 0;
        for (std::size_t index = digits_.size(); index-- > 0;) {
            remainder = (remainder * kBase + digits_[index]) % divisor;
        }
        return static_cast<std::uint32_t>(remainder);
    }

    std::string to_string() const {
        std::ostringstream stream;
        stream << digits_.back();
        for (std::size_t index = digits_.size() - 1; index-- > 0;) {
            stream << std::setw(9) << std::setfill('0') << digits_[index];
        }
        return stream.str();
    }

private:
    static constexpr u64 kBase = UINT64_C(1000000000);
    std::vector<std::uint32_t> digits_;

    void trim() {
        while (digits_.size() > 1 && digits_.back() == 0) {
            digits_.pop_back();
        }
    }
};

using FactorMap = std::map<std::uint32_t, int>;

void add_factors(FactorMap& factors, u64 value, int multiplicity) {
    require(value > 0, "cannot factor zero");
    for (u64 prime = 2; prime * prime <= value; ++prime) {
        while (value % prime == 0) {
            factors[static_cast<std::uint32_t>(prime)] += multiplicity;
            value /= prime;
        }
    }
    if (value > 1) {
        factors[static_cast<std::uint32_t>(value)] += multiplicity;
    }
}

BigUInt power_small(u64 base, int exponent) {
    BigUInt result(1);
    for (int count = 0; count < exponent; ++count) {
        require(base <= std::numeric_limits<std::uint32_t>::max(),
                "power base exceeds BigUInt small multiplier");
        result.multiply_small(static_cast<std::uint32_t>(base));
    }
    return result;
}

void multiply_factors(BigUInt& value, const FactorMap& factors) {
    for (const auto& entry : factors) {
        for (int count = 0; count < entry.second; ++count) {
            value.multiply_small(entry.first);
        }
    }
}

FactorMap factor_difference(const FactorMap& larger,
                            const FactorMap& smaller) {
    FactorMap difference;
    for (const auto& entry : larger) {
        const auto found = smaller.find(entry.first);
        const int exponent = entry.second
                           - (found == smaller.end() ? 0 : found->second);
        require(exponent >= 0, "negative factor-map difference");
        if (exponent > 0) {
            difference[entry.first] = exponent;
        }
    }
    return difference;
}

FactorMap factor_maximum(const std::vector<FactorMap>& maps) {
    FactorMap maximum;
    for (const FactorMap& factors : maps) {
        for (const auto& entry : factors) {
            maximum[entry.first] = std::max(maximum[entry.first], entry.second);
        }
    }
    return maximum;
}

struct BigFraction {
    BigUInt numerator;
    FactorMap denominator_factors;
};

void reduce_fraction(BigFraction& fraction) {
    for (auto& entry : fraction.denominator_factors) {
        while (entry.second > 0
               && fraction.numerator.modulo_small(entry.first) == 0) {
            require(fraction.numerator.divide_small(entry.first) == 0,
                    "BigFraction exact division drift");
            --entry.second;
        }
    }
}

std::string denominator_string(const FactorMap& factors) {
    BigUInt denominator(1);
    multiply_factors(denominator, factors);
    return denominator.to_string();
}

BigFraction sum_positive_fractions(const std::vector<BigFraction>& fractions) {
    std::vector<FactorMap> denominator_maps;
    denominator_maps.reserve(fractions.size());
    for (const BigFraction& fraction : fractions) {
        denominator_maps.push_back(fraction.denominator_factors);
    }
    const FactorMap common = factor_maximum(denominator_maps);
    BigUInt numerator(0);
    for (const BigFraction& fraction : fractions) {
        BigUInt scaled = fraction.numerator;
        multiply_factors(scaled,
                         factor_difference(common,
                                           fraction.denominator_factors));
        numerator.add(scaled);
    }
    BigFraction result{numerator, common};
    reduce_fraction(result);
    return result;
}

BigFraction subtract_fractions(const BigFraction& left,
                               const BigFraction& right) {
    const FactorMap common = factor_maximum(
        {left.denominator_factors, right.denominator_factors}
    );
    BigUInt left_scaled = left.numerator;
    BigUInt right_scaled = right.numerator;
    multiply_factors(left_scaled,
                     factor_difference(common, left.denominator_factors));
    multiply_factors(right_scaled,
                     factor_difference(common, right.denominator_factors));
    require(left_scaled.compare(right_scaled) > 0,
            "expected positive fraction difference");
    left_scaled.subtract(right_scaled);
    BigFraction result{left_scaled, common};
    reduce_fraction(result);
    return result;
}

std::pair<u64, u64> reduced_ratio(u64 numerator, u64 denominator) {
    const u64 divisor = std::gcd(numerator, denominator);
    return {numerator / divisor, denominator / divisor};
}

class SemanticLedger {
public:
    void add(const std::string& key, const std::string& value) {
        const std::string row = key + "=" + value;
        rows_.push_back(row);
        for (const unsigned char character : row + "\n") {
            hash_ ^= static_cast<u64>(character);
            hash_ *= UINT64_C(1099511628211);
        }
    }

    void add(const std::string& key, u64 value) {
        add(key, std::to_string(value));
    }

    const std::vector<std::string>& rows() const {
        return rows_;
    }

    u64 hash() const {
        return hash_;
    }

private:
    std::vector<std::string> rows_;
    u64 hash_ = UINT64_C(14695981039346656037);
};

struct PoolAudit {
    u64 speeds_tested = 0;
    u64 pool_speeds = 0;
    u64 outside_speeds = 0;
    i64 minimum_fourth_band_defect = std::numeric_limits<i64>::max();
    int minimum_defect_m = 0;
};

PoolAudit audit_complete_pools() {
    PoolAudit audit;
    for (int m = kMinM; m <= kMaxM; ++m) {
        const int s = 12 * m + 1;
        const Rational lower(1, 14 * m);
        const Rational upper(8, 7 * s);
        const Rational width = upper - lower;
        const Rational expected_width(4 * m - 1, 14 * m * s);
        require(width == expected_width,
                "width mismatch at m=" + std::to_string(m));

        const std::vector<int> pool = proposed_pool(m);
        const int final_upper = block_upper(m, 2);
        std::vector<int> reconstructed;
        for (int speed = 1; speed <= final_upper + 64; ++speed) {
            const bool safe = speed_safe_on_interval(speed, lower, upper);
            const bool expected = in_proposed_pool(speed, m);
            require(safe == expected,
                    "pool mismatch at m=" + std::to_string(m)
                    + ", speed=" + std::to_string(speed));
            if (safe) {
                reconstructed.push_back(speed);
                ++audit.pool_speeds;
            } else {
                ++audit.outside_speeds;
            }
            ++audit.speeds_tested;
        }
        require(reconstructed == pool,
                "reconstructed pool vector mismatch at m="
                + std::to_string(m));

        // For a hypothetical safe band k>=3, endpoint containment would
        // require this defect to be <=0.  It is already 1 at (m,k)=(2,3)
        // and increases in both variables.
        const i64 defect = static_cast<i64>(16) * (14 * 3 + 1) * m
                         - static_cast<i64>(14 * 3 + 13) * s;
        if (defect < audit.minimum_fourth_band_defect) {
            audit.minimum_fourth_band_defect = defect;
            audit.minimum_defect_m = m;
        }
        require(defect > 0, "fourth band not excluded");
        require(56 * m - 14 > 0,
                "higher-band defect is not increasing in band index");
    }
    require(audit.minimum_fourth_band_defect == 1
            && audit.minimum_defect_m == 2,
            "sharp fourth-band defect drift");
    return audit;
}

struct M1Audit {
    std::size_t pool_size = 0;
    u64 family_count = 0;
    i64 fifth_band_defect = 0;
};

M1Audit audit_m1_boundary() {
    const Rational lower(1, 14);
    const Rational upper(8, 91);
    const std::vector<int> pool = proposed_pool_m1();
    require(block_lower(1, 0) == 1 && block_upper(1, 0) == 10,
            "m=1 first block drift");
    require(block_lower(1, 1) == 15 && block_upper(1, 1) == 21,
            "m=1 second block drift");
    require(block_lower(1, 2) == 29 && block_upper(1, 2) == 33,
            "m=1 third block drift");
    require(block_lower(1, 3) == 43 && block_upper(1, 3) == 44,
            "m=1 fourth block drift");
    require(pool.size() == 24, "m=1 pool-size drift");
    for (int speed = 1; speed <= block_upper(1, 3) + 64; ++speed) {
        require(speed_safe_on_interval(speed, lower, upper)
                == in_proposed_pool_m1(speed),
                "m=1 pool reconstruction mismatch at speed="
                + std::to_string(speed));
    }
    const i64 fifth_defect = static_cast<i64>(16) * (14 * 4 + 1)
                           - static_cast<i64>(14 * 4 + 13) * 13;
    require(fifth_defect == 15, "m=1 fifth-band defect drift");
    require(56 - 14 > 0, "m=1 later-band monotonicity drift");
    const u64 family_count = choose_u64(23, 10);
    require(family_count == UINT64_C(1144066),
            "m=1 family count drift");
    return M1Audit{pool.size(), family_count, fifth_defect};
}

struct EndpointAudit {
    u64 small_tail_endpoint_gates = 0;
    u64 large_tail_room_gates = 0;
    u64 tooth_gates = 0;
};

EndpointAudit audit_endpoint_and_tooth() {
    EndpointAudit audit;
    for (int m = 1; m <= kMaxM; ++m) {
        const int s = 12 * m + 1;
        const Rational lower(1, 14 * m);
        const Rational upper(8, 7 * s);
        const Rational width(4 * m - 1, 14 * m * s);
        const Rational tooth_left(6, 7 * s);
        const Rational tooth_right(8, 7 * s);
        const Rational exact_left_slack(1, 14 * m * s);
        require(lower - tooth_left == exact_left_slack,
                "left tooth slack mismatch");
        require(upper == tooth_right,
                "right tooth wall mismatch");
        require(upper - lower == width,
                "tooth width mismatch");
        ++audit.tooth_gates;

        const Rational upper_lift = (lower + Rational(1)) / 2;
        for (int tail = 1; tail <= 12 * m - 1; tail += 2) {
            require(circle_gap_safe(tail, upper_lift),
                    "small-tail endpoint failure at m="
                    + std::to_string(m));
            ++audit.small_tail_endpoint_gates;
        }

        for (int tail : {s, s + 2, 3 * s}) {
            const Rational remaining_room(4 * m - 1, 14 * m * tail);
            require(remaining_room <= width,
                    "large-tail tooth room failure");
            ++audit.large_tail_room_gates;
        }

        // At the sharp b=s endpoint, at least one physical lift is exactly
        // safe.  The second lift is the equality-wall control.
        const Rational lower_physical = upper / 2;
        const Rational upper_physical = (upper + Rational(1)) / 2;
        require(circle_gap_safe(s, lower_physical)
                && circle_gap_safe(s, upper_physical),
                "sharp tail endpoint is not closed-safe");
    }
    return audit;
}

void add_carrier_walls(std::vector<Rational>& walls,
                       int tail,
                       const Rational& lower,
                       const Rational& upper) {
    const i64 maximum_index = floor_nonnegative(upper * tail) + 2;
    for (i64 index = 0; index <= maximum_index; ++index) {
        const i64 plus_numerator = 7 * index + 1;
        const Rational plus_wall(plus_numerator, 7 * tail);
        if (lower <= plus_wall && plus_wall <= upper) {
            walls.push_back(plus_wall);
        }
        const i64 minus_numerator = 7 * index - 1;
        if (minus_numerator > 0) {
            const Rational minus_wall(minus_numerator, 7 * tail);
            if (lower <= minus_wall && minus_wall <= upper) {
                walls.push_back(minus_wall);
            }
        }
    }
}

bool direct_tail_pair_witness(int first,
                              int second,
                              const Rational& lower,
                              const Rational& upper,
                              u64& phase_gates) {
    std::vector<Rational> walls{lower, upper};
    add_carrier_walls(walls, first, lower, upper);
    add_carrier_walls(walls, second, lower, upper);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::vector<Rational> candidates = walls;
    for (std::size_t index = 0; index + 1 < walls.size(); ++index) {
        candidates.push_back((walls[index] + walls[index + 1]) / 2);
    }
    for (const Rational& y : candidates) {
        const Rational first_lift = y / 2;
        const Rational second_lift = (y + Rational(1)) / 2;
        phase_gates += 2;
        const bool first_safe = circle_gap_safe(first, first_lift)
                             && circle_gap_safe(second, first_lift);
        const bool second_safe = circle_gap_safe(first, second_lift)
                              && circle_gap_safe(second, second_lift);
        if (first_safe || second_safe) {
            return true;
        }
    }
    return false;
}

struct PairAudit {
    u64 pairs = 0;
    u64 phase_gates = 0;
};

PairAudit audit_direct_tail_pairs() {
    PairAudit audit;
    for (int m = 1; m <= kDirectPairMaxM; ++m) {
        const Rational lower(1, 14 * m);
        const Rational upper(8, 7 * (12 * m + 1));
        for (int first = 1; first <= kDirectTailMax; first += 2) {
            for (int second = first + 2;
                 second <= kDirectTailMax;
                 second += 2) {
                require(direct_tail_pair_witness(first, second, lower, upper,
                                                 audit.phase_gates),
                        "direct tail-pair failure at m=" + std::to_string(m)
                        + ", tails=" + std::to_string(first) + ","
                        + std::to_string(second));
                ++audit.pairs;
            }
        }
    }
    return audit;
}

struct M7Audit {
    std::size_t pool_size = 0;
    std::size_t free_size = 0;
    u64 family_count = 0;
    u64 all_m7_eleven_bodies = 0;
    std::pair<u64, u64> fraction_in_pool;
    std::pair<u64, u64> fraction_in_ambient;
};

M7Audit audit_m7_family() {
    constexpr int m = 7;
    const std::vector<int> pool = proposed_pool(m);
    const std::vector<int> anchors{7, 120, 126, 143};
    require(block_lower(m, 0) == 7 && block_upper(m, 0) == 69,
            "m=7 first block drift");
    require(block_lower(m, 1) == 105 && block_upper(m, 1) == 143,
            "m=7 second block drift");
    require(block_lower(m, 2) == 203 && block_upper(m, 2) == 217,
            "m=7 third block drift");
    require(pool.size() == 117, "m=7 pool-size drift");
    for (const int anchor : anchors) {
        require(std::binary_search(pool.begin(), pool.end(), anchor),
                "m=7 anchor outside pool");
    }

    int anchor_gcd = 0;
    for (const int anchor : anchors) {
        anchor_gcd = std::gcd(anchor_gcd, anchor);
    }
    require(anchor_gcd == 1, "anchor gcd is not primitive");
    for (int divisor = 1; divisor <= 14; ++divisor) {
        bool covered = false;
        for (const int anchor : anchors) {
            if (anchor % divisor == 0) {
                covered = true;
                break;
            }
        }
        require(covered,
                "anchor divisor coverage failure at d="
                + std::to_string(divisor));
    }

    const int forced_minimum = 7;
    const int forced_maximum = 143;
    require(27 * (13 * forced_minimum - forced_maximum)
            < 4 * forced_minimum * forced_maximum,
            "THM-4148 gate unexpectedly passes");
    require(16 * forced_maximum > 156 * forced_minimum + 13,
            "THM-4151 gate unexpectedly passes");
    require(!circle_gap_safe(2 * 120, Rational(1, 12)),
            "anchor 120 does not kill x=1/12");

    M7Audit audit;
    audit.pool_size = pool.size();
    audit.free_size = pool.size() - anchors.size();
    audit.family_count = choose_u64(static_cast<unsigned>(audit.free_size), 7);
    require(audit.family_count == UINT64_C(38620298376),
            "m=7 family count drift");
    audit.all_m7_eleven_bodies = choose_u64(116, 10);
    audit.fraction_in_pool = reduced_ratio(audit.family_count,
                                           choose_u64(117, 11));
    audit.fraction_in_ambient = reduced_ratio(audit.family_count,
                                              choose_u64(217, 11));
    return audit;
}

int available_pool_count(int m, int ambient_maximum) {
    int count = 0;
    for (int band = 0; band <= 2; ++band) {
        const int lower = block_lower(m, band);
        const int upper = std::min(block_upper(m, band), ambient_maximum);
        if (lower <= upper) {
            count += upper - lower + 1;
        }
    }
    return count;
}

int available_m1_pool_count(int ambient_maximum) {
    int count = 0;
    for (int band = 0; band <= 3; ++band) {
        const int lower = block_lower(1, band);
        const int upper = std::min(block_upper(1, band), ambient_maximum);
        if (lower <= upper) {
            count += upper - lower + 1;
        }
    }
    return count;
}

u64 uniform_family_count(int ambient_maximum) {
    u64 count = 0;
    for (int m = 2; m <= ambient_maximum; ++m) {
        const int available = available_pool_count(m, ambient_maximum);
        if (available >= 11) {
            const u64 addition = choose_u64(static_cast<unsigned>(available - 1),
                                            10);
            require(count <= std::numeric_limits<u64>::max() - addition,
                    "uniform family count overflow");
            count += addition;
        }
    }
    return count;
}

u64 exceptional_m1_family_count(int ambient_maximum) {
    const int available = available_m1_pool_count(ambient_maximum);
    return available >= 11
         ? choose_u64(static_cast<unsigned>(available - 1), 10) : 0;
}

struct FamilyCountRow {
    int ambient_maximum;
    u64 uniform_count;
    u64 exceptional_count;
    u64 combined_count;
};

std::vector<FamilyCountRow> audit_finite_family_counts() {
    const std::vector<std::pair<int, u64>> expected{
        {20, UINT64_C(75582)},
        {40, UINT64_C(812850987)},
        {80, UINT64_C(3595550244611)},
        {120, UINT64_C(397529462747261)},
        {160, UINT64_C(10616582432233990)},
        {200, UINT64_C(132777517674540845)}
    };
    std::vector<FamilyCountRow> rows;
    for (const auto& control : expected) {
        const u64 uniform = uniform_family_count(control.first);
        require(uniform == control.second,
                "finite uniform family count drift at N="
                + std::to_string(control.first));
        const u64 exceptional = exceptional_m1_family_count(control.first);
        require(uniform <= std::numeric_limits<u64>::max() - exceptional,
                "combined family count overflow");
        rows.push_back(FamilyCountRow{control.first, uniform, exceptional,
                                      uniform + exceptional});
    }
    require(rows[0].exceptional_count == choose_u64(15, 10),
            "N=20 exceptional count drift");
    require(rows[1].exceptional_count == choose_u64(21, 10),
            "N=40 exceptional count drift");
    for (std::size_t index = 2; index < rows.size(); ++index) {
        require(rows[index].exceptional_count == UINT64_C(1144066),
                "stabilized exceptional count drift");
    }
    return rows;
}

BigFraction density_term(const Rational& left_value,
                         const Rational& right_value,
                         const Rational& absolute_slope) {
    require(left_value.numerator >= 0 && right_value.numerator >= 0,
            "negative density endpoint");
    require(absolute_slope.numerator > 0,
            "nonpositive density slope");

    BigUInt right_term = power_small(
        static_cast<u64>(right_value.numerator), 11
    );
    BigUInt left_term = power_small(
        static_cast<u64>(left_value.numerator), 11
    );
    for (int count = 0; count < 11; ++count) {
        right_term.multiply_small(
            static_cast<std::uint32_t>(left_value.denominator)
        );
        left_term.multiply_small(
            static_cast<std::uint32_t>(right_value.denominator)
        );
    }
    if (right_term.compare(left_term) >= 0) {
        right_term.subtract(left_term);
    } else {
        left_term.subtract(right_term);
        right_term = left_term;
    }
    right_term.multiply_small(
        static_cast<std::uint32_t>(absolute_slope.denominator)
    );

    FactorMap denominator;
    add_factors(denominator,
                static_cast<u64>(left_value.denominator), 11);
    add_factors(denominator,
                static_cast<u64>(right_value.denominator), 11);
    add_factors(denominator,
                static_cast<u64>(absolute_slope.numerator), 1);
    return BigFraction{right_term, denominator};
}

struct DensityAudit {
    BigFraction density;
    BigFraction old_density;
    BigFraction gain;
    std::string density_string;
    std::string old_density_string;
    std::string gain_string;
};

DensityAudit audit_asymptotic_density() {
    // If x=m/N, the available-label length is the following continuous
    // piecewise-linear function.  The density is 11 integral f(x)^10 dx.
    // On a linear segment this equals |f(R)^11-f(L)^11|/|slope|.
    const std::vector<BigFraction> terms{
        density_term(Rational(0), Rational(21, 41), Rational(63, 4)),
        density_term(Rational(21, 41), Rational(14, 29), Rational(15)),
        density_term(Rational(14, 29), Rational(56, 81), Rational(14)),
        density_term(Rational(56, 81), Rational(7, 12), Rational(25, 4)),
        density_term(Rational(7, 12), Rational(35, 39), Rational(35, 4)),
        density_term(Rational(35, 39), Rational(0), Rational(1))
    };
    BigFraction density = sum_positive_fractions(terms);
    BigFraction old_density{power_small(35, 10), FactorMap{}};
    add_factors(old_density.denominator_factors, 39, 10);
    reduce_fraction(old_density);
    BigFraction gain = subtract_fractions(density, old_density);

    const std::string density_string = density.numerator.to_string() + "/"
        + denominator_string(density.denominator_factors);
    const std::string old_density_string = old_density.numerator.to_string()
        + "/" + denominator_string(old_density.denominator_factors);
    const std::string gain_string = gain.numerator.to_string() + "/"
        + denominator_string(gain.denominator_factors);

    require(density_string
            == "848953086913769850118498851618778832628468542103282298683365532079/"
               "2481088067163593416217816176836483026480276818419826456353950662656",
            "asymptotic density rational drift");
    require(old_density_string
            == "2758547353515625/8140406085191601",
            "old density rational drift");
    require(gain_string
            == "59367779473724913286745294455473885142894795652409271/"
               "17997353908972062863364250852276005638984892348209823744",
            "density gain rational drift");
    return DensityAudit{density, old_density, gain, density_string,
                        old_density_string, gain_string};
}

}  // namespace

int main() {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    try {
        const PoolAudit pool = audit_complete_pools();
        const M1Audit m1 = audit_m1_boundary();
        const EndpointAudit endpoint = audit_endpoint_and_tooth();
        const PairAudit pairs = audit_direct_tail_pairs();
        const M7Audit m7 = audit_m7_family();
        const std::vector<FamilyCountRow> family_rows =
            audit_finite_family_counts();
        const DensityAudit density = audit_asymptotic_density();

        std::ostringstream uniform_counts;
        std::ostringstream exceptional_counts;
        std::ostringstream combined_counts;
        for (std::size_t index = 0; index < family_rows.size(); ++index) {
            if (index != 0) {
                uniform_counts << ',';
                exceptional_counts << ',';
                combined_counts << ',';
            }
            const FamilyCountRow& row = family_rows[index];
            uniform_counts << row.ambient_maximum << ':' << row.uniform_count;
            exceptional_counts << row.ambient_maximum << ':'
                               << row.exceptional_count;
            combined_counts << row.ambient_maximum << ':' << row.combined_count;
        }

        SemanticLedger ledger;
        ledger.add("status", "ACCEPT_CANDIDATE_AS_STATED");
        ledger.add("uniform_m_range", "2..1000");
        ledger.add("three_blocks", "[m,floor(13s/16)];[15m,floor(27s/16)];[29m,floor(41s/16)]");
        ledger.add("interval", "[1/(14m),8/(7(12m+1))]");
        ledger.add("width", "(4m-1)/(14m(12m+1))");
        ledger.add("speed_gates", pool.speeds_tested);
        ledger.add("pool_speed_gates", pool.pool_speeds);
        ledger.add("outside_speed_gates", pool.outside_speeds);
        ledger.add("fourth_band_min_defect",
                   static_cast<u64>(pool.minimum_fourth_band_defect));
        ledger.add("m1_blocks", "[1,10];[15,21];[29,33];[43,44]");
        ledger.add("m1_pool_size", static_cast<u64>(m1.pool_size));
        ledger.add("m1_fifth_band_defect",
                   static_cast<u64>(m1.fifth_band_defect));
        ledger.add("m1_eleven_body_count", m1.family_count);
        ledger.add("small_tail_endpoint_gates", endpoint.small_tail_endpoint_gates);
        ledger.add("large_tail_room_gates", endpoint.large_tail_room_gates);
        ledger.add("tooth_gates", endpoint.tooth_gates);
        ledger.add("direct_tail_pair_box", "m=1..100;odd tails<=101");
        ledger.add("direct_tail_pairs", pairs.pairs);
        ledger.add("direct_phase_gates", pairs.phase_gates);
        ledger.add("m7_blocks", "[7,69];[105,143];[203,217]");
        ledger.add("m7_pool_size", static_cast<u64>(m7.pool_size));
        ledger.add("m7_anchors", "7,120,126,143");
        ledger.add("m7_free_size", static_cast<u64>(m7.free_size));
        ledger.add("m7_family_count", m7.family_count);
        ledger.add("m7_all_eleven_bodies_containing_7", m7.all_m7_eleven_bodies);
        ledger.add("m7_family_fraction_in_pool",
                   std::to_string(m7.fraction_in_pool.first) + "/"
                   + std::to_string(m7.fraction_in_pool.second));
        ledger.add("m7_family_fraction_in_ambient_1_to_217",
                   std::to_string(m7.fraction_in_ambient.first) + "/"
                   + std::to_string(m7.fraction_in_ambient.second));
        ledger.add("finite_uniform_family_counts", uniform_counts.str());
        ledger.add("finite_m1_exception_counts", exceptional_counts.str());
        ledger.add("finite_combined_family_counts", combined_counts.str());
        ledger.add("asymptotic_density", density.density_string);
        ledger.add("old_density", density.old_density_string);
        ledger.add("density_gain", density.gain_string);
        ledger.add("supplementary_pool_label_density", "21/41");
        ledger.add("scope", "LRC14_OPEN;ENTRY_OPEN;M1_FINITE_SLICE_DOES_NOT_CHANGE_DENSITY");

        std::cout << "THM4158_INDEPENDENT_CPP_AUDIT\n";
        for (const std::string& row : ledger.rows()) {
            std::cout << row << '\n';
        }
        std::ostringstream hash_stream;
        hash_stream << std::hex << std::setw(16) << std::setfill('0')
                    << ledger.hash();
        std::cout << "semantic_fnv64=" << hash_stream.str() << '\n';
        std::cout << "result=PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "result=FAIL\nreason=" << error.what() << '\n';
        return 1;
    }
}
