#!/usr/bin/env python3
"""Exact, dependency-free checks for the finite fourth-power budget note."""
from collections import Counter
from itertools import product
from math import factorial
import json


class CheckFailure(RuntimeError):
    pass


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise CheckFailure(message)


def nonnegative_power_values(p, limit):
    values = []
    base = 0
    while base ** p <= limit:
        values.append(base ** p)
        base += 1
    return values


def ordered_dp(values, slots, limit):
    """Count exact ordered choices with repetition by convolution."""
    counts = [0] * (limit + 1)
    counts[0] = 1
    for _ in range(slots):
        following = [0] * (limit + 1)
        for total, multiplicity in enumerate(counts):
            if multiplicity:
                for value in values:
                    if total + value <= limit:
                        following[total + value] += multiplicity
        counts = following
    return counts


def minimum_terms(values, limit):
    """Unbounded-coin DP; exclude zero to avoid cyclic recurrences."""
    counts = [0] + [limit + 1] * limit
    for target in range(1, limit + 1):
        counts[target] = min(
            (counts[target - value] + 1 for value in values if 0 < value <= target),
            default=limit + 1,
        )
    return counts


def canonical_minimum(n, B):
    q, r = divmod(n, B)
    return q + r


def count_pairs(n, B, k):
    q, r = divmod(n, B)
    if q + r > k:
        return []
    T = min(q, (k - q - r) // (B - 1))
    return [
        (r + B * t, q - t, k - q - r - (B - 1) * t)
        for t in range(T + 1)
    ]


def ordered_formula(n, B, k, signed=False):
    result = 0
    for a, b, c in count_pairs(n, B, k):
        coefficient = factorial(k) // (factorial(a) * factorial(b) * factorial(c))
        result += coefficient * (2 ** (a + b) if signed else 1)
    return result


def cutoff_count(B, k, N):
    return sum(
        1 + min(B - 1, k - q, N - B * q)
        for q in range(min(k, N // B) + 1)
    )


def full_count(B, k):
    if k <= B - 1:
        return (k + 1) * (k + 2) // 2
    return B * (k + 1) - B * (B - 1) // 2


def runs(values):
    result = []
    for value in sorted(values):
        if result and result[-1][1] + 1 == value:
            result[-1][1] = value
        else:
            result.append([value, value])
    return result


def exact_fourth_power_case():
    expected = [0, 1, 2, 3, 4, 16, 17, 18, 19, 32, 33, 34, 48, 49]
    tuples = list(product(range(3), repeat=4))
    direct = Counter(sum(base ** 4 for base in t) for t in tuples)
    through50 = {n: multiplicity for n, multiplicity in sorted(direct.items()) if n <= 50}
    require(list(through50) == expected, "direct fourth-power attainable values")
    require(len(through50) == 14, "14 nonnegative attained values")
    require(len([n for n in through50 if n > 0]) == 13, "13 positive attained values")
    missing = set(range(51)) - set(through50)
    require(len(missing) == 37, "37 missing nonnegative-range values")
    require(set(range(1, 51)) - set(through50) == missing, "positive-range complement unchanged")
    require(runs(missing) == [[5, 15], [20, 31], [35, 47], [50, 50]], "exact gap runs")
    require([hi - lo + 1 for lo, hi in runs(missing)] == [11, 12, 13, 1], "gap lengths")
    dp = ordered_dp([0, 1, 16], 4, 80)
    require(all(dp[n] == direct[n] for n in range(81)), "direct enumeration versus convolution")
    require(all(dp[n] == ordered_formula(n, 16, 4) for n in range(81)), "multinomial coefficients")
    expected_counts = [1, 4, 6, 4, 1, 4, 12, 12, 4, 6, 12, 6, 4, 4]
    require(list(through50.values()) == expected_counts, "exact table coefficients")
    require(sum(through50.values()) == 80, "80 ordered nonnegative tuples")
    require(sum(c for n, c in through50.items() if n > 0) == 79, "79 tuples for positive targets")
    require(len(direct) == 15 and max(direct) == 64, "full fourth-power palette")
    require(len(set(range(65)) - set(direct)) == 50, "full missing count")
    require(runs(set(range(65)) - set(direct))[-1] == [50, 63], "singleton extends to full gap")
    minima = minimum_terms(nonnegative_power_values(4, 80), 80)
    require(all(minima[n] == canonical_minimum(n, 16) for n in range(81)), "minimum terms through 80")
    require(minima[37] == 7 and minima[50] == 5, "37 and 50 exact minimum terms")
    fourth_residues = {base ** 4 % 16 for base in range(16)}
    require(fourth_residues == {0, 1}, "all fourth-power residues modulo 16")
    residue_sums = {sum(t) % 16 for t in product(fourth_residues, repeat=4)}
    require(residue_sums == {0, 1, 2, 3, 4}, "four-summand residues")
    require(37 % 16 not in residue_sums and 50 % 16 in residue_sums, "local versus budget obstruction")
    positive = Counter(sum(base ** 4 for base in t) for t in product([1, 2], repeat=4))
    require(sorted(n for n in positive if n <= 50) == [4, 19, 34, 49], "strict-positive convention")
    require(sum(c for n, c in positive.items() if n <= 50) == 15, "strict-positive tuple count")
    signed = Counter(sum(base ** 4 for base in t) for t in product(range(-2, 3), repeat=4))
    require(sorted(n for n in signed if n <= 50) == expected, "negative bases preserve attained set")
    require(sum(c for n, c in signed.items() if n <= 50) == 609, "signed ordered tuple count")
    require(all(signed[n] == ordered_formula(n, 16, 4, signed=True) for n in range(81)), "signed multinomials")
    minima81 = minimum_terms(nonnegative_power_values(4, 81), 81)
    require(minima81[81] == 1 and canonical_minimum(81, 16) == 6, "81 breaks unguarded minimum")
    actual81 = ordered_dp(nonnegative_power_values(4, 81), 4, 81)
    require(actual81[81] == 4 and canonical_minimum(81, 16) > 4, "81 breaks unguarded existence")
    require(cutoff_count(16, 4, 49) - 1 == 13 and 49 - 13 == 36, "cutoff 49 sensitivity")
    require(cutoff_count(16, 4, 51) - 1 == 13 and 51 - 13 == 38, "cutoff 51 sensitivity")
    return {
        "attained_nonnegative_through_50": expected,
        "attained_positive_count": 13,
        "missing_count": 37,
        "missing_intervals": runs(missing),
        "ordered_nonnegative_tuple_count": 80,
        "ordered_signed_tuple_count": 609,
        "minimum_terms_37": minima[37],
        "minimum_terms_50": minima[50],
        "first_new_fourth_power": {"n": 81, "actual_minimum": minima81[81], "restricted_minimum": 6},
    }


def general_power_checks():
    cases = []
    for p in [2, 4, 6]:
        B = 2 ** p
        limit = 3 ** p - 1
        values = nonnegative_power_values(p, limit)
        require(values == [0, 1, B], "exact palette below third power")
        minima = minimum_terms(values, limit)
        require(all(minima[n] == canonical_minimum(n, B) for n in range(limit + 1)), "general minimum theorem")
        budgets = sorted({0, 1, 2, 4, B - 1, B, B + 2})
        for k in budgets:
            dp = ordered_dp(values, k, limit)
            for n in range(limit + 1):
                pairs = count_pairs(n, B, k)
                brute_pairs = [(a, b, k - a - b) for b in range(k + 1)
                               for a in [n - B * b] if 0 <= a and a + b <= k]
                require(sorted(pairs) == sorted(brute_pairs), "collision parameterization")
                require(bool(dp[n]) == (canonical_minimum(n, B) <= k), "budget existence iff")
                require(dp[n] == ordered_formula(n, B, k), "general multinomial formula")
                positive_exists = (n >= k and (n - k) % (B - 1) == 0
                                   and (n - k) // (B - 1) <= k)
                require(positive_exists == any(c == 0 for _, _, c in pairs), "strict-positive iff")
            attained = {n for n in range(limit + 1) if dp[n]}
            for N in sorted({0, 1, min(B, limit), min(B * k, limit), limit}):
                require(cutoff_count(B, k, N) == len([n for n in attained if n <= N]), "general cutoff count")
        require(len(count_pairs(B, B, B - 1)) == 1, "no collision below B slots")
        require(len(count_pairs(B, B, B)) == 2, "first collision at B slots")
        first_new = 3 ** p
        require(canonical_minimum(first_new, B) > 1, "new denomination strictly improves minimum")
        cases.append({"p": p, "B": B, "maximum_n": limit, "budgets": budgets})
    return cases


def palette_and_gap_checks():
    cases = 0
    for B in range(2, 25):
        for k in range(0, B + 4):
            attained = {a + B * b for b in range(k + 1) for a in range(k - b + 1)}
            require(len(attained) == full_count(B, k), "full palette piecewise formula")
            require(len(attained) == sum(min(B, j + 1) for j in range(k + 1)), "canonical full count")
            require(cutoff_count(B, k, B * k) == len(attained), "full versus cutoff formula")
            if k < B:
                pair_count = (k + 1) * (k + 2) // 2
                require(len(attained) == pair_count, "triangle injectivity")
                gaps = [B - k + b - 1 for b in range(k)]
                observed_gaps = [B * (b + 1) - (B * b + k - b) - 1 for b in range(k)]
                require(gaps == observed_gaps and all(g >= 0 for g in gaps), "gap formula")
                missing_count = B * k + 1 - len(attained)
                require(missing_count == sum(gaps), "gaps exhaust complement")
                require(missing_count == k * (B - k - 1) + k * (k - 1) // 2, "summed gap formula")
            else:
                require(len(attained) < (k + 1) * (k + 2) // 2, "triangle overcounts after collision")
                require(B * k + 1 - len(attained) == (B - 1) * (B - 2) // 2, "post-collision complement plateau")
            cases += 1
    return cases


def hostile_controls():
    """Reject tempting false statements explicitly, without assert."""
    require(not (50 % 16 <= 4 and canonical_minimum(50, 16) <= 4), "residue sufficiency is false")
    require(cutoff_count(16, 4, 50) != 13, "13 is not the count including zero")
    require(sum(ordered_dp([0, 1, 16], 4, 50)) != 14, "representations are not distinct values")
    require(full_count(4, 4) != (4 + 1) * (4 + 2) // 2, "unguarded triangular count is false")
    require(3 < 4 and 4 * 3 >= 3 ** 2, "k < B does not enforce third-power range")
    square_dp = ordered_dp(nonnegative_power_values(2, 10), 3, 10)
    require(square_dp[10] > 0 and canonical_minimum(10, 4) > 3, "square cutoff counterexample")
    require(3 ** 6 % (2 ** 6) == 25, "global zero-one residue claim does not generalize to sixth powers")
    wrong_gaps = [16 - 4 + b for b in range(4)]
    require(wrong_gaps != [11, 12, 13, 14], "gap off-by-one rejected")


def main():
    fourth = exact_fourth_power_case()
    general = general_power_checks()
    palette_cases = palette_and_gap_checks()
    hostile_controls()
    print(json.dumps({
        "status": "PASS",
        "checks": CHECKS,
        "arithmetic": "exact integers only",
        "assert_statements_used": False,
        "fourth_power_case": fourth,
        "general_power_cases": general,
        "full_palette_parameter_cases": palette_cases,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
