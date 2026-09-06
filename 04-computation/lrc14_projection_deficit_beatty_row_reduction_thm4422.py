#!/usr/bin/env python3
"""Exact audit for THM-4422 projection-deficit and Beatty-row reductions.

Everything is integer or Fraction arithmetic. The scoped reductions are
proved analytically; this verifier does not assert the remaining projection
inequality or LRC(14).
"""

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import gcd
import sys

TARGET = Q(6, 77)
RUNNER_UP = Q(12, 161)
AP_MAX = Q(24, 343)
CHECKS = 0


def need(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(payload)


def strict_floor(num, den):
    return (num - 1) // den


def ceil_div(num, den):
    return -((-num) // den)


def egcd(a, b):
    if b == 0:
        return (a, 1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)


def bounds(w):
    return tuple(strict_floor(3 * (sum(w) - w[i]), 14) for i in range(3))


def carriers_box(w):
    """Definition-level carrier enumeration in ker(w)."""
    a, b, c = w
    bx, by, bz = bounds(w)
    result = set()
    for x in range(-bx, bx + 1):
        for y in range(-by, by + 1):
            nz = -(a * x + b * y)
            if nz % c:
                continue
            z = nz // c
            C = (x, y, z)
            if abs(z) <= bz and all(v % 3 for v in C):
                result.add(C)
    return result


def projection_sums(w, carriers):
    """Return (S_1,S_2,S_3) in coordinate, not pair-name, order."""
    q = Q(3, 7 * w[2])
    out = []
    for i in range(3):
        j, k = tuple(u for u in range(3) if u != i)
        total = Q(0)
        for C in carriers:
            p = 3 * (w[j] + w[k]) - 14 * abs(C[i])
            need(p > 0, ("live-roof", w, C, i, p))
            total += min(q, Q(p, 14 * w[j] * w[k]))
        out.append(total)
    return tuple(out)


def normalized_deficits(w, carriers):
    """Return E_i with S_i=q(|Lambda|-E_i)."""
    a, b, c = w
    q = Q(3, 7 * c)
    result = []
    for i in range(3):
        j, k = tuple(u for u in range(3) if u != i)
        total = Q(0)
        for C in carriers:
            p = 3 * (w[j] + w[k]) - 14 * abs(C[i])
            total += max(Q(0), 1 - Q(c * p, 6 * w[j] * w[k]))
        result.append(total)
    sums = projection_sums(w, carriers)
    need(all(sums[i] == q * (len(carriers) - result[i]) for i in range(3)),
         ("deficit-identity", w, sums, result))
    need(
        (min(sums) <= TARGET)
        == (max(result) >= Q(len(carriers), 1) - Q(2 * c, 11)),
        ("deficit-target-equivalence", w, sums, result),
    )
    return tuple(result)


def first_coordinate_rows(w):
    """Reconstruct carriers by the exact one-dimensional C_1 rows."""
    a, b, c = w
    d = gcd(b, c)
    need(gcd(a, d) == 1, ("primitive-row", w, d))
    b0, c0 = b // d, c // d
    g, u, v = egcd(b0, c0)
    need(g == 1 and b0 * u + c0 * v == 1, ("bezout", w, b0, c0, u, v))
    bx, by, bz = bounds(w)
    result = set()
    row_counts = {}
    uniform_bound = ceil_div(d * (a + c), 7 * c)
    for n in range(-(bx // d), bx // d + 1):
        if n % 3 == 0:
            continue
        y0, z0 = -a * n * u, -a * n * v
        lo = max(ceil_div(-by - y0, c0), ceil_div(z0 - bz, b0))
        hi = min((by - y0) // c0, (z0 + bz) // b0)
        residues = [
            r for r in range(3)
            if (y0 + c0 * r) % 3 and (z0 - b0 * r) % 3
        ]
        need(len(residues) == 1, ("unique-live-residue", w, n, residues))
        tau = residues[0]
        first = lo + ((tau - lo) % 3)
        ts = range(first, hi + 1, 3) if first <= hi else ()
        count = 0
        for t in ts:
            C = (d * n, y0 + c0 * t, z0 - b0 * t)
            need(sum(w[i] * C[i] for i in range(3)) == 0, ("row-kernel", w, C))
            result.add(C)
            count += 1
        row_counts[n] = count
        need(count <= uniform_bound, ("row-multiplicity-bound", w, n, count, uniform_bound))
    return result, row_counts, uniform_bound


def predicted_slope_two(m, epsilon):
    if epsilon == 1:
        limit_num = 3 * (m + 1)
        base = (-1, -2, 1)
    else:
        limit_num = 3 * m
        base = (1, -2, 1)
    ts = tuple(t for t in range(1, strict_floor(limit_num, 14) + 1) if t % 3)
    positive = {tuple(t * x for x in base) for t in ts}
    return positive | {tuple(-x for x in C) for C in positive}, len(ts)


def ap_predicted_carriers(w):
    a, b, c = w
    need(c == 2 * b - a, ("not-ap", w))
    K = strict_floor(3 * b, 14)
    positive = {
        (k, -2 * k, k)
        for k in range(1, K + 1)
        if k % 3
    }
    return positive | {tuple(-x for x in C) for C in positive}


def ap_middle_projection(w):
    """Closed formula for S_2 on c=2b-a."""
    a, b, c = w
    K = strict_floor(3 * b, 14)
    T = sum(
        (min(Q(1), Q(3 * b - 14 * k, 3 * a))
         for k in range(1, K + 1) if k % 3),
        Q(0),
    )
    return Q(6, 7 * c) * T


def arithmetic_progression_family():
    """Finite base and wide exact control for the all-height AP proof."""
    base = []
    wide_rows = 0
    wide_best = (Q(0), None)
    for c in range(1, 500, 2):
        if c % 3 == 0:
            continue
        for b in range(1, c, 2):
            a = 2 * b - c
            if not (0 < a < b and a % 3 and b % 3 and gcd(a, b) == 1):
                continue
            w = (a, b, c)
            formula = ap_middle_projection(w)
            if c < 67:
                direct = carriers_box(w)
                predicted = ap_predicted_carriers(w)
                need(direct == predicted, ("ap-carrier-collapse", w, direct ^ predicted))
                actual = projection_sums(w, direct)[1]
                need(actual == formula, ("ap-middle-formula", w, actual, formula))
                base.append((actual, w))
            if formula > wide_best[0]:
                wide_best = (formula, w)
            wide_rows += 1
    need(len(base) == 49, ("ap-finite-base-count", len(base)))
    need(max(base) == (AP_MAX, (1, 25, 49)), ("ap-finite-base-max", max(base)))
    need(wide_best == (AP_MAX, (1, 25, 49)), ("ap-wide-max", wide_best))
    need(Q(3, 49) + Q(4, 7 * 67) < AP_MAX < TARGET,
         ("ap-analytic-tail", Q(3, 49) + Q(4, 7 * 67), AP_MAX, TARGET))
    return len(base), wide_rows, wide_best


def primitive_direction(C):
    g = gcd(gcd(abs(C[0]), abs(C[1])), abs(C[2]))
    D = tuple(x // g for x in C)
    if next(x for x in D if x) < 0:
        D = tuple(-x for x in D)
    return D


def short_relation_data(w):
    """Return (doubled coordinate, ray, cap A, roof B) for a norm-four relation."""
    a, b, c = w
    if c == 2 * a + b:
        return 0, (2, 1, -1), Q(3 * a, 14), Q(3 * (a + b), 14)
    if c == 2 * b - a:
        return 1, (1, -2, 1), Q(3 * (b - a), 14), Q(3 * b, 14)
    if c == a + 2 * b:
        return 1, (1, 2, -1), Q(3 * b, 14), Q(3 * (a + b), 14)
    return None


def short_relation_formula(w):
    r, ray, A, B = short_relation_data(w)
    K = strict_floor(B.numerator, B.denominator)
    W = B - A
    T = sum(
        (min(Q(1), (B - k) / W) for k in range(1, K + 1) if k % 3),
        Q(0),
    )
    return Q(6, 7 * w[2]) * T, r, ray, K


def norm_four_relation_families():
    """Finite base and wide controls for all three signed (1,1,2) relations."""
    base = []
    wide = []
    type_counts = Counter()
    for c in range(1, 500, 2):
        if c % 3 == 0:
            continue
        for a in range(1, c, 2):
            if a % 3 == 0:
                continue
            for b in range(a + 2, c, 2):
                if b % 3 == 0 or gcd(gcd(a, b), c) != 1:
                    continue
                w = (a, b, c)
                data = short_relation_data(w)
                if data is None:
                    continue
                formula, r, ray, K = short_relation_formula(w)
                type_counts[ray] += 1
                wide.append((formula, w, r))
                if c < 67:
                    direct = carriers_box(w)
                    predicted = {
                        tuple(k * x for x in ray)
                        for k in range(-K, K + 1)
                        if k and k % 3
                    }
                    need(direct == predicted,
                         ("norm-four-carrier-collapse", w, ray, direct ^ predicted))
                    actual = projection_sums(w, direct)[r]
                    need(actual == formula,
                         ("norm-four-selected-formula", w, r, actual, formula))
                    base.append((actual, w, r))
    need(len(base) == 105, ("norm-four-base-count", len(base)))
    need(max(base) == (TARGET, (1, 5, 11), 1), ("norm-four-base-max", max(base)))
    ordered = sorted(wide, reverse=True)
    need(ordered[:3] == [
        (TARGET, (1, 5, 11), 1),
        (RUNNER_UP, (1, 11, 23), 1),
        (AP_MAX, (1, 25, 49), 1),
    ], ("norm-four-wide-leaders", ordered[:3]))
    need(Q(3, 49) + Q(4, 7 * 67) < AP_MAX,
         ("norm-four-tail-below-third", Q(3, 49) + Q(4, 7 * 67), AP_MAX))
    return len(base), len(wide), type_counts, ordered[:3]


def averaging_stopper():
    named = {
        (1, 5, 7): (Q(8, 245), Q(6, 49), Q(4, 35)),
        (5, 11, 17): (Q(122, 1309), Q(8, 119), Q(12, 119)),
    }
    for w, expected in named.items():
        actual = projection_sums(w, carriers_box(w))
        need(actual == expected, ("named-projections", w, actual, expected))

    # First row forces alpha_1 >= 49/110; second forces alpha_2 >= 10/17.
    force_1, force_2 = Q(49, 110), Q(10, 17)
    need(force_1 + force_2 == Q(1933, 1870) > 1,
         ("incompatible-fixed-weights", force_1, force_2))
    return named, force_1, force_2


def main():
    sys.stdout.reconfigure(newline="\n")
    # General row formula, hostile-audited on every eligible triple to height 79.
    values = tuple(v for v in range(1, 80, 2) if v % 3)
    triples = 0
    rows = 0
    max_multiplicity = 0
    sparse_automatic = 0
    dense_deficit = 0
    dense_directions = Counter()
    dense_winners = Counter()
    dense_gcd_profiles = Counter()
    dense_exceptions = []
    for w in combinations(values, 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        direct = carriers_box(w)
        rebuilt, counts, _ = first_coordinate_rows(w)
        need(direct == rebuilt, ("row-reconstruction", w, direct ^ rebuilt))
        deficits = normalized_deficits(w, direct)
        if Q(len(direct), 1) <= Q(2 * w[2], 11):
            sparse_automatic += 1
        else:
            need(max(deficits) >= Q(len(direct), 1) - Q(2 * w[2], 11),
                 ("bounded-dense-deficit", w, len(direct), deficits))
            dense_deficit += 1
            directions = {primitive_direction(C) for C in direct}
            winner = max(range(3), key=lambda i: deficits[i])
            dense_gcd_profiles[(gcd(w[1], w[2]), gcd(w[0], w[2]), gcd(w[0], w[1]))] += 1
            if len(directions) == 1:
                direction = next(iter(directions))
                dense_directions[direction] += 1
                dense_winners[(direction, winner)] += 1
            else:
                dense_exceptions.append((w, tuple(sorted(directions)), winner))
        triples += 1
        rows += len(counts)
        max_multiplicity = max(max_multiplicity, max(counts.values(), default=0))

    family_rows = 0
    family_equality = []
    best_after_equality = (Q(0), None)
    sample = {}
    for epsilon, start in ((1, 5), (-1, 7)):
        for m in range(start, 500, 6):
            w = (1, m, 2 * m + epsilon)
            direct = carriers_box(w)
            predicted, n = predicted_slope_two(m, epsilon)
            need(direct == predicted, ("slope-two-carriers", w, direct ^ predicted))
            sums = projection_sums(w, direct)
            cap_bound = Q(6 * n, 7 * w[2])
            need(all(value <= cap_bound for value in sums),
                 ("slope-two-cap", w, sums, cap_bound))
            need(min(sums) <= TARGET, ("slope-two-target", w, sums))
            if min(sums) == TARGET:
                family_equality.append(w)
            elif min(sums) > best_after_equality[0]:
                best_after_equality = (min(sums), w)
            if w in ((1, 5, 11), (1, 11, 23), (1, 25, 49)):
                sample[w] = (n, sums, cap_bound)
            family_rows += 1
    need(family_equality == [(1, 5, 11)], ("family-equality", family_equality))
    need(best_after_equality == (RUNNER_UP, (1, 11, 23)),
         ("family-runner-up", best_after_equality))

    named, force_1, force_2 = averaging_stopper()
    ap_base, ap_rows, ap_best = arithmetic_progression_family()
    norm4_base, norm4_rows, norm4_types, norm4_leaders = norm_four_relation_families()
    need(dense_gcd_profiles == {(1, 1, 1): 114},
         ("dense-gcd-classification", dense_gcd_profiles))
    need(dense_directions == {(2, 1, -1): 55, (1, -2, 1): 52, (1, 2, -1): 6},
         ("dense-direction-classification", dense_directions))
    need(dense_winners == {
        ((2, 1, -1), 0): 55,
        ((1, -2, 1), 1): 52,
        ((1, 2, -1), 1): 6,
    }, ("dense-winner-classification", dense_winners))
    need(len(dense_exceptions) == 1 and dense_exceptions[0][0] == (19, 23, 29),
         ("dense-exception", dense_exceptions))
    exception_w = (19, 23, 29)
    exception_carriers = carriers_box(exception_w)
    expected_exception_carriers = {
        (1, 8, -7), (10, -7, -1), (11, 1, -8),
        (-1, -8, 7), (-10, 7, 1), (-11, -1, 8),
    }
    need(exception_carriers == expected_exception_carriers,
         ("dense-exception-carriers", exception_carriers))
    exception_sums = projection_sums(exception_w, exception_carriers)
    need(exception_sums == (Q(156, 4669), Q(192, 3857), Q(3840, 88711)),
         ("dense-exception-projections", exception_sums))
    print("LRC14 PROJECTION DEFICIT + BEATTY-ROW REDUCTION THM-4422")
    print("status=PASS PROVED_REDUCTIONS+FINITE_EXACT_CONTROLS; universal_inequality=OPEN; LRC14=OPEN")
    print("row_formula_height=79 triples=%d first_coordinate_rows=%d max_row_multiplicity=%d"
          % (triples, rows, max_multiplicity))
    print("carrier_count_split sparse_automatic=%d dense_needing_deficit=%d"
          % (sparse_automatic, dense_deficit))
    print("dense_H79 gcd_profiles=%s directions=%s winners=%s exception=%s"
          % (dict(dense_gcd_profiles), dict(dense_directions), dict(dense_winners), dense_exceptions))
    print("dense_exception_carriers=%s projections=%s"
          % (tuple(sorted(exception_carriers)), tuple(str(x) for x in exception_sums)))
    print("slope_two_family_rows_through_c499=%d equality=%s runner_up=%s"
          % (family_rows, family_equality, best_after_equality))
    print("arithmetic_progression_family c=2b-a finite_base=%d formula_rows_c499=%d best=%s"
          % (ap_base, ap_rows, ap_best))
    print("norm_four_relation_families finite_base_c67=%d formula_rows_c499=%d types=%s leaders=%s"
          % (norm4_base, norm4_rows, dict(norm4_types), norm4_leaders))
    for w in sorted(sample):
        print("family_control", w, "positive_t_count", sample[w][0],
              "S", tuple(str(v) for v in sample[w][1]), "cap", sample[w][2])
    print("fixed_weight_stopper_rows=")
    for w, values in named.items():
        print(" ", w, tuple(str(v) for v in values))
    print("forced_alpha1=%s forced_alpha2=%s sum=%s"
          % (force_1, force_2, force_1 + force_2))
    print("checks=%d verdict=PASS" % CHECKS)


if __name__ == "__main__":
    main()
