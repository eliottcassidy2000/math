#!/usr/bin/env python3
"""Exact hostile audit for the exponent-swap / LRC pair bridge.

The executable deliberately uses only integer arithmetic and Fraction.  Its
primary path reconstructs owner-cell lengths from the sorted equality clocks;
the closed formulas are evaluated on a separate path.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def distance_to_integer(x):
    x %= 1
    return min(x, 1 - x)


def pair_candidate_points(m, n):
    points = {Fraction(0)}
    for speed in (m, n):
        for j in range(2 * speed):
            points.add(Fraction(j, 2 * speed))
    for clock in (n - m, n + m):
        for j in range(clock):
            points.add(Fraction(j, clock))
    return points


def exact_pair_maximum(m, n):
    values = []
    for t in pair_candidate_points(m, n):
        value = min(distance_to_integer(m * t), distance_to_integer(n * t))
        values.append((value, t))
    maximum = max(value for value, _ in values)
    maximizers = tuple(sorted(t for value, t in values if value == maximum))
    return maximum, maximizers


def exact_owner_measure(m, n):
    """Measure of {t: ||m t|| < ||n t||} from its exact cell decomposition."""
    equality = {Fraction(j, n - m) for j in range(n - m)}
    equality.update(Fraction(j, n + m) for j in range(n + m))
    cuts = sorted(equality)
    total = Fraction(0)
    for index, left in enumerate(cuts):
        right = cuts[index + 1] if index + 1 < len(cuts) else Fraction(1)
        midpoint = (left + right) / 2
        if distance_to_integer(m * midpoint) < distance_to_integer(n * midpoint):
            total += right - left
    return total, len(cuts)


def closed_pair_data(m, n):
    require(0 < m < n and gcd(m, n) == 1, "primitive ordered pair required")
    both_odd = (m % 2 == 1 and n % 2 == 1)
    if both_odd:
        maximum = Fraction(1, 2)
        owner_measure = Fraction(1, 2)
        tie_count = 2 * n - 2
        maximizers = (Fraction(1, 2),)
    else:
        q = m + n
        maximum = Fraction(q - 1, 2 * q)
        owner_measure = Fraction(1, 2) + Fraction(1, 2 * (n * n - m * m))
        tie_count = 2 * n - 1
        residue = ((q - 1) // 2) * pow(m, -1, q) % q
        maximizers = tuple(sorted((Fraction(residue, q), Fraction((-residue) % q, q))))
    return maximum, owner_measure, tie_count, maximizers


def closed_scaled_pair_data(u, v):
    require(0 < u < v, "ordered distinct positive speeds required")
    common = gcd(u, v)
    m, n = u // common, v // common
    maximum, measure, _, core_maximizers = closed_pair_data(m, n)
    same_level = valuation_two(u) == valuation_two(v)
    tie_count = 2 * v - (2 * common if same_level else common)
    maximizers = set()
    for core_time in core_maximizers:
        for sheet in range(common):
            maximizers.add((core_time + sheet) / common)
    return maximum, measure, tie_count, tuple(sorted(maximizers))


def exact_multi_maximum(speeds):
    speeds = tuple(sorted(set(speeds)))
    points = {Fraction(0)}
    for speed in speeds:
        for j in range(2 * speed):
            points.add(Fraction(j, 2 * speed))
    for i, left in enumerate(speeds):
        for right in speeds[i + 1 :]:
            for clock in (right - left, right + left):
                for j in range(clock):
                    points.add(Fraction(j, clock))
    values = []
    for t in points:
        value = min(distance_to_integer(speed * t) for speed in speeds)
        values.append((value, t))
    maximum = max(value for value, _ in values)
    maximizers = tuple(sorted(t for value, t in values if value == maximum))
    return maximum, maximizers


def factor_small(n):
    result = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def exponent_core(a, b):
    fa = factor_small(a)
    fb = factor_small(b)
    primes = set(fa) | set(fb)
    signed = {p: b * fa.get(p, 0) - a * fb.get(p, 0) for p in primes}
    left = 1
    right = 1
    for p, exponent in signed.items():
        if exponent > 0:
            left *= p**exponent
        elif exponent < 0:
            right *= p ** (-exponent)
    return left, right, signed


def valuation_two(n):
    value = 0
    while n % 2 == 0:
        value += 1
        n //= 2
    return value


def labelled_pair_summary(left, right):
    common = gcd(left, right)
    x, y = left // common, right // common
    m, n = sorted((x, y))
    maximum, smaller_measure, _, _ = closed_pair_data(m, n)
    if x < y:
        directed_measure = smaller_measure
    else:
        directed_measure = 1 - smaller_measure
    return maximum, directed_measure


def decode_labelled_ratio(maximum, directed_measure):
    """Recover left/right from a nonzero-imbalance pair summary."""
    defect = Fraction(1, 2) - maximum
    imbalance = directed_measure - Fraction(1, 2)
    require(defect > 0 and imbalance != 0, "mixed-valuation summary required")
    q = Fraction(1, 2 * defect)
    r = defect / abs(imbalance)
    m = (q - r) / 2
    n = (q + r) / 2
    require(m.denominator == n.denominator == 1, "decoded core is nonintegral")
    if imbalance > 0:
        return Fraction(m, n)
    return Fraction(n, m)


def recover_row_from_reciprocal_owner_tree(speeds):
    """Recover a primitive mixed-parity row from N-1 cross-v2 edge summaries."""
    speeds = tuple(speeds)
    require(gcd_many(speeds) == 1, "primitive row required")
    levels = [valuation_two(speed) for speed in speeds]
    require(len(set(levels)) >= 2, "at least two 2-adic levels required")
    root = 0
    bridge = next(index for index, level in enumerate(levels) if level != levels[root])
    edges = []
    for index in range(1, len(speeds)):
        if levels[index] != levels[root]:
            edges.append((root, index))
        else:
            edges.append((bridge, index))
    require(len(edges) == len(speeds) - 1, "tree edge count changed")
    require(all(levels[i] != levels[j] for i, j in edges), "blind edge entered tree")

    adjacency = [[] for _ in speeds]
    for i, j in edges:
        maximum, measure = labelled_pair_summary(speeds[i], speeds[j])
        ratio = decode_labelled_ratio(maximum, measure)
        adjacency[i].append((j, ratio))       # speed_i / speed_j
        adjacency[j].append((i, 1 / ratio))

    coordinates = {root: Fraction(1)}
    stack = [root]
    while stack:
        i = stack.pop()
        for j, ratio_i_over_j in adjacency[i]:
            if j in coordinates:
                continue
            coordinates[j] = coordinates[i] / ratio_i_over_j
            stack.append(j)
    require(len(coordinates) == len(speeds), "owner tree disconnected")
    denominator_lcm = 1
    for value in coordinates.values():
        denominator_lcm = denominator_lcm * value.denominator // gcd(denominator_lcm, value.denominator)
    integer_row = [int(coordinates[index] * denominator_lcm) for index in range(len(speeds))]
    common = gcd_many(integer_row)
    return tuple(value // common for value in integer_row), tuple(edges)


def gcd_many(values):
    result = 0
    for value in values:
        result = gcd(result, value)
    return result


def audit_primitive_pairs(limit):
    audited = 0
    mixed = 0
    odd = 0
    for m in range(1, limit + 1):
        for n in range(m + 1, limit + 1):
            if gcd(m, n) != 1:
                continue
            audited += 1
            actual_maximum, actual_maximizers = exact_pair_maximum(m, n)
            actual_measure, actual_ties = exact_owner_measure(m, n)
            maximum, measure, ties, maximizers = closed_pair_data(m, n)
            require(actual_maximum == maximum, f"maximum mismatch {(m, n)}")
            require(actual_maximizers == maximizers, f"maximizer mismatch {(m, n)}")
            require(actual_measure == measure, f"owner measure mismatch {(m, n)}")
            require(actual_ties == ties, f"tie count mismatch {(m, n)}")
            defect = Fraction(1, 2) - maximum
            imbalance = measure - Fraction(1, 2)
            require(defect == (n - m) * imbalance, f"reciprocity mismatch {(m, n)}")
            if imbalance:
                recovered_sum = Fraction(1, 2 * defect)
                recovered_gap = defect / imbalance
                recovered_m = (recovered_sum - recovered_gap) / 2
                recovered_n = (recovered_sum + recovered_gap) / 2
                require((recovered_m, recovered_n) == (m, n), f"recovery mismatch {(m, n)}")
                mixed += 1
            else:
                require(m % 2 == n % 2 == 1, f"unexpected zero imbalance {(m, n)}")
                odd += 1
    return audited, mixed, odd


def audit_scaled_pairs(limit):
    audited = 0
    repeated_levels = 0
    cross_levels = 0
    for u in range(1, limit + 1):
        for v in range(u + 1, limit + 1):
            audited += 1
            actual_maximum, actual_maximizers = exact_pair_maximum(u, v)
            actual_measure, actual_ties = exact_owner_measure(u, v)
            maximum, measure, ties, maximizers = closed_scaled_pair_data(u, v)
            require(actual_maximum == maximum, f"scaled maximum mismatch {(u, v)}")
            require(actual_maximizers == maximizers, f"scaled maximizer mismatch {(u, v)}")
            require(actual_measure == measure, f"scaled owner mismatch {(u, v)}")
            require(actual_ties == ties, f"scaled tie mismatch {(u, v)}")
            common = gcd(u, v)
            if valuation_two(u) == valuation_two(v):
                repeated_levels += 1
                require((maximum, measure) == (Fraction(1, 2), Fraction(1, 2)),
                        f"repeated-v2 edge not blind {(u, v)}")
                require(len(maximizers) == common, f"odd sheet count mismatch {(u, v)}")
            else:
                cross_levels += 1
                require(Fraction(1, 2) - maximum == Fraction(common, 2 * (u + v)),
                        f"scaled defect mismatch {(u, v)}")
                require(measure - Fraction(1, 2) == Fraction(common * common, 2 * (v * v - u * u)),
                        f"scaled imbalance mismatch {(u, v)}")
                require(len(maximizers) == 2 * common, f"mixed sheet count mismatch {(u, v)}")
    return audited, cross_levels, repeated_levels


def audit_square_wave_correlations(limit):
    """Independent common-grid audit of the Fourier normalization."""
    audited = 0
    for r in range(1, limit + 1):
        for s in range(1, limit + 1):
            audited += 1
            common = gcd(r, s)
            period_grid = r * s // common
            signed_sum = 0
            for cell in range(period_grid):
                midpoint_numerator = 2 * cell + 1
                sign_r = 1 if ((r * midpoint_numerator) // (2 * period_grid)) % 2 == 0 else -1
                sign_s = 1 if ((s * midpoint_numerator) // (2 * period_grid)) % 2 == 0 else -1
                signed_sum += sign_r * sign_s
            actual = Fraction(signed_sum, period_grid)
            rr, ss = r // common, s // common
            expected = Fraction(common * common, r * s) if rr % 2 == ss % 2 == 1 else Fraction(0)
            require(actual == expected, f"square-wave normalization mismatch {(r, s)}")
    return audited


def audit_exponent_pairs(limit):
    audited = 0
    unequal = 0
    equal = []
    same_odd = 0
    mixed = 0
    representative = []
    for a in range(2, limit + 1):
        for b in range(a + 1, limit + 1):
            audited += 1
            left_value = a**b
            right_value = b**a
            common = gcd(left_value, right_value)
            left, right, signed = exponent_core(a, b)
            require((left, right) == (left_value // common, right_value // common),
                    f"valuation core mismatch {(a, b)}")
            require(gcd(left, right) == 1, f"nonprimitive core {(a, b)}")
            swapped_left, swapped_right, swapped_signed = exponent_core(b, a)
            require((swapped_left, swapped_right) == (right, left), f"swap mismatch {(a, b)}")
            require(all(swapped_signed[p] == -signed[p] for p in set(signed) | set(swapped_signed)),
                    f"valuation sign mismatch {(a, b)}")
            if left == right:
                equal.append((a, b))
                require((left, right) == (1, 1), f"bad equality core {(a, b)}")
                continue
            unequal += 1
            m, n = sorted((left, right))
            q = m + n
            require(gcd(m * n, q) == 1, f"ABC triple not primitive {(a, b)}")
            delta_two = signed.get(2, 0)
            require((delta_two == 0) == (m % 2 == n % 2 == 1),
                    f"2-adic parity gate mismatch {(a, b)}")
            maximum, measure, _, _ = closed_pair_data(m, n)
            if delta_two == 0:
                same_odd += 1
                require((maximum, measure) == (Fraction(1, 2), Fraction(1, 2)),
                        f"odd shell mismatch {(a, b)}")
            else:
                mixed += 1
                defect = Fraction(1, 2) - maximum
                require(defect == Fraction(1, 2 * q), f"sum-clock mismatch {(a, b)}")
            if (a, b) in {(2, 3), (2, 5), (3, 5), (3, 6), (6, 12)}:
                representative.append((a, b, m, n, q, delta_two, maximum, measure))
    require(equal == [(2, 4)], f"unexpected exponent equalities through {limit}: {equal}")
    return audited, unequal, equal, same_odd, mixed, representative


def audit_row_recovery(limit, size):
    primitive_mixed = 0
    all_odd_safe = 0
    for row in combinations(range(1, limit + 1), size):
        if gcd_many(row) != 1:
            continue
        if all(speed % 2 for speed in row):
            all_odd_safe += 1
            require(all(distance_to_integer(Fraction(speed, 2)) == Fraction(1, 2)
                        for speed in row), f"all-odd antipode failed {row}")
            continue
        primitive_mixed += 1
        recovered, edges = recover_row_from_reciprocal_owner_tree(row)
        require(recovered == row, f"row recovery mismatch {row} -> {recovered}")
        require(len(edges) == size - 1, f"nonminimal tree {row}")
    return primitive_mixed, all_odd_safe


def main():
    pair_data = audit_primitive_pairs(80)
    scaled_data = audit_scaled_pairs(60)
    correlation_data = audit_square_wave_correlations(40)
    exponent_data = audit_exponent_pairs(18)
    row_data = audit_row_recovery(16, 4)

    odd_a = closed_pair_data(1, 3)
    odd_b = closed_pair_data(1, 5)
    require(odd_a[:2] == odd_b[:2] == (Fraction(1, 2), Fraction(1, 2)),
            "odd-shell hostile summaries differ")
    hostile_a = exact_multi_maximum((1, 2, 3))
    hostile_b = exact_multi_maximum((1, 2, 5))
    require(hostile_a == (Fraction(1, 4), (Fraction(1, 4), Fraction(3, 4))),
            "AP3 hostile changed")
    require(hostile_b == (Fraction(1, 3), (Fraction(1, 3), Fraction(2, 3))),
            "odd-pair extension hostile changed")

    print("LRC EXPONENT/OWNER RECIPROCITY EXACT AUDIT")
    print(f"primitive pairs m<n<=80: total={pair_data[0]}, mixed={pair_data[1]}, odd={pair_data[2]}")
    print("all exact pair maxima, complete maximizer sets, owner-cell measures, tie counts, and mixed-pair recovery: PASS")
    print(f"all scaled pairs 0<u<v<=60: total={scaled_data[0]}, cross-v2={scaled_data[1]}, repeated-v2={scaled_data[2]}")
    print("all scaled equality clocks, complete lifted maximizers, owner measures, and v2-level conditions: PASS")
    print(f"independent square-wave correlations 1<=r,s<=40: total={correlation_data}, Fourier normalization: PASS")
    print(f"exponent pairs 2<=a<b<=18: total={exponent_data[0]}, unequal={exponent_data[1]}, equal={exponent_data[2]}")
    print(f"exponent-core parity: same-odd={exponent_data[3]}, mixed={exponent_data[4]}")
    for row in exponent_data[5]:
        a, b, m, n, q, delta_two, maximum, measure = row
        print(f"  (a,b)=({a},{b}) core=({m},{n}) q={q} Delta2={delta_two} M={maximum} owner={measure}")
    print("hostile pair summaries: {1,3} and {1,5} both (M,owner)=(1/2,1/2)")
    print(f"hostile extensions: M(1,2,3)={hostile_a[0]}, M(1,2,5)={hostile_b[0]}")
    print(f"primitive 4-rows in [1,16]: reconstructed mixed={row_data[0]}, all-odd antipode controls={row_data[1]}")
    print("VERDICT: exact pair reconstruction holds only in the mixed-parity shell; owner locations remain mandatory for multi-speed LRC.")


if __name__ == "__main__":
    main()
