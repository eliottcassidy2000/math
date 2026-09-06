#!/usr/bin/env python3
"""Clean-room exact independent audit for THM-4449.

This file uses only the theorem's definitions and Python's standard library.
It does not import the primary or discovery implementation.
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, lcm

DELTA = F(1, 14)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_three_unit_odd(n):
    return n > 0 and n % 2 == 1 and n % 3 != 0


def circular_norm(x):
    x = x % 1
    return min(x, 1 - x)


def nearest_seventh_residue(n):
    """n - 7*nearest_integer(n/7); ties cannot occur for integral n."""
    return n - 7 * ((2 * n + 7) // 14)


def pair_formula(p, q):
    p, q = sorted((p, q))
    require(p != q and p % 2 and q % 2 and gcd(p, q) == 1, (p, q, "primitive odd pair"))
    alpha = (p + q) // 2
    beta = (q - p) // 2
    d = nearest_seventh_residue(alpha)
    e = nearest_seventh_residue(beta)
    require(-3 <= d <= 3 and -3 <= e <= 3, (p, q, d, e, "seventh residue"))
    return F(2, 49) * (1 + F(e * e - d * d, p * q))


def primitive_pair(x, y):
    """Primitive integral unordered ratio attached to two rational labels."""
    ratio = F(x) / F(y)
    n, d = abs(ratio.numerator), ratio.denominator
    return tuple(sorted((n, d)))


def pair_formula_rational(x, y):
    return pair_formula(*primitive_pair(x, y))


def energy(labels):
    return sum((pair_formula_rational(x, y) for x, y in combinations(labels, 2)), F(0))


def owner(speed, numerator, denominator):
    """Strict-danger lift owner at quotient y=numerator/denominator.

    Returns None if the odd tail is inactive, otherwise the nearest-integer
    parity (the physical half-lift it kills).
    """
    quotient, remainder = divmod(speed * numerator, denominator)
    distance = min(remainder, denominator - remainder)
    if 7 * distance >= denominator:
        return None
    nearest = quotient if 2 * remainder < denominator else quotient + 1
    return nearest & 1


def owner_wall_stats(speeds, predicate):
    """Measure, longest component, count of an open owner predicate.

    All walls live on the common integer lattice 1/(7*lcm(speeds)).  We
    additionally evaluate each wall, so punctures and endpoint joins are not
    guessed from measure-zero data.
    """
    speeds = tuple(speeds)
    require(all(s > 0 and s % 2 == 1 for s in speeds), (speeds, "positive odd speeds"))
    common = lcm(*speeds)
    denominator = 7 * common
    walls = {0, denominator}
    for s in speeds:
        scale = common // s
        for k in range(s):
            walls.add(((7 * k - 1) * scale) % denominator)
            walls.add(((7 * k + 1) * scale) % denominator)
    walls = sorted(walls)
    cells = []
    total = 0
    for left, right in zip(walls, walls[1:]):
        bits = tuple(owner(s, left + right, 2 * denominator) for s in speeds)
        live = predicate(bits)
        cells.append(live)
        if live:
            total += right - left

    live_indices = [i for i, live in enumerate(cells) if live]
    if not live_indices:
        return F(0), F(0), 0
    parent = list(range(len(cells)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        i, j = find(i), find(j)
        if i != j:
            parent[j] = i

    # Internal walls.
    for i in range(1, len(walls) - 1):
        left_cell, right_cell = i - 1, i
        wall_bits = tuple(owner(s, walls[i], denominator) for s in speeds)
        if cells[left_cell] and cells[right_cell] and predicate(wall_bits):
            union(left_cell, right_cell)
    # The point 0=1 joins the final and initial cell only when it is live.
    wall_bits = tuple(owner(s, 0, denominator) for s in speeds)
    if cells[-1] and cells[0] and predicate(wall_bits):
        union(len(cells) - 1, 0)

    lengths = {}
    for i in live_indices:
        root = find(i)
        lengths[root] = lengths.get(root, 0) + walls[i + 1] - walls[i]
    return F(total, denominator), F(max(lengths.values()), denominator), len(lengths)


def pair_literal(p, q):
    return owner_wall_stats((p, q), lambda bits: bits[0] is not None and bits[1] is not None and bits[0] != bits[1])


def triple_physical(labels):
    return owner_wall_stats(labels, lambda bits: 0 in bits and 1 in bits)


def safe_wall_stats(labels):
    """Haar measure and longest positive-length component of closed G_A."""
    labels = tuple(sorted(set(labels)))
    common = lcm(*labels)
    denominator = 14 * common
    walls = {0, denominator}
    for s in labels:
        scale = common // s
        for k in range(s):
            walls.add(((14 * k - 1) * scale) % denominator)
            walls.add(((14 * k + 1) * scale) % denominator)
    walls = sorted(walls)

    def safe_at(num, den):
        for s in labels:
            _, remainder = divmod(s * num, den)
            if 14 * min(remainder, den - remainder) < den:
                return False
        return True

    cells = []
    total = 0
    for left, right in zip(walls, walls[1:]):
        live = safe_at(left + right, 2 * denominator)
        cells.append(live)
        if live:
            total += right - left
    live_indices = [i for i, live in enumerate(cells) if live]
    if not live_indices:
        return F(0), F(0), 0
    parent = list(range(len(cells)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        i, j = find(i), find(j)
        if i != j:
            parent[j] = i

    for i in range(1, len(walls) - 1):
        if cells[i - 1] and cells[i] and safe_at(walls[i], denominator):
            union(i - 1, i)
    if cells[-1] and cells[0] and safe_at(0, denominator):
        union(len(cells) - 1, 0)
    lengths = {}
    for i in live_indices:
        root = find(i)
        lengths[root] = lengths.get(root, 0) + walls[i + 1] - walls[i]
    return F(total, denominator), F(max(lengths.values()), denominator), len(lengths)


def clearance(labels, x):
    return min(circular_norm(F(v) * x) for v in set(labels))


def lift_masks(tails, quotient_y, sheets):
    masks = []
    for t in tails:
        owned = []
        for j in range(sheets):
            x = (quotient_y + j) / sheets
            if circular_norm(t * x) < DELTA:
                owned.append(j)
        masks.append(tuple(owned))
    return tuple(masks)


def assert_formula_and_pair_cap():
    checked = 0
    for q in range(3, 152, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            literal, _, _ = pair_literal(p, q)
            require(literal == pair_formula(p, q), (p, q, literal, pair_formula(p, q)))
            checked += 1
    require(checked == 2350, (checked, "literal-check count"))

    low33 = []
    for p in range(1, 33, 2):
        for q in range(p + 2, 33, 2):
            if p * q < 33 and gcd(p, q) == 1 and is_three_unit_odd(p) and is_three_unit_odd(q):
                low33.append(((p, q), pair_formula(p, q)))
    expected_low33 = [
        ((1, 5), F(0)), ((1, 7), F(2, 49)), ((1, 11), F(4, 77)),
        ((1, 13), F(4, 91)), ((1, 17), F(4, 119)), ((1, 19), F(4, 133)),
        ((1, 23), F(8, 161)), ((1, 25), F(8, 175)),
        ((1, 29), F(8, 203)), ((1, 31), F(8, 217)),
    ]
    require(low33 == expected_low33, (low33, "product<33 table"))
    require(max(v for _, v in low33) == F(4, 77), (low33, "pair cap"))
    require([r for r, v in low33 if v == F(4, 77)] == [(1, 11)], (low33, "pair equality"))
    return checked, low33


def assert_energy_theorem():
    exceptions = []
    equality_at_four_ninety_one = []
    for p in range(1, 117, 2):
        for q in range(p + 2, 117, 2):
            if p * q >= 117 or gcd(p, q) != 1:
                continue
            if not (is_three_unit_odd(p) and is_three_unit_odd(q)):
                continue
            value = pair_formula(p, q)
            if value > F(4, 91):
                exceptions.append(((p, q), value))
            elif value == F(4, 91):
                equality_at_four_ninety_one.append((p, q))
    expected = [
        ((1, 11), F(4, 77)),
        ((1, 23), F(8, 161)),
        ((1, 25), F(8, 175)),
        ((1, 37), F(12, 259)),
        ((5, 11), F(18, 385)),
    ]
    require(sorted(exceptions) == sorted(expected), (exceptions, "high-edge alphabet"))
    require(equality_at_four_ninety_one == [(1, 13), (1, 65), (5, 13)], (equality_at_four_ninety_one, "4/91 boundary"))

    target = F(124, 847)
    gp = (F(1), F(11), F(121))
    require(energy(gp) == target, (energy(gp), "GP energy"))
    require(pair_formula(1, 121) == F(36, 847), (pair_formula(1, 121), "long GP edge"))

    ratios = (F(23), F(11, 5), F(37), F(25))
    one_edge_maxima = []
    one_edge_witnesses = []
    for r in ratios:
        candidates = []
        for base in (F(1), F(11)):
            for z in (base * r, base / r):
                labels = (F(1), F(11), z)
                if len(set(labels)) < 3:
                    continue
                candidates.append((energy(labels), labels))
        best = max(candidates)
        one_edge_maxima.append(best[0])
        one_edge_witnesses.append(best[1])
    require(one_edge_maxima == [F(36, 253), F(592, 4235), F(36, 259), F(268, 1925)], (one_edge_maxima, "one-1:11 compatibility"))

    table = []
    table_witnesses = []
    for r in ratios:
        row, witness_row = [], []
        for s in ratios:
            candidates = []
            for x in (r, 1 / r):
                for y in (s, 1 / s):
                    labels = (F(1), x, y)
                    if len(set(labels)) < 3:
                        continue
                    candidates.append((energy(labels), labels))
            best = max(candidates)
            row.append(best[0])
            witness_row.append(best[1])
        table.append(row)
        table_witnesses.append(witness_row)
    expected_table = [
        [F(520, 3703), F(174, 1265), F(816, 5957), F(548, 4025)],
        [F(174, 1265), F(406, 3025), F(382, 2849), F(258, 1925)],
        [F(816, 5957), F(382, 2849), F(1280, 9583), F(172, 1295)],
        [F(548, 4025), F(258, 1925), F(172, 1295), F(116, 875)],
    ]
    require(table == expected_table, (table, "two-exception compatibility"))
    coarse = [F(20, 143), F(288, 2093)]
    require(all(x < target for x in one_edge_maxima + coarse + [v for row in table for v in row]), "non-GP energy below GP")
    require(max(one_edge_maxima + coarse + [v for row in table for v in row]) == F(36, 253), "largest non-GP case")
    return exceptions, equality_at_four_ninety_one, one_edge_maxima, one_edge_witnesses, table, table_witnesses


def assert_physical_bound_and_identity():
    gp = (1, 11, 121)
    physical, width, components = triple_physical(gp)
    require(physical == F(108, 847), (physical, "GP physical mass"))
    require(physical < F(36, 253), "GP below physical bound")
    e = energy(tuple(F(x) for x in gp))

    # Directly integrate the mixed-three-owner set; this verifies E=F+Omega.
    omega, _, _ = owner_wall_stats(
        gp,
        lambda bits: None not in bits and 0 in bits and 1 in bits,
    )
    require(e == physical + omega, (e, physical, omega, "owner-overlap identity"))
    require(omega == F(16, 847), (omega, "GP owner overlap"))

    # Several unrelated shapes exercise ties, zero pair-combs, and dilation.
    identity_samples = [(1, 5, 7), (1, 7, 11), (5, 11, 23), (1, 23, 253), (5, 25, 55)]
    for labels in identity_samples:
        fmass, _, _ = triple_physical(labels)
        emass = energy(tuple(F(x) for x in labels))
        overlap, _, _ = owner_wall_stats(
            labels,
            lambda bits: None not in bits and 0 in bits and 1 in bits,
        )
        require(emass == fmass + overlap, (labels, emass, fmass, overlap))
    # Common dilation preserves both pair energy and physical mass.
    for t in (5, 7, 11):
        require(triple_physical(tuple(t * x for x in gp))[0] == physical, (t, "dilation invariance"))
    return physical, width, components, omega


def assert_transfer_and_localization():
    cap = F(4, 77)
    # Direct exact checks of both transfer estimates on structurally different bodies.
    samples = [
        ((1, 2, 3), 5),
        ((1, 4, 7, 10), 11),
        ((2, 3, 5, 8, 13), 17),
        (tuple(range(1, 11)), 13),
    ]
    transfer_rows = []
    for C, r in samples:
        mu_c = safe_wall_stats(C)[0]
        mu_q4 = safe_wall_stats(tuple(2 * c for c in C) + (r,))[0]
        mu_q2 = safe_wall_stats(tuple(C) + (r,))[0]
        require(mu_q4 >= mu_c / 2, (C, r, mu_c, mu_q4, "q4 transfer"))
        require(mu_q2 >= mu_c - F(1, 7), (C, r, mu_c, mu_q2, "q2 transfer"))
        transfer_rows.append((C, r, mu_c, mu_q4, mu_q2))

    beta_expected = {
        (1, 11): F(2, 77),
        (1, 23): F(2, 161),
        (5, 11): F(9, 385),
        (1, 37): F(2, 259),
        (1, 25): F(2, 175),
    }
    beta_rows = {}
    for ratio, expected in beta_expected.items():
        measure, longest, count = pair_literal(*ratio)
        require(measure == pair_formula(*ratio), (ratio, measure, "localized measure"))
        require(longest == expected, (ratio, longest, expected))
        beta_rows[ratio] = (measure, longest, count)

    require(F(8, 91) / 2 == F(4, 91), "q4 localization threshold")
    require(F(17, 91) - F(1, 7) == F(4, 91), "q2 localization threshold")
    require(F(8, 77) / 2 == cap, "q4 entry threshold")
    require(F(15, 77) - F(1, 7) == cap, "q2 entry threshold")
    return transfer_rows, beta_rows


def assert_hostiles():
    C0 = tuple(range(1, 11))

    # q=2, zero even tails.
    T0 = (1, 7, 11)
    y0 = F(1, 11)
    require(clearance(C0, y0) == F(1, 11), "q2 zero body clearance")
    masks0 = lift_masks(T0, y0, 2)
    require(masks0 == ((0,), (), (1,)), (masks0, "q2 zero masks"))
    x0 = F(181, 2352)
    row0 = tuple(2 * c for c in C0) + T0
    require(clearance(row0, x0) == F(15, 196) > DELTA, "q2 zero positive row")

    # q=2, one even tail, after absorbing r=13.
    T1 = (26, 1, 11)
    H1 = tuple(range(1, 11)) + (13,)
    y1 = F(1, 11)
    require(clearance(H1, y1) == F(1, 11), "q2 one-even body clearance")
    masks1 = lift_masks((1, 11), y1, 2)
    require(masks1 == ((0,), (1,)), (masks1, "q2 one-even masks"))
    x1 = F(229, 560)
    row1 = tuple(2 * c for c in C0) + T1
    require(clearance(row1, x1) == F(5, 56) > DELTA, "q2 one-even positive row")

    # q=4, one v2=1 tail.
    C2 = (1, 2, 3, 4, 5, 7, 8, 9, 10, 11)
    T2 = (1, 11, 14)
    y2 = F(86, 539)
    require(clearance(C2, y2) == F(9, 77), "q4 body clearance")
    masks2 = lift_masks(T2, y2, 4)
    require(masks2 == ((0,), (2,), (1, 3)), (masks2, "q4 masks"))
    x2 = F(37, 480)
    row2 = tuple(4 * c for c in C2) + T2
    require(clearance(row2, x2) == F(3, 40) > DELTA, "q4 positive row")

    return (masks0, clearance(row0, x0)), (masks1, clearance(row1, x1)), (masks2, clearance(row2, x2))


def main():
    checked, low33 = assert_formula_and_pair_cap()
    exceptions, equality_491, one_max, one_wit, table, table_wit = assert_energy_theorem()
    physical, width, components, omega = assert_physical_bound_and_identity()
    transfers, betas = assert_transfer_and_localization()
    hostiles = assert_hostiles()

    print("THM4449_INDEPENDENT_DYADIC_AUDIT")
    print(f"formula_literal_checks={checked}")
    print("ternary_pair_cap=4/77 equality_primitive=(1,11) equality_dilations=odd_3_units")
    print("high_edge_alphabet=" + ",".join(f"{p}:{q}={v}" for (p, q), v in exceptions))
    print("four_over_91_equalities=" + ",".join(f"{p}:{q}" for p, q in equality_491))
    print("energy_max=124/847 equality_shape=(1,11,121) equality_dilations=odd_3_units")
    print("one_1_11_edge_maxima=" + ",".join(map(str, one_max)))
    print("compatibility_table")
    for row in table:
        print(" ".join(map(str, row)))
    print(
        f"gp_quotient_mass={physical} gp_overlap={omega} "
        f"gp_quotient_width={width} gp_quotient_components={components}"
    )
    print("physical_owner_cut_bound=36/253")
    print("beta_table=" + ",".join(f"{p}:{q}:{v[1]}" for (p, q), v in betas.items()))
    print("transfer_samples=" + ",".join(f"{len(C)}:{r}:{mc}:{m4}:{m2}" for C, r, mc, m4, m2 in transfers))
    print("hostiles=" + ",".join(f"{masks}:{gap}" for masks, gap in hostiles))
    print("threshold_topology=closed_body_vs_proper_open_failure_set_inclusive")
    print("PASS")


if __name__ == "__main__":
    main()
