from fractions import Fraction as F
from math import gcd
from functools import reduce
from hashlib import sha256


R_D = F(1, 14)
R_E = F(1, 7)


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def frac(x):
    return x - (x.numerator // x.denominator)


def b2(x):
    u = frac(x)
    return u * u - u + F(1, 6)


def equal_formula(p, q):
    return F(1, 49) + (b2(F(q - p, 14)) - b2(F(q + p, 14))) / (p * q)


def mixed_formula(p, q):
    return F(2, 49) + (b2(F(q - 2 * p, 14)) - b2(F(q + 2 * p, 14))) / (p * q)


def is_danger(s, x, radius):
    u = frac(s * x)
    return min(u, 1 - u) < radius


def walls_for_frequency(s, radius):
    ans = {F(0), F(1)}
    for k in range(s):
        ans.add(F(k, s) + radius / s)
        ans.add(F(k + 1, s) - radius / s)
    return ans


def integrate(specs, predicate):
    walls = {F(0), F(1)}
    for s, radius in specs:
        walls.update(walls_for_frequency(s, radius))
    walls = sorted(walls)
    total = F(0)
    for a, b in zip(walls, walls[1:]):
        if predicate((a + b) / 2):
            total += b - a
    return total


def literal_equal(p, q):
    return integrate(
        ((p, R_D), (q, R_D)),
        lambda x: is_danger(p, x, R_D) and is_danger(q, x, R_D),
    )


def literal_mixed(p, q):
    return integrate(
        ((p, R_D), (q, R_E)),
        lambda x: is_danger(p, x, R_D) and is_danger(q, x, R_E),
    )


def ratio_name(pair):
    p, q = pair
    return str(p) if q == 1 else f"{p}/{q}"


def enumerate_level(allowed_q, level, cutoff, formula):
    below = []
    equal = []
    checks = 0
    for p in range(1, cutoff + 1):
        for q in range(1, cutoff // p + 1):
            if p == q or gcd(p, q) != 1 or not allowed_q(q):
                continue
            value = formula(p, q)
            checks += 1
            if value < level:
                below.append((p, q))
            elif value == level:
                equal.append((p, q))
    key = lambda z: F(z[0], z[1])
    return sorted(below, key=key), sorted(equal, key=key), checks


def reduced_ratio(c, r):
    d = gcd(c, r)
    return c // d, r // d


def gcd_list(xs):
    return reduce(gcd, xs)


def safe_measure(speeds):
    specs = tuple((s, R_D) for s in speeds)
    return integrate(specs, lambda x: all(not is_danger(s, x, R_D) for s in speeds))


def union_overlap_with_r(C, r):
    specs = tuple((s, R_D) for s in tuple(C) + (r,))
    return integrate(
        specs,
        lambda x: is_danger(r, x, R_D) and any(is_danger(c, x, R_D) for c in C),
    )


def safe_intersect_E(C, r):
    specs = tuple((s, R_D) for s in C) + ((r, R_E),)
    return integrate(
        specs,
        lambda y: all(not is_danger(c, y, R_D) for c in C) and is_danger(r, y, R_E),
    )


def affine_danger(s, lift, y):
    return is_danger(1, F(s, 2) * (y + lift), R_D)


def affine_walls(s, lift):
    ans = {F(0), F(1)}
    # s(y+lift)/2=n +/- 1/14.  The deliberately generous range is exact
    # after retaining only walls in the unit interval.
    for n in range(-s - 2, 2 * s + 3):
        for sign in (-1, 1):
            y = F(2, s) * (n + sign * R_D) - lift
            if 0 < y < 1:
                ans.add(y)
    return ans


def cross_comb_measure_and_components(p, q):
    walls = {F(0), F(1)}
    for s in (p, q):
        for lift in (0, 1):
            walls.update(affine_walls(s, lift))
    walls = sorted(walls)
    def covered(y):
        covered0 = affine_danger(p, 0, y) or affine_danger(q, 0, y)
        covered1 = affine_danger(p, 1, y) or affine_danger(q, 1, y)
        return covered0 and covered1

    live = []
    for a, b in zip(walls, walls[1:]):
        y = (a + b) / 2
        if covered(y):
            live.append((a, b))
    merged = []
    for a, b in live:
        if merged and merged[-1][1] == a and covered(a):
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return sum((b - a for a, b in merged), F(0)), tuple(merged)


def main():
    # Formula checks are independent literal wall integrations.
    eq_literal_checks = 0
    mixed_literal_checks = 0
    for p in range(1, 57):
        for q in range(1, 57 // p + 1):
            if gcd(p, q) == 1:
                need(equal_formula(p, q) == literal_equal(p, q), f"equal formula {p}/{q}")
                eq_literal_checks += 1
    for p in range(1, 50):
        for q in range(1, 50 // p + 1):
            if gcd(p, q) == 1:
                need(mixed_formula(p, q) == literal_mixed(p, q), f"mixed formula {p}/{q}")
                mixed_literal_checks += 1

    classes = (
        ("g=1", lambda q: gcd(q, 6) == 1, F(1, 63), 55, 56),
        ("g=2", lambda q: q % 3 != 0, F(1, 70), 40, 41),
        ("g=3", lambda q: q % 2 == 1, F(1, 70), 40, 41),
        ("g=6", lambda q: True, F(1, 77), 33, 34),
    )
    derived = {}
    for name, allowed, level, cutoff, first_excluded in classes:
        below, equal, checks = enumerate_level(allowed, level, cutoff, equal_formula)
        need(F(1, 49) - F(1, 4 * first_excluded) > level, f"cutoff {name}")
        derived[name] = (below, equal, checks)
        print(
            f"equal-radius {name}: level={level}; cutoff={cutoff}; "
            f"below={','.join(map(ratio_name, below))}; "
            f"equal={','.join(map(ratio_name, equal))}; enumerated={checks}"
        )

    expected = {
        "g=1": (
            ((1, 13), (1, 11), (2, 11), (3, 11), (10, 1), (11, 1), (12, 1), (13, 1)),
            ((9, 5), (9, 1), (27, 1)),
        ),
        "g=2": (
            ((1, 13), (1, 11), (2, 11), (3, 11), (11, 2), (11, 1), (12, 1), (13, 1)),
            ((1, 10), (3, 10), (10, 1)),
        ),
        "g=3": (
            ((1, 13), (1, 11), (2, 11), (3, 11), (11, 3), (11, 1), (12, 1), (13, 1)),
            ((10, 3), (10, 1)),
        ),
        "g=6": (
            ((1, 13), (1, 12), (12, 1), (13, 1)),
            ((1, 11), (2, 11), (3, 11), (11, 3), (11, 2), (11, 1)),
        ),
    }
    for name, (want_below, want_equal) in expected.items():
        got_below, got_equal, _ = derived[name]
        need(tuple(got_below) == want_below, f"below atlas {name}")
        need(tuple(got_equal) == want_equal, f"equality atlas {name}")

    controls = (
        ("g=1", 715, (55, 9295, 8580, 65, 130, 195, 7865, 7150, 1287, 6435), F(1, 63)),
        ("g=2", 1430, (110, 130, 260, 390, 7865, 15730, 17160, 18590, 143, 429), F(1, 70)),
        ("g=3", 429, (33, 39, 78, 117, 1573, 4719, 5148, 5577, 1430, 4290), F(1, 70)),
        ("g=6", 1716, (132, 143, 20592, 22308, 156, 312, 468, 6292, 9438, 18876), F(1, 77)),
    )
    for name, r, C, level in controls:
        need(len(C) == len(set(C)) == 10 and r not in C, f"control distinctness {name}")
        need(gcd_list(C) == 1, f"control primitivity {name}")
        vals = [equal_formula(*reduced_ratio(c, r)) for c in C]
        need(max(vals) == level and all(v <= level for v in vals), f"control sharpness {name}")
        print(
            f"control {name}: r={r}; gcdC={gcd_list(C)}; "
            f"below={sum(v < level for v in vals)}; equal={sum(v == level for v in vals)}; max={max(vals)}"
        )

    losses = {
        "g=1": F(1, 7) - F(1, 63),
        "g=2/3": F(1, 7) - F(1, 70),
        "g=6": F(1, 7) - F(1, 77),
    }
    need(losses == {"g=1": F(8, 63), "g=2/3": F(9, 70), "g=6": F(10, 77)}, "q2 losses")
    caps = (F(4, 63), F(4, 77), F(4, 91))
    expected_thresholds = {
        "g=1": (F(4, 21), F(124, 693), F(20, 117)),
        "g=2/3": (F(121, 630), F(139, 770), F(157, 910)),
        "g=6": (F(134, 693), F(2, 11), F(174, 1001)),
    }
    for name, loss in losses.items():
        got = tuple(loss + cap for cap in caps)
        need(got == expected_thresholds[name], f"q2 thresholds {name}")
        print(f"q2 {name}: loss={loss}; thresholds={','.join(map(str, got))}")

    # The five primitive cross-comb masses are reconstructed directly from
    # the two physical lifts, not imported as constants.
    five_rays = ((1, 11), (1, 23), (5, 11), (1, 37), (1, 25))
    expected_pair_caps = (F(4, 77), F(8, 161), F(18, 385), F(12, 259), F(8, 175))
    pair_caps = []
    for ray, want in zip(five_rays, expected_pair_caps):
        mass, components = cross_comb_measure_and_components(*ray)
        need(mass == want, f"cross-comb mass {ray}")
        pair_caps.append(mass)
        print(
            f"cross-comb {ray[0]}:{ray[1]}: mass={mass}; components={len(components)}; "
            f"widths={','.join(str(b-a) for a,b in components)}"
        )
    ratio_entries = tuple(F(8, 63) + x for x in pair_caps)
    expected_entries = (F(124, 693), F(256, 1449), F(86, 495), F(404, 2331), F(272, 1575))
    need(ratio_entries == expected_entries, "ratio-specific q2 entries")
    print(f"q2 five-ray entries={','.join(map(str, ratio_entries))}")

    # Literal set integrations test both transfer identities, including their
    # strict-boundary convention, on unrelated small bodies.
    identity_tests = (
        ((1, 4, 7), 3),
        ((2, 5, 9, 11), 7),
        ((1, 2, 3, 5), 11),
    )
    for C, r in identity_tests:
        gc = safe_measure(C)
        q2_left = safe_measure(tuple(C) + (r,))
        q2_right = gc - F(1, 7) + union_overlap_with_r(C, r)
        need(q2_left == q2_right, f"q2 identity C={C},r={r}")
        q4_left = safe_measure(tuple(2 * c for c in C) + (r,))
        q4_right = gc - F(1, 2) * safe_intersect_E(C, r)
        need(q4_left == q4_right, f"q4 identity C={C},r={r}")
        print(f"identities C={C},r={r}: q2={q2_left}; q4={q4_left}")

    mixed_level = F(1, 28)
    mixed_below, mixed_equal, mixed_checks = enumerate_level(
        lambda q: gcd(q, 6) == 1, mixed_level, 49, mixed_formula
    )
    need(F(2, 49) - F(1, 4 * 50) > mixed_level, "mixed cutoff")
    need(
        tuple(mixed_below) == ((1, 25), (1, 13), (1, 11), (2, 11), (5, 1), (6, 1), (13, 1)),
        "mixed below atlas",
    )
    need(tuple(mixed_equal) == ((4, 5), (4, 1), (12, 1), (20, 1)), "mixed equality atlas")
    print(
        f"mixed: level={mixed_level}; cutoff=49; below={','.join(map(ratio_name,mixed_below))}; "
        f"equal={','.join(map(ratio_name,mixed_equal))}; enumerated={mixed_checks}"
    )

    mixed_r = 3575
    mixed_C = (21450, 325, 17875, 650, 275, 46475, 143, 2860, 14300, 42900)
    need(len(mixed_C) == len(set(mixed_C)) == 10 and mixed_r not in mixed_C, "mixed control distinctness")
    need(gcd_list(mixed_C) == 1 and gcd(mixed_r, 6) == 1, "mixed control class/primitivity")
    mixed_vals = [mixed_formula(*reduced_ratio(c, mixed_r)) for c in mixed_C]
    need(max(mixed_vals) == mixed_level and all(v <= mixed_level for v in mixed_vals), "mixed control")
    print(
        f"mixed control: r={mixed_r}; gcdC={gcd_list(mixed_C)}; "
        f"below={sum(v < mixed_level for v in mixed_vals)}; "
        f"equal={sum(v == mixed_level for v in mixed_vals)}; max={max(mixed_vals)}"
    )

    # Transfer arithmetic from the mixed tenth statistic:
    # mu(E_r)=2/7 and union_C D_c captures at least 1/28 of E_r.
    e_mass = F(2, 7)
    max_gc_intersect_e = e_mass - mixed_level
    need(max_gc_intersect_e == F(1, 4), "q4 intersection cap")
    need(F(1, 2) * max_gc_intersect_e == F(1, 8), "q4 additive loss")
    print(
        f"q4 transfer: muE={e_mass}; GCcap={max_gc_intersect_e}; additive_loss={F(1,8)}; "
        "floor=max(muGC/2,muGC-1/8)"
    )

    print(
        f"status=PASS;equal_literal_checks={eq_literal_checks};"
        f"mixed_literal_checks={mixed_literal_checks};q2_identity_checks={len(identity_tests)};"
        f"q4_identity_checks={len(identity_tests)}"
    )


if __name__ == "__main__":
    main()
