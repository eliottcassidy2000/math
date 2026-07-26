#!/usr/bin/env python3
"""Exact companion for THM-2402.

The executable checks rational one-sheet and excess examples, signed
boundary pushforward, the common-core arithmetic, and the septimal
endpoint congruence.  It uses only the Python standard library and keeps
all assertions active under ``python -O`` by calling ``require``.
"""

from collections import defaultdict
from fractions import Fraction


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def mod_one(x):
    return x % 1


def boundary(intervals):
    """Distributional boundary of reduced nonwrapping half-open intervals."""
    out = defaultdict(int)
    for left, right in intervals:
        require(Fraction(0) <= left < right <= Fraction(1), "bad interval")
        out[mod_one(left)] += 1
        out[mod_one(right)] -= 1
    return {x: c for x, c in out.items() if c}


def push_boundary(q, signed):
    out = defaultdict(int)
    for x, coefficient in signed.items():
        out[mod_one(q * x)] += coefficient
    return {z: c for z, c in out.items() if c}


def in_intervals(x, intervals):
    x = mod_one(x)
    return any(left <= x < right for left, right in intervals)


def root_count(q, intervals, z):
    return sum(
        in_intervals((z + j) / q, intervals)
        for j in range(q)
    )


def verify_piecewise_counts(q, intervals, expected, critical):
    points = sorted(set(Fraction(x) for x in critical) | {Fraction(0), Fraction(1)})
    for left, right in zip(points, points[1:]):
        if left == right:
            continue
        z = (left + right) / 2
        require(root_count(q, intervals, z) == expected(z), "root-count mismatch")


def fmt_signed(signed):
    return ", ".join(
        f"{x}:{coefficient:+d}" for x, coefficient in sorted(signed.items())
    ) or "0"


def main():
    print("THM-2402 ORBIT EQUALITY AND ENDPOINT-PAIRING AUDIT")
    print()

    # A nonconstant one-sheet section over H=[0,1/2).
    q = 5
    H = [(Fraction(0), Fraction(1, 2))]
    E = [
        (Fraction(0), Fraction(1, 20)),
        (Fraction(1, 4), Fraction(3, 10)),
    ]
    mu_H = Fraction(1, 2)
    mu_E = sum(right - left for left, right in E)
    require(q * mu_E == mu_H, "one-sheet mass identity")
    verify_piecewise_counts(
        q,
        E,
        lambda z: int(Fraction(0) <= z < Fraction(1, 2)),
        [Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(1)],
    )
    dE = boundary(E)
    dH = boundary(H)
    pushed_dE = push_boundary(q, dE)
    require(pushed_dE == dH, "one-sheet signed-boundary pushforward")
    require(
        push_boundary(q, {Fraction(1, 20): -1, Fraction(1, 4): 1}) == {},
        "interior opposite-sign pairing",
    )

    print("NONCONSTANT ONE-SHEET SECTION (q=5)")
    print(f"  mu(E)={mu_E}, q*mu(E)=mu(H)={mu_H}")
    print(f"  D(E)={fmt_signed(dE)}")
    print(f"  (T_q)_*D(E)=D(H)={fmt_signed(dH)}")
    print("  interior switch z=1/4 pairs exit x=1/20 with entry x=1/4")
    print()

    # Upper multiplicity-support bound is sharp: exactly one extra sheet.
    E_upper = E + [(Fraction(17, 40), Fraction(9, 20))]
    mu_E_upper = sum(right - left for left, right in E_upper)
    excess_upper = q * mu_E_upper - mu_H
    B_upper = Fraction(1, 8)
    require(excess_upper == B_upper, "upper support bound should be sharp")
    verify_piecewise_counts(
        q,
        E_upper,
        lambda z: (
            2
            if Fraction(1, 8) <= z < Fraction(1, 4)
            else int(Fraction(0) <= z < Fraction(1, 2))
        ),
        [
            Fraction(0),
            Fraction(1, 8),
            Fraction(1, 4),
            Fraction(1, 2),
            Fraction(1),
        ],
    )
    defect_boundary = {
        z: push_boundary(q, boundary(E_upper)).get(z, 0) - dH.get(z, 0)
        for z in set(push_boundary(q, boundary(E_upper))) | set(dH)
    }
    defect_boundary = {z: c for z, c in defect_boundary.items() if c}
    require(
        defect_boundary
        == {Fraction(1, 8): 1, Fraction(1, 4): -1},
        "strict defect boundary",
    )

    # Lower multiplicity-support bound is sharp: all five sheets above J.
    J = (Fraction(0), Fraction(1, 8))
    E_lower = [(Fraction(0), Fraction(1, 5))]
    E_lower += [
        ((Fraction(j) + J[0]) / q, (Fraction(j) + J[1]) / q)
        for j in range(1, q)
    ]
    mu_E_lower = sum(right - left for left, right in E_lower)
    excess_lower = q * mu_E_lower - 1
    B_lower = J[1] - J[0]
    require(excess_lower == (q - 1) * B_lower, "lower support bound sharp")
    verify_piecewise_counts(
        q,
        E_lower,
        lambda z: q if J[0] <= z < J[1] else 1,
        [Fraction(0), Fraction(1, 8), Fraction(1)],
    )

    print("SHARP EXCESS-SUPPORT BOUNDS")
    print(
        "  one extra sheet: "
        f"X={excess_upper}, mu(B_2)={B_upper}=X"
    )
    print(
        "  all five sheets: "
        f"X={excess_lower}, mu(B_2)={B_lower}=X/(q-1)"
    )
    print(f"  strict-example D(N_E-1_H)={fmt_signed(defect_boundary)}")
    print()

    # THM-2396 common-core specialization.
    common_core = Fraction(1) - 2 * Fraction(1, 7) + Fraction(1, 91)
    floor = common_core / 49
    require(common_core == Fraction(66, 91), "common-core mass")
    require(floor == Fraction(66, 4459), "common-core clean floor")
    require(floor == Fraction(396, 26754), "universal-denominator form")

    print("COMMON-CORE q=49 SPECIALIZATION")
    print(f"  mu(H_0)=1-2/7+1/91={common_core}")
    print(f"  equality floor delta=mu(H_0)/49={floor}")
    print(
        "  strict surplus X=49*delta-66/91 obeys "
        "X/48 <= mu(B_2) <= X"
    )
    print()

    # Exhaust the elementary septimal endpoint congruence on a finite
    # positive-control box.  The proof in the theorem is algebraic and
    # unbounded; this scan catches sign, normalization, and guard-width bugs.
    tested_pairs = 0
    matched_pairs = 0
    unit_matches = 0
    primitive_unit_failures = 0
    for v in range(1, 21):
        for w in range(1, 21):
            for L in (1, 2):
                for M in (1, 2):
                    for a in range(v):
                        for b in range(w):
                            for eps in (-1, 1):
                                for eta in (-1, 1):
                                    tested_pairs += 1
                                    x = Fraction(14 * a + eps * L, 14 * v)
                                    y = Fraction(14 * b + eta * M, 14 * w)
                                    displacement = 49 * (x - y)
                                    if displacement.denominator != 1:
                                        continue
                                    matched_pairs += 1
                                    m = displacement.numerator
                                    lhs = 49 * (
                                        w * (14 * a + eps * L)
                                        - v * (14 * b + eta * M)
                                    )
                                    rhs = 14 * m * v * w
                                    require(lhs == rhs, "endpoint congruence")
                                    if v % 7 and w % 7:
                                        unit_matches += 1
                                        if m % 7:
                                            primitive_unit_failures += 1
    require(primitive_unit_failures == 0, "septimal unit conclusion")

    # Sharp boundary: one seven-divisible speed permits primitive displacement.
    v, w, L, M, a, b, eps, eta = 7, 1, 1, 1, 0, 0, -1, -1
    raw_x = Fraction(14 * a + eps * L, 14 * v)
    raw_y = Fraction(14 * b + eta * M, 14 * w)
    x, y = mod_one(raw_x), mod_one(raw_y)
    m = 49 * (raw_x - raw_y)
    require(x == Fraction(97, 98) and y == Fraction(13, 14), "hostile endpoints")
    require(m == 3 and m % 7, "primitive hostile displacement")

    print("SEPTIMAL ENDPOINT CONGRUENCE")
    print(f"  finite positive-control endpoint pairs tested={tested_pairs}")
    print(f"  common T_49 image pairs={matched_pairs}")
    print(f"  both speeds seven-units={unit_matches}")
    print("  primitive unit-speed failures=0")
    print(
        "  sharp hostile: "
        "v=7,w=1,x=97/98,y=13/14,49(x-y)=3"
    )
    print()

    print("SCOPE")
    print("  measurable equality does not require BV or rational endpoints")
    print("  signed endpoint pairing requires a finite-interval/BV boundary")
    print("  THM-2399 is pointwise one-orbit sharpness, not global equality")
    print("  primitive paired-event existence and canonical target typing remain open")
    print()
    print("PASS")


if __name__ == "__main__":
    main()
