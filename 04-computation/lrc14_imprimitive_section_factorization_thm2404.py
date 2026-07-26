#!/usr/bin/env python3
"""Exact companion for THM-2404.

Checks the two-stage section factorization, imprimitive/primitive
intermediate-boundary controls, and the sharp common-core length gate.
Only standard-library exact rational arithmetic is used.
"""

from collections import defaultdict
from fractions import Fraction


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def mod_one(x):
    return x % 1


def boundary(intervals):
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
    return sum(in_intervals((z + j) / q, intervals) for j in range(q))


def in_danger(v, x):
    residue = mod_one(v * x)
    return min(residue, 1 - residue) < Fraction(1, 14)


def check_section(q, E, H, critical):
    points = sorted(set(critical) | {Fraction(0), Fraction(1)})
    for left, right in zip(points, points[1:]):
        if left == right:
            continue
        z = (left + right) / 2
        expected = int(any(a <= z < b for a, b in H))
        require(root_count(q, E, z) == expected, "section count mismatch")


def fmt_boundary(signed):
    return ", ".join(
        f"{x}:{coefficient:+d}" for x, coefficient in sorted(signed.items())
    ) or "0"


def main():
    print("THM-2404 IMPRIMITIVE SECTION FACTORIZATION AUDIT")
    print()

    H = [(Fraction(0), Fraction(1, 2))]
    cut = Fraction(1, 4)

    # Sheet labels 9 and 23 have the same residue 2 modulo seven.
    E_imprimitive = [
        (Fraction(9, 49), Fraction(37, 196)),
        (Fraction(93, 196), Fraction(47, 98)),
    ]
    Y_imprimitive = [(Fraction(2, 7), Fraction(5, 14))]
    check_section(
        49,
        E_imprimitive,
        H,
        [Fraction(0), cut, Fraction(1, 2), Fraction(1)],
    )
    check_section(
        7,
        Y_imprimitive,
        H,
        [Fraction(0), Fraction(1, 2), Fraction(1)],
    )
    check_section(
        7,
        E_imprimitive,
        Y_imprimitive,
        [
            Fraction(0),
            Fraction(2, 7),
            Fraction(9, 28),
            Fraction(5, 14),
            Fraction(1),
        ],
    )
    require(23 - 9 == 14 and (23 - 9) % 7 == 0, "imprimitive switch")
    dY_imprimitive = boundary(Y_imprimitive)
    require(cut not in push_boundary(7, dY_imprimitive), "no interior quotient boundary")
    require(
        all(not in_intervals(y, [(Fraction(0), Fraction(1, 14)),
                                 (Fraction(13, 14), Fraction(1))])
            for y in (Fraction(2, 7), Fraction(9, 28), Fraction(5, 14))),
        "intermediate q-star safety control",
    )

    print("IMPRIMITIVE SWITCH CONTROL")
    print("  T_49 sheet labels: 9 -> 23, displacement=14")
    print("  intermediate residue: 2 -> 2 mod 7")
    print("  Y=[2/7,5/14) has no boundary above z=1/4")
    print(f"  D(Y)={fmt_boundary(dY_imprimitive)}")
    print()

    # Labels 9 and 10 differ modulo seven.  The intermediate section
    # acquires an exit/entry pair over the same interior base point.
    E_primitive = [
        (Fraction(9, 49), Fraction(37, 196)),
        (Fraction(41, 196), Fraction(3, 14)),
    ]
    Y_primitive = [
        (Fraction(2, 7), Fraction(9, 28)),
        (Fraction(13, 28), Fraction(1, 2)),
    ]
    check_section(
        49,
        E_primitive,
        H,
        [Fraction(0), cut, Fraction(1, 2), Fraction(1)],
    )
    check_section(
        7,
        Y_primitive,
        H,
        [Fraction(0), cut, Fraction(1, 2), Fraction(1)],
    )
    check_section(
        7,
        E_primitive,
        Y_primitive,
        [
            Fraction(0),
            Fraction(2, 7),
            Fraction(9, 28),
            Fraction(13, 28),
            Fraction(1, 2),
            Fraction(1),
        ],
    )
    require(10 - 9 == 1 and (10 - 9) % 7, "primitive switch")
    dY_primitive = boundary(Y_primitive)
    interior_pair = {
        Fraction(9, 28): dY_primitive[Fraction(9, 28)],
        Fraction(13, 28): dY_primitive[Fraction(13, 28)],
    }
    require(interior_pair == {Fraction(9, 28): -1, Fraction(13, 28): 1},
            "primitive quotient signs")
    require(push_boundary(7, interior_pair) == {}, "opposite quotient pairing")

    print("PRIMITIVE SWITCH CONTROL")
    print("  T_49 sheet labels: 9 -> 10, displacement=1")
    print("  intermediate residue: 2 -> 3 mod 7")
    print("  quotient exit 9/28 and entry 13/28 both lie over z=1/4")
    print(f"  interior D(Y)={fmt_boundary(interior_pair)}")
    print()

    # Normalized common core D_1^c intersect D_13^c consists of eleven
    # equal rational gaps.  Pullback by V gives 11V components.
    core_components = [
        (Fraction(14 * r + 1, 182), Fraction(14 * r + 13, 182))
        for r in range(1, 12)
    ]
    require(len(core_components) == 11, "core component count")
    require(
        all(right - left == Fraction(6, 91) for left, right in core_components),
        "core component length",
    )
    require(
        sum(right - left for left, right in core_components) == Fraction(66, 91),
        "core mass",
    )
    critical = {Fraction(0), Fraction(1), Fraction(1, 14), Fraction(13, 14)}
    for k in range(13):
        critical.add(mod_one(Fraction(14 * k - 1, 182)))
        critical.add(mod_one(Fraction(14 * k + 1, 182)))
    critical = sorted(critical)
    for left, right in zip(critical, critical[1:]):
        if left == right:
            continue
        x = (left + right) / 2
        formula_safe = in_intervals(x, core_components)
        direct_safe = not in_danger(1, x) and not in_danger(13, x)
        require(formula_safe == direct_safe, "normalized core reconstruction")

    V = 1
    u = 91
    lifted_length = Fraction(6, 637 * V)
    unit_safe_gap = Fraction(6, 7 * u)
    require(lifted_length == unit_safe_gap, "sharp length boundary")

    # At u=91V every branch-s lift of normalized core gap r is exactly
    # the D_u-safe gap indexed by k=r+13s.
    exact_fits = 0
    for r, (left, right) in enumerate(core_components, start=1):
        for s in range(7):
            lift = ((left + s) / 7, (right + s) / 7)
            k = r + 13 * s
            safe_gap = (
                Fraction(14 * k + 1, 14 * u),
                Fraction(14 * k + 13, 14 * u),
            )
            require(lift == safe_gap, "sharp aligned gap fit")
            exact_fits += 1
    require(exact_fits == 77, "sharp fit census")

    # One unit beyond the boundary makes every lifted component longer
    # than every connected safe gap, independently of alignment.
    u_hostile = 92
    require(
        Fraction(6, 637 * V) > Fraction(6, 7 * u_hostile),
        "strict length obstruction",
    )
    require(7 * u_hostile > 637 * V, "q-star/c3 comparison")

    print("COMMON-CORE LENGTH GATE")
    print("  Hbar has 11 components, each length 6/91, total 66/91")
    print("  H_0=T_V^-1(Hbar): 11V components, each length 6/(91V)")
    print("  one T_7 root lift length=6/(637V)")
    print("  every D_u-safe component length=6/(7u)")
    print("  imprimitive equality therefore forces u<=91V, i.e. q_*<=c_3")
    print("  boundary u=91V: exact aligned fits=77 for V=1")
    print("  hostile u=92,V=1: q_*=644>c_3=637, so primitive switch is forced")
    print()

    # A primitive C7 boundary has no automatic C13 target content.
    # Tensoring with a root-constant thirteen-vector preserves the boundary
    # but its centered/nonzero-character part is identically zero.
    target_gate = [Fraction(1) for _ in range(13)]
    target_mean = sum(target_gate) / 13
    centered_target = [value - target_mean for value in target_gate]
    require(centered_target == [Fraction(0)] * 13, "root-constant target hostile")

    print("CROSS-PRIME TARGET HOSTILE")
    print("  primitive C_7 boundary tensor root-constant C_13 gate")
    print("  signed boundary survives; all 12 nonzero target colours vanish")
    print("  a nonconstant lawful 7x13 boundary table is an additional sidecar")
    print()

    print("SCOPE")
    print("  every one-sheet T_49 section factors through Y=T_7(E)")
    print("  imprimitive iff the intermediate sheet is componentwise constant")
    print("  q_*=7u safety forces Y subset D_u^c")
    print("  q_*>c_3 forces a primitive q_*-typed switch only on the equality face")
    print("  canonical owner/target transport and LRC(14) remain open")
    print()
    print("PASS")


if __name__ == "__main__":
    main()
