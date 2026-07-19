#!/usr/bin/env python3
"""Independent exact referee for the sharpened THM-1213 peel-first gate.

The generic 6/49 endpoint discrepancy in THM-1097 is sharp when the
fractional window length is unknown.  After peeling the first proportional
comb, however, the remaining q-comb sees the known scaled length 6q/(7a).
This referee keeps that fractional part and the integer tooth-incidence
floor.  It checks the resulting scale-free gate, the r >= 25/4 envelope,
the two equality cases, and an asymptotically sharp family immediately below
25/4.  All decisions use integers or fractions.Fraction.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import gcd


HEIGHT = 64
EXHAUSTIVE_CONE_HEIGHT = 400
CHAMPION_DENOMINATOR = 120
BOUNDARY_DENOMINATOR = 4096


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def floor_fraction(x: F) -> int:
    return x.numerator // x.denominator


def epsilon_fraction(q: int, a: int) -> F:
    """Exact excess over density 1/7 on a window of scaled length 6q/(7a)."""
    x = F(6 * q, 7 * a)
    theta = x - floor_fraction(x)
    return min(theta, F(1, 7)) - theta / 7


def epsilon_numerator(q: int, a: int) -> int:
    """Return E with epsilon=E/(49a), avoiding Fraction in large audits."""
    remainder = (6 * q) % (7 * a)
    if remainder <= a:
        return 6 * remainder
    return 7 * a - remainder


def incidence(q: int, a: int) -> int:
    """Phase-independent tooth count on a closed interval of length 6/(7ma)."""
    return (6 * q + 8 * a) // (7 * a)


def exact_gate_numerator(shape: tuple[int, int, int, int]) -> int:
    """Positive iff the exact mass/incidence component gate is strict."""
    a, b, c, d = shape
    qs = b, c, d
    product = b * c * d
    components = 1 + sum(incidence(q, a) for q in qs)
    # This is (new_mass_score-components) * 7*a*b*c*d.
    return (
        24 * d * product
        - d * sum(epsilon_numerator(q, a) * (product // q) for q in qs)
        - 7 * a * components * product
    )


def exact_gate_margin(shape: tuple[int, int, int, int]) -> F:
    a, b, c, d = shape
    mass_score = F(24 * d, 7 * a) - sum(
        F(7 * d, q) * epsilon_fraction(q, a) for q in (b, c, d)
    )
    return mass_score - 1 - sum(incidence(q, a) for q in (b, c, d))


def hybrid_q(shape: tuple[int, int, int, int]) -> F:
    a, b, c, d = shape
    return F(6 * (3 * d - b - c), a) - F(6 * d, b) - F(6 * d, c) - 37


def primitive(shape: tuple[int, int, int, int]) -> bool:
    a, b, c, d = shape
    return gcd(gcd(a, b), gcd(c, d)) == 1


def phi_pair(numerator: int, denominator: int) -> tuple[int, int] | None:
    if numerator < denominator:
        return None
    if 7 * numerator >= 13 * denominator:
        return 6, 7
    out_numerator = 7 * numerator - denominator
    out_denominator = 14 * denominator
    common = gcd(out_numerator, out_denominator)
    return out_numerator // common, out_denominator // common


def exact_transfer_gate(shape: tuple[int, int, int, int]) -> bool:
    numerator, denominator = 6, 7
    for old, new in zip(shape, shape[1:]):
        output = phi_pair(new * numerator, old * denominator)
        if output is None:
            return False
        numerator, denominator = output
    return 7 * numerator > denominator


def thm1148_residual(shape: tuple[int, int, int, int]) -> bool:
    a, b, c, d = shape
    return not (8 * a > 7 * d or 2 * d > a + b + c or exact_transfer_gate(shape))


def local_formula_audit() -> int:
    rows = 0
    for a in range(1, 301):
        for q in range(a + 1, 10 * a + 1):
            epsilon = epsilon_fraction(q, a)
            require(epsilon == F(epsilon_numerator(q, a), 49 * a), "epsilon formula")
            require(0 <= epsilon <= F(6, 49), "sharp discrepancy range")
            x = F(6 * q, 7 * a)
            # A closed centre interval of length x+1/7 contains at most
            # floor(x+1/7)+1 integers.  The +1 is exactly the 8/7 formula,
            # including the case in which x+1/7 is itself an integer.
            require(incidence(q, a) == floor_fraction(x + F(1, 7)) + 1, "incidence floor")
            rows += 1
    return rows


def height64_census() -> dict[str, object]:
    residual = 0
    old_gate = 0
    exact_gate = 0
    clean_cone = 0
    tightest: tuple[F, tuple[int, int, int, int]] | None = None
    for d in range(4, HEIGHT + 1):
        for a, b, c in combinations(range(1, d), 3):
            shape = a, b, c, d
            if not primitive(shape) or not thm1148_residual(shape):
                continue
            residual += 1
            old = hybrid_q(shape) > 0
            numerator = exact_gate_numerator(shape)
            exact = numerator > 0
            require(exact == (exact_gate_margin(shape) > 0), "cleared exact gate")
            if old:
                old_gate += 1
                require(exact, "exact gate must dominate Q_hyb")
            if exact:
                exact_gate += 1
                row = exact_gate_margin(shape), shape
                tightest = row if tightest is None else min(tightest, row)
            if 4 * d >= 25 * a:
                clean_cone += 1
                require(2 * d <= a + b + c, "cone row left Q4 residual")
                require(exact, "25/4 cone failure in height-64 census")

    require(residual == 95_336, "height-64 residual census moved")
    require(old_gate == 484, "old Q_hyb count moved")
    require(exact_gate == 4_028, "exact discrepancy count moved")
    require(clean_cone == 920, "25/4 cone count moved")
    require(tightest == (F(1, 1120), (11, 32, 35, 39)), "tight exact row moved")
    return {
        "residual": residual,
        "old": old_gate,
        "exact": exact_gate,
        "clean": clean_cone,
        "tightest": tightest,
    }


def kappa(q: int, a: int, d: int) -> F:
    """Sawtooth cost kappa_r(z), with r=d/a and z=q/a."""
    return F(incidence(q, a)) + F(d * epsilon_numerator(q, a), 7 * a * q)


def kappa_at_d(a: int, d: int) -> F:
    return kappa(d, a, d)


def champion_bound(a: int, d: int) -> F:
    """The exact upper-envelope formula on 25/4 <= d/a < 49/6."""
    r = F(d, a)
    require(F(25, 4) <= r < F(49, 6), "champion domain")
    if r < F(41, 6):
        return 6 + r / 7
    if r <= F(246, 35):
        return 7 + F(6, 287) * r
    if r < F(43, 6):
        return kappa_at_d(a, d)
    if r < 8:
        return 7 + F(36, 301) * r
    return 8 + r / 56


def champion_lattice_audit() -> int:
    rows = 0
    for a in range(1, CHAMPION_DENOMINATOR + 1):
        d_min = (25 * a + 3) // 4
        d_max = (49 * a - 1) // 6
        for d in range(d_min, d_max + 1):
            if not (25 * a <= 4 * d and 6 * d < 49 * a):
                continue
            champion = champion_bound(a, d)
            for q in range(d - a + 1, d):
                require(kappa(q, a, d) <= champion, "champion envelope failure")
                rows += 1
    return rows


def sufficient_margin_pieces() -> tuple[tuple[str, F, F], ...]:
    """Return each affine lower bound at the two ends of its domain."""
    pieces = (
        ("[25/4,41/6]", F(25, 4), F(41, 6), lambda r: 4 * r - 25),
        ("[41/6,7]", F(41, 6), F(7), lambda r: F(1218, 287) * r - 28),
        ("[7,246/35]", F(7), F(246, 35), lambda r: 14 - F(504, 287) * r),
        ("[246/35,43/6]", F(246, 35), F(43, 6), lambda r: 86 - 12 * r),
        ("[43/6,8]", F(43, 6), F(8), lambda r: F(1218, 301) * r - 29),
        ("[8,49/6]", F(8), F(49, 6), lambda r: F(119, 28) * r - 32),
    )
    answer = []
    for name, left, right, formula in pieces:
        left_margin = formula(left)
        right_margin = formula(right)
        require(left_margin >= 0 and right_margin >= 0, "negative affine envelope margin")
        answer.append((name, left_margin, right_margin))
    expected = (
        ("[25/4,41/6]", F(0), F(7, 3)),
        ("[41/6,7]", F(1), F(70, 41)),
        ("[7,246/35]", F(70, 41), F(58, 35)),
        ("[246/35,43/6]", F(58, 35), F(0)),
        ("[43/6,8]", F(0), F(145, 43)),
        ("[8,49/6]", F(2), F(65, 24)),
    )
    require(tuple(answer) == expected, "piecewise margin ledger moved")
    return tuple(answer)


def exhaustive_cone_audit() -> int:
    """All Q4-residual ordered shapes in the new cone through d=400."""
    rows = 0
    for d in range(4, EXHAUSTIVE_CONE_HEIGHT + 1):
        a_min = (6 * d + 48) // 49
        a_max = (4 * d) // 25
        for a in range(max(1, a_min), a_max + 1):
            if not (25 * a <= 4 * d and 6 * d <= 49 * a):
                continue
            # Put u=d-b and v=d-c.  The Q4-residual condition is u+v<=a.
            for v in range(1, a + 1):
                for u in range(v + 1, a - v + 1):
                    shape = a, d - u, d - v, d
                    require(shape[0] < shape[1] < shape[2] < shape[3], "deficit ordering")
                    require(2 * d <= sum(shape[:3]), "deficit residual")
                    require(exact_gate_numerator(shape) > 0, "exhaustive 25/4 cone counterexample")
                    rows += 1
    return rows


def nearest_integer(x: F) -> set[int]:
    floor = floor_fraction(x)
    return {floor - 2, floor - 1, floor, floor + 1, floor + 2, floor + 3}


def boundary_hunt() -> int:
    """Adversarial rational-cell audit far beyond the height-64 atlas."""
    r_boundaries = (F(25, 4), F(41, 6), F(7), F(246, 35), F(43, 6), F(8), F(49, 6))
    z_boundaries = (F(17, 3), F(35, 6), F(6), F(41, 6), F(7), F(43, 6), F(8))
    rows = 0
    for a in range(1, BOUNDARY_DENOMINATOR + 1):
        d_candidates: set[int] = set()
        for boundary in r_boundaries:
            d_candidates.update(nearest_integer(boundary * a))
        for d in d_candidates:
            if not (a < d and 25 * a <= 4 * d and 6 * d <= 49 * a):
                continue
            q_candidates = {d - a + 1, d - a + 2, d - 2, d - 1}
            for boundary in z_boundaries:
                q_candidates.update(nearest_integer(boundary * a))
            legal_q = sorted(q for q in q_candidates if a < q < d)
            for b, c in combinations(legal_q, 2):
                if 2 * d > a + b + c:
                    continue
                require(exact_gate_numerator((a, b, c, d)) > 0, "boundary counterexample")
                rows += 1
    return rows


def sharp_family_audit() -> tuple[F, F, int]:
    """Check exact threshold strictness and the family immediately below it."""
    rows = 0
    for n in range(3, 10_001):
        at_threshold = 4 * n, 24 * n, 24 * n + 1, 25 * n
        below = 4 * n, 24 * n, 24 * n + 1, 25 * n - 1
        require(2 * at_threshold[3] <= sum(at_threshold[:3]), "threshold residual")
        require(2 * below[3] <= sum(below[:3]), "below-threshold residual")
        threshold_margin = exact_gate_margin(at_threshold)
        below_margin = exact_gate_margin(below)
        require(threshold_margin == F(25, 4 * (24 * n + 1)), "threshold family formula")
        require(below_margin == -F(71 * n + 5, 4 * n * (24 * n + 1)), "sharp family formula")
        require(threshold_margin > 0 and below_margin < 0, "sharp-family sign")
        rows += 1
    return exact_gate_margin((12, 72, 73, 75)), exact_gate_margin((12, 72, 73, 74)), rows


def main() -> None:
    formula_rows = local_formula_audit()
    census = height64_census()
    champion_rows = champion_lattice_audit()
    margins = sufficient_margin_pieces()
    cone_rows = exhaustive_cone_audit()
    boundary_rows = boundary_hunt()
    threshold_margin, below_margin, sharp_rows = sharp_family_audit()

    print("THM-1213 exact fractional-discrepancy / 25/4 cone independent referee")
    print("arithmetic=integers and fractions.Fraction")
    print("epsilon(theta)=min(theta,1/7)-theta/7; theta_q={6q/(7a)}")
    print("nu_a(q)=floor(6q/(7a)+8/7)")
    print("gate: 24d/(7a)-7d*sum_q epsilon(theta_q)/q > 1+sum_q nu_a(q)")
    print(f"local epsilon/incidence rows={formula_rows}")
    print(
        "height-64 residual/exact telemetry="
        f"{census['residual']}/{census['exact']}; old_Q={census['old']}; "
        f"clean_25/4={census['clean']}; tightest={census['tightest']}"
    )
    print(f"champion rational-lattice rows (denominator<={CHAMPION_DENOMINATOR})={champion_rows}")
    print(f"piecewise sufficient margins (left,right)={margins}")
    print(f"exhaustive cone rows through d<={EXHAUSTIVE_CONE_HEIGHT}={cone_rows}")
    print(f"adversarial rational-boundary rows through a<={BOUNDARY_DENOMINATOR}={boundary_rows}")
    print(
        "sharp family (4n,24n,24n+1,25n[-1]) rows="
        f"{sharp_rows}; n=3 margins threshold/below={threshold_margin}/{below_margin}"
    )
    print("Tournament Analysis:")
    print("  vertices=phase obligations z; labels=(nu_a(z),epsilon_a(z)); gauge=r=d/a")
    print("  observable=kappa_r(z); switch=cost equality; tie path=sorted rational breakpoints")
    print("  champion path=6 -> 41/6 -> moving endpoint r -> 43/6 -> 8")
    print("  naked numerical-order tournament is transitive and destroys both metric labels")
    print("VERDICT: exact gate is valid; every Q4-residual shape with 4d>=25a closes")


if __name__ == "__main__":
    main()
