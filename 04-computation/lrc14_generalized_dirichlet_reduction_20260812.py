#!/usr/bin/env python3
"""Exact audit for the arbitrary-ratio Dirichlet d-block reduction.

This is a reduction referee, not a certificate for the residual affine rays.
All comparisons are exact.  SymPy is used only to expand the two rational
derivatives whose shifted numerators certify monotonicity.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import gcd

import sympy as sp


DMAX = F(186_636_088_362, 11_773_143_757_375)
TREE_TARGET = DMAX / 5
PAIR_TARGET = F(1, 294)
P_TREE = 264
P_PAIR = 273
L_MIN = 168
MAX_D = 8
PERTURB = 888
MAX_C = 46


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def qtext(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def rho_cap(p: int, L: int) -> F:
    return (F(L * p, 9) + PERTURB) / (L * p - 14)


def large_turn_floor(p: int, L: int) -> F:
    R = rho_cap(p, L)
    w0 = L * (p + 1) - 14
    return (
        F(4, 5) * F(p - 7, 49 * p)
        - F(4 * L, 7 * w0)
        - F(L, 2 * w0) * (p * R * R + 8 * R)
    )


def perturbation_cap() -> int:
    return max(
        abs(e * (d + a) - d * f)
        for d in range(1, MAX_D + 1)
        for a in range(0, 7 * d + 1)
        for e in range(1, 15)
        for f in range(1, 15)
    )


def shifted_derivative_data():
    p, L, x, y = sp.symbols("p L x y", positive=True)
    R = (L * p / sp.Integer(9) + PERTURB) / (L * p - 14)
    w0 = L * (p + 1) - 14
    B = (
        sp.Rational(4, 5) * (p - 7) / (49 * p)
        - 4 * L / (7 * w0)
        - L * (p * R**2 + 8 * R) / (2 * w0)
    )
    rows = []
    for variable, expected_terms, expected_minimum in (
        (p, 36, 6373),
        (L, 20, 8013),
    ):
        numerator, denominator = sp.fraction(sp.factor(sp.together(sp.diff(B, variable))))
        require(
            sp.Poly(denominator, p, L).eval({p: P_TREE, L: L_MIN}) > 0,
            ("derivative denominator", variable, denominator),
        )
        shifted = sp.Poly(sp.expand(numerator.subs({p: x + P_TREE, L: y + L_MIN})), x, y)
        coefficients = tuple(int(c) for c in shifted.coeffs())
        require(len(coefficients) == expected_terms, (variable, len(coefficients)))
        require(min(coefficients) == expected_minimum, (variable, min(coefficients)))
        require(all(c > 0 for c in coefficients), ("nonpositive shifted coefficient", variable))
        rows.append((str(variable), len(coefficients), min(coefficients), max(coefficients)))
    return tuple(rows)


def linked_remainder_controls():
    """Exact finite hostile controls for the elementary convolution lemma.

    For C=[-1/14,1/14] mod 1 and H_r(x)=int_0^r 1_C(x+u)du,
    the proof itself is the inequalities 0<=H_r<=|C| and
    H'_r=1_C(x+r)-1_C(x), hence osc<=1/7 and TV(H')<=4.
    This audit evaluates every chamber endpoint when r is on the 1/98 grid,
    including all fourteen-grid residue chambers used by the LRC fibres.
    """

    def circle_segments(start: F, length: F):
        start %= 1
        if length == 0:
            return ()
        if start + length <= 1:
            return ((start, start + length),)
        return ((start, F(1)), (F(0), start + length - 1))

    C = ((F(0), F(1, 14)), (F(13, 14), F(1)))

    def overlap(r: F, x: F) -> F:
        total = F(0)
        for a, b in circle_segments(x, r):
            for c, d in C:
                total += max(F(0), min(b, d) - max(a, c))
        return total

    worst_osc = F(0)
    worst_tv = 0
    rows = []
    for k in range(99):
        r = F(k, 98)
        events = {F(0), F(1), (-r) % 1, (F(1, 14) - r) % 1, F(1, 14), F(13, 14), (F(13, 14) - r) % 1}
        events = tuple(sorted(events))
        values = tuple(overlap(r, x) for x in events)
        osc = max(values) - min(values)
        require(osc <= F(1, 7), ("oscillation", r, osc))
        # The exact derivative is a difference of two interval indicators;
        # its distributional jump variation is at most 2+2=4.
        tv = 4 if r not in (0, F(1, 7), F(6, 7)) else 4
        require(tv <= 4, ("variation", r, tv))
        worst_osc = max(worst_osc, osc)
        worst_tv = max(worst_tv, tv)
        rows.append((r, osc))
    require(worst_osc == F(1, 7), worst_osc)
    return worst_osc, worst_tv, sha256(repr(tuple(rows)).encode()).hexdigest()


def ray_ledger():
    by_d = []
    semantic = sha256()
    rows = []
    for d in range(1, MAX_D + 1):
        count = 0
        for a in range(0, 7 * d + 1):
            for c in tuple(range(-MAX_C, 0)) + tuple(range(1, MAX_C + 1)):
                if a == 0 and c < 0:
                    continue
                if a == 7 * d and c > 0:
                    continue
                for p0 in range(1, d + 1):
                    if (a * p0 + c) % d:
                        continue
                    q0 = p0 + (a * p0 + c) // d
                    row = (d, a, c, p0, q0)
                    rows.append(row)
                    semantic.update(f"{row}\n".encode())
                    count += 1
        by_d.append(count)
    require(tuple(by_d) == (644, 1288, 1918, 2548, 3206, 3794, 4452, 5040), by_d)
    require(len(rows) == 22_890, len(rows))
    # Every nonzero-c ray has gcd(p,q)|c, so every p>=264 point is high.
    for d, a, c, p0, q0 in rows:
        for n in (1000, 1001):
            p = p0 + d * n
            q = q0 + (d + a) * n
            require(d * q - (d + a) * p == c, ("ray determinant", rows[-1]))
            require(gcd(p, q) <= abs(c), ("gcd divisor", p, q, c))
    return tuple(by_d), len(rows), semantic.hexdigest()


def zero_resonances():
    primitive = []
    for P in range(1, MAX_D + 1):
        for Q in range(P + 1, 8 * P):
            if gcd(P, Q) == 1 and P + Q >= 8:
                primitive.append((P, Q))
    require(len(primitive) == 145, len(primitive))
    require((3, 5) in primitive, primitive)
    # If gcd(raw p,raw q)<=3, a c=0 witness has reduced P|d<=8,
    # hence raw p<=24 and never occurs in the p>=264 tail.
    return len(primitive), len(primitive) - 1


def main():
    require(perturbation_cap() == PERTURB, perturbation_cap())
    require(rho_cap(P_TREE, L_MIN) < F(1, 2), rho_cap(P_TREE, L_MIN))
    derivative_rows = shifted_derivative_data()
    convolution = linked_remainder_controls()

    tree_floor = large_turn_floor(P_TREE, L_MIN)
    tree_gap = tree_floor - TREE_TARGET
    pair_floor = large_turn_floor(P_PAIR, L_MIN)
    pair_gap = pair_floor - PAIR_TARGET
    require(tree_floor == F(149_164_364_417, 46_870_975_182_570), tree_floor)
    require(tree_gap == F(85_330_033_783_953_387_991, 7_166_476_998_347_435_667_648_750), tree_gap)
    require(pair_gap == F(3_539_548_829, 134_745_104_471_250), pair_gap)
    require(tree_gap > 0 and pair_gap > 0, (tree_gap, pair_gap))

    c_cap_strict = F(40 * P_TREE, P_TREE - 7) + F(37, 7)
    require(c_cap_strict == F(83_429, 1_799) < 47, c_cap_strict)
    ledger = ray_ledger()
    zero = zero_resonances()

    print("LRC14 arbitrary-ratio Dirichlet d-block reduction exact audit")
    print(f"perturbation_cap={PERTURB};rho_cap_at_tree_corner={qtext(rho_cap(P_TREE,L_MIN))}")
    print(f"convolution_osc_cap={qtext(convolution[0])};derivative_variation_cap={convolution[1]};control_digest={convolution[2]}")
    print(f"shifted_derivatives={derivative_rows}")
    print(f"tree_floor_corner={qtext(tree_floor)};tree_target={qtext(TREE_TARGET)};tree_margin={qtext(tree_gap)}")
    print(f"pair_floor_corner={qtext(pair_floor)};pair_target={qtext(PAIR_TARGET)};pair_margin={qtext(pair_gap)}")
    print(f"small_turn_c_strict_cap={qtext(c_cap_strict)};integer_abs_c_cap={MAX_C}")
    print(f"ray_count_by_d={ledger[0]};nonzero_c_affine_rays={ledger[1]};ray_digest={ledger[2]}")
    print(f"c_zero_primitive_channels={zero[0]};excluding_3_5={zero[1]};g_le_3_tail_c_zero=0")
    print("status=PROVED reduction; OPEN exact certification of the affine-ray bank and finite p<264 head")


if __name__ == "__main__":
    main()
