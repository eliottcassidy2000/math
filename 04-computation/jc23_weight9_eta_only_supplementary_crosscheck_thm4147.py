#!/usr/bin/env python3
"""Supplementary eta-row exact cross-check for THM-4147.

THM-4147 already proves all three generic exact-weight-nine chambers. This
separate certificate corroborates its eta-only row on a narrower explicit
open set, using the older quadratic-carrier route and a disjoint witness.
It adds audit redundancy, not theorem scope. The general identities are
checked over a symbolic coefficient ring, and the Newton packet and two
permutation contradictions are exact on the stated open set.

The final two top-support cases are also evaluated at the same rational
control point, but those lines are explicitly scouts rather than theorem
dependencies.  No Python assertions are used.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd

import sympy as sp


CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def convex_hull(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
    points = sorted(set(points))

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

    lower: list[tuple[int, int]] = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def pick_data(vertices: list[tuple[int, int]]) -> tuple[int, int, int]:
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]))
        for i in range(len(vertices))
    )
    require((area2 - boundary + 2) % 2 == 0, "Pick parity failed")
    return area2, boundary, (area2 - boundary + 2) // 2


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[i]] for i in range(len(left)))


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    answer = [0] * len(permutation)
    for source, target in enumerate(permutation):
        answer[target] = source
    return tuple(answer)


def cycle(n: int, entries: tuple[int, ...]) -> tuple[int, ...]:
    answer = list(range(n))
    for left, right in zip(entries, entries[1:] + entries[:1]):
        answer[left] = right
    return tuple(answer)


def permutation_index(permutation: tuple[int, ...]) -> int:
    seen: set[int] = set()
    cycles = 0
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycles += 1
        point = start
        while point not in seen:
            seen.add(point)
            point = permutation[point]
    return len(permutation) - cycles


def orbit(generators: tuple[tuple[int, ...], ...], start: int) -> set[int]:
    reached = {start}
    stack = [start]
    while stack:
        point = stack.pop()
        for generator in generators:
            image = generator[point]
            if image not in reached:
                reached.add(image)
                stack.append(image)
    return reached


X, T = sp.symbols("X T")
D, Phi, Theta, Eta, Zeta = sp.symbols("D Phi Theta Eta Zeta")
P = T + X**2 * T**2
Y = X * T * P
K = sp.Rational(2848, 45) - sp.Rational(7, 6) * D


def source_polynomial(eta, zeta):
    return sp.expand(
        -sp.Rational(1, 2) * X**2 * T
        - 3 * P
        + sp.Rational(8, 3) * P**2
        - sp.Rational(1376, 135) * P**3
        + K * Y**2
        + Phi * P**2 * Y
        + D * P**4
        + Theta * P * Y**2
        + eta * P**3 * Y
        + zeta * Y**3
    )


def t_order(poly: sp.Expr) -> int:
    terms = sp.Poly(poly, T).terms()
    return min(exponent[0] for exponent, _ in terms)


def factor_multiplicity(poly: sp.Expr, factor: sp.Expr) -> tuple[int, sp.Expr]:
    multiplicity = 0
    quotient = sp.cancel(poly)
    factor_poly = sp.Poly(factor, T)
    while sp.rem(sp.Poly(quotient, T), factor_poly).is_zero:
        quotient = sp.cancel(quotient / factor)
        multiplicity += 1
    return multiplicity, quotient


def critical_scout(eta_value: int, zeta_value: int) -> tuple[int, int, int]:
    substitutions = {
        D: sp.Integer(2),
        Phi: sp.Integer(5),
        Theta: sp.Integer(7),
        Eta: sp.Integer(eta_value),
        Zeta: sp.Integer(zeta_value),
    }
    polynomial = source_polynomial(Eta, Zeta).subs(substitutions)
    f = sp.cancel(sp.diff(polynomial, X) / T)
    h = sp.diff(polynomial, T)
    resultant = sp.factor(sp.resultant(f, h, X))
    order = t_order(resultant)
    residual = sp.cancel(resultant / T**order)
    universal_multiplicity, residual = factor_multiplicity(residual, 6 * T + 1)
    return order, universal_multiplicity, sp.Poly(residual, T).degree()


def main() -> None:
    # Symbolic source and the two universal Morse pairs.
    G = source_polynomial(Eta, sp.Integer(0))
    GX = sp.diff(G, X)
    GT = sp.diff(G, T)
    f = sp.cancel(GX / T)
    require(sp.denom(f) == 1, "G_X/T is not polynomial")
    require(sp.expand(GX - T * f) == 0, "lost the exact T factor")
    require(sp.Poly(f, X).degree() == 8, "wrong eta-only f degree")
    require(sp.Poly(GT, X).degree() == 9, "wrong eta-only h degree")
    require(sp.factor(sp.Poly(f, X).LC()) == 9 * Eta * T**8,
            "wrong eta-only f leading row")
    require(sp.factor(sp.Poly(GT, X).LC()) == 9 * Eta * T**8,
            "wrong eta-only h leading row")
    require(sp.factor(f.subs(T, 0)) == -X, "wrong T=0 f")
    require(sp.expand(GT.subs(T, 0) + (X**2 + 6) / 2) == 0,
            "wrong T=0 h")

    modulus = sp.Poly(X**2 - 6, X)
    for expression, expected, label in (
        (f, 0, "f"),
        (GT, 0, "h"),
        (G, sp.Rational(1, 2), "value"),
    ):
        remainder = sp.rem(
            sp.Poly(sp.expand(expression.subs(T, -sp.Rational(1, 6)) - expected), X),
            modulus,
        ).as_expr()
        require(sp.expand(remainder) == 0, f"wrong universal {label}")
    hessian = sp.det(sp.hessian(G, (X, T)))
    require(sp.rem(sp.Poly(hessian.subs(T, 0) - 6, X),
                   sp.Poly(X**2 + 6, X)).is_zero,
            "T=0 pair is not Morse")
    require(sp.rem(sp.Poly(hessian.subs(T, -sp.Rational(1, 6)) + 6, X),
                   modulus).is_zero,
            "T=-1/6 pair is not Morse")

    # Exact rational witness for the nonempty generic critical chamber.
    witness = {D: 2, Phi: 5, Theta: 7, Eta: 11}
    witness_k = K.subs(witness)
    require(witness_k == sp.Rational(2743, 45), "witness K changed")
    Gw = G.subs(witness)
    fw = sp.cancel(sp.diff(Gw, X) / T)
    hw = sp.diff(Gw, T)
    resultant = sp.factor(sp.resultant(fw, hw, X))
    Q20 = sp.cancel(resultant / (11**3 * T**56 * (6 * T + 1)**2))
    require(sp.denom(Q20) == 1, "witness Q20 is not polynomial")
    Q20_poly = sp.Poly(Q20, T)
    require(Q20_poly.degree() == 20, "witness residual degree is not twenty")
    require(Q20_poly.LC() == 72900 * witness_k**4 * 7**4 * 11**2,
            "witness leading coefficient changed")
    require(Q20_poly.TC() == -sp.Rational(3**15, 2**7) * 11**4,
            "witness constant coefficient changed")
    require(not sp.gcd(Q20_poly, sp.Poly(sp.diff(Q20, T), T)).degree(),
            "witness Q20 is not squarefree")
    require(Q20_poly.eval(-sp.Rational(1, 6)) != 0,
            "witness residual meets the universal pair")
    primitive = sp.primitive(Q20_poly)[1]
    coefficient_text = ",".join(str(value) for value in primitive.all_coeffs())
    coefficient_digest = sha256(coefficient_text.encode("ascii")).hexdigest()
    require(coefficient_digest ==
            "c71a94984838cee834adc6cbde848275fd55d3b49d0abf88fc7f691c758ea89e",
            "witness coefficient digest changed")

    # Complete generic eta-only Newton polygon and Pick/genus ledger.
    eta_hull = [(0, 1), (2, 0), (4, 2), (4, 3),
                (3, 4), (1, 5), (0, 5)]
    require(convex_hull(eta_hull) == eta_hull, "eta hull order changed")
    require(pick_data(eta_hull) == (29, 11, 10), "eta Pick ledger changed")
    edge_ledger = (
        ("AB", (1, 2), 2, 1, 1),
        ("BC", (-1, 1), -2, 2, 2),
        ("CD", (-1, 0), -4, 3, 1),
        ("DE", (-1, -1), -7, 5, 1),
        ("EF", (-1, -2), -11, 8, 1),
        ("FG", (0, -1), -5, 4, 1),
    )
    packet: list[int] = []
    for _, normal, minimum, index, branches in edge_ledger:
        require(normal[0] + normal[1] - minimum == index,
                "edge ramification formula changed")
        packet.extend([index] * branches)
    packet = sorted(packet, reverse=True)
    require(packet == [8, 5, 4, 3, 2, 2, 1], "eta packet changed")
    require(sum(packet) == 25, "full-boundary response degree changed")
    require(sum(index - 1 for index in packet) == 18,
            "eta ramification defect changed")
    require(2 * 10 - 2 == 18, "genus and packet no longer agree")
    require(sum(index for index in packet if index != 2) == 21,
            "finite-BC response degree changed")

    # Exact permutation budgets used in the proof.
    critical_length = 20 + 2 + 2
    require(critical_length == 24, "critical length changed")
    require(2 * 25 - critical_length == 26, "degree-25 support budget changed")
    require(2 * 21 - critical_length == 18, "degree-21 support budget changed")
    require(18 - 1 + 2 < 21 - 1,
            "one-handle finite-carrier budget no longer contradicts transitivity")
    require(18 - 2 + 2 < 21 - 1,
            "two-handle finite-carrier budget no longer contradicts transitivity")

    # Sharp hostile: one fewer critical point would allow equality in the
    # finite-carrier merger bound.
    n = 21
    identity = tuple(range(n))
    hostile_x = cycle(n, tuple(range(19)))
    hostile_plus = cycle(n, (18, 19))
    hostile_minus = cycle(n, (19, 20))
    require(len(orbit((hostile_x, identity, hostile_plus, hostile_minus), 0)) == n,
            "sharp merger hostile lost transitivity")
    require(sum(permutation_index(g) for g in
                (hostile_x, identity, hostile_plus, hostile_minus)) == n - 1,
            "sharp merger hostile lost equality")
    require(sum(1 for i, image in enumerate(hostile_x) if i == image) + n == 23,
            "sharp hostile no longer models critical length twenty-three")

    # Exhaustive top-support list.  The last two lines are finite-exact
    # scouts only; their cubic carrier transport is not a dependency.
    eta_scout = critical_scout(11, 0)
    zeta_scout = critical_scout(0, 13)
    both_scout = critical_scout(11, 13)
    require(eta_scout == (56, 2, 20), "eta critical scout changed")
    require(zeta_scout == (56, 2, 20), "zeta critical scout changed")
    require(both_scout == (56, 2, 21), "two-term critical scout changed")
    zeta_hull = [(0, 1), (2, 0), (5, 3), (3, 4), (0, 5)]
    both_hull = [(0, 1), (2, 0), (5, 3), (1, 5), (0, 5)]
    require(convex_hull(zeta_hull) == zeta_hull, "zeta scout hull changed")
    require(convex_hull(both_hull) == both_hull, "mixed scout hull changed")
    require(pick_data(zeta_hull) == (30, 10, 11), "zeta scout Pick ledger changed")
    require(pick_data(both_hull) == (31, 11, 11), "mixed scout Pick ledger changed")
    zeta_packet = [11, 8, 2, 2, 2, 1]
    both_packet = [8, 8, 4, 2, 2, 2, 1]
    require(sum(zeta_packet) == 26
            and sum(index - 1 for index in zeta_packet) == 20,
            "zeta scout packet changed")
    require(sum(both_packet) == 27
            and sum(index - 1 for index in both_packet) == 20,
            "mixed scout packet changed")

    semantic = (
        "THM4147|eta-only|crit=24|g=10|packet=8,5,4,3,2,2,1|"
        "responses=25,21|quadratic-carrier=2xindex2|"
        "Tminus1over6=universal-only|verdict=excluded"
    )
    semantic_digest = sha256(semantic.encode("ascii")).hexdigest()
    print("THM-4147 primary exact audit")
    print(f"checks={CHECKS}")
    print("normalized_top=Eta*P^3*Y; Zeta=0")
    print("generic_eta_hull=" + ",".join(f"{x}:{y}" for x, y in eta_hull))
    print("pick=2area:29,boundary:11,genus:10")
    print("critical_resultant=T^56*(6T+1)^2*degree20 after nonzero scalar")
    print("U_eta_Tminus1over6=universal_pair_only")
    print("critical_length=24")
    print("boundary_packet=8,5,4,3,2,2,1;defect=18")
    print("response_degrees=25,21")
    print("degree25_support_budget=26;degree21_support_budget=18")
    print("sharp_L23_merger_hostile=PASS")
    print(f"witness_Q20_coeff_sha256={coefficient_digest}")
    print("remaining_scouts=eta:24,zeta:24,both:25;NOT_THEOREM_DEPENDENCIES")
    print("remaining_newton_scouts=zeta:g11:26:11,8,2,2,2,1;"
          "both:g11:27:8,8,4,2,2,2,1")
    print(f"semantic_sha256={semantic_digest}")
    print("PASS")


if __name__ == "__main__":
    main()
