#!/usr/bin/env python3
"""
Exact symbolic closure for the GMC(2) span-6 {+3,+1,-1,-3} stratum.

This strengthens death-star-S73's numerical B2 residual check.  For

    P = a Z^3 + b Z + c W + d W^3       (W = Zbar)

the Wick moments are computed exactly.  In the all-nonzero B2 case,
E[P^2]=0 gives bc = -6ad.  Writing u = ad, x = a c^3, y = b^3 d, the
fourth moment reduces to

    E[P^4] = 24*(x + y + 54 u^2).

The sixth moment reduces to

    E[P^6] = 38880*u*(x+y) + 2566080*u^3,

so on E[P^4]=0 it becomes 466560*u^3, nonzero when a,d are nonzero.

The file is deliberately dependency-free: it uses only exact integer arithmetic.
It also records the required Tournament Analysis/assumption challenge for the
repo's research workflow.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import product
from math import factorial


NAMES = ("a", "b", "c", "d")


def moment_terms(m: int) -> dict[tuple[int, int, int, int], int]:
    """Return the exact Wick moment polynomial E[P^m]."""
    # (Z-power, W-power) contributions for a,b,c,d.
    terms = ((3, 0), (1, 0), (0, 1), (0, 3))
    out: dict[tuple[int, int, int, int], int] = defaultdict(int)
    for combo in product(range(m + 1), repeat=4):
        if sum(combo) != m:
            continue
        zpow = sum(combo[i] * terms[i][0] for i in range(4))
        wpow = sum(combo[i] * terms[i][1] for i in range(4))
        if zpow != wpow:
            continue
        coeff = factorial(m) * factorial(zpow)
        for e in combo:
            coeff //= factorial(e)
        out[combo] += coeff
    return dict(out)


def monomial_string(exp: tuple[int, int, int, int]) -> str:
    pieces = []
    for name, power in zip(NAMES, exp):
        if power == 1:
            pieces.append(name)
        elif power > 1:
            pieces.append(f"{name}^{power}")
    return "".join(pieces) or "1"


def polynomial_string(poly: dict[tuple[int, int, int, int], int]) -> str:
    if not poly:
        return "0"
    parts = []
    for exp, coeff in sorted(poly.items()):
        mon = monomial_string(exp)
        if mon == "1":
            parts.append(str(coeff))
        else:
            parts.append(f"{coeff}{mon}")
    return " + ".join(parts)


def reduce_b2_monomial(exp: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    """Reduce a charge-zero monomial using bc=-6ad.

    The output is (factor, u_power, x_power, y_power), meaning

        factor * u^u_power * x^x_power * y^y_power,

    where u=ad, x=a c^3, y=b^3 d.
    """
    ia, ib, ic, id_ = exp
    pairs = min(ib, ic)
    factor = (-6) ** pairs
    ia += pairs
    id_ += pairs
    ib -= pairs
    ic -= pairs

    if ib == 0 and ic == 0:
        if ia != id_:
            raise ValueError(f"not a pure u power after reduction: {exp}")
        return factor, ia, 0, 0

    if ib == 0:
        if ic % 3 != 0:
            raise ValueError(f"remaining c power not divisible by 3: {exp}")
        xp = ic // 3
        if ia < xp or ia - xp != id_:
            raise ValueError(f"cannot rewrite as x^r u^k: {exp}")
        return factor, id_, xp, 0

    if ic == 0:
        if ib % 3 != 0:
            raise ValueError(f"remaining b power not divisible by 3: {exp}")
        yp = ib // 3
        if id_ < yp or id_ - yp != ia:
            raise ValueError(f"cannot rewrite as y^r u^k: {exp}")
        return factor, ia, 0, yp

    raise ValueError(f"unreduced b and c powers remain: {exp}")


def reduce_b2(poly: dict[tuple[int, int, int, int], int]) -> dict[tuple[int, int, int], int]:
    """Reduce a polynomial to the B2 basis u^i x^j y^k."""
    out: dict[tuple[int, int, int], int] = defaultdict(int)
    for exp, coeff in poly.items():
        factor, up, xp, yp = reduce_b2_monomial(exp)
        out[(up, xp, yp)] += coeff * factor
    return {k: v for k, v in out.items() if v}


def basis_string(basis: tuple[int, int, int]) -> str:
    names = ("u", "x", "y")
    pieces = []
    for name, power in zip(names, basis):
        if power == 1:
            pieces.append(name)
        elif power > 1:
            pieces.append(f"{name}^{power}")
    return "*".join(pieces) or "1"


def reduced_string(poly: dict[tuple[int, int, int], int]) -> str:
    if not poly:
        return "0"
    parts = []
    for basis, coeff in sorted(poly.items()):
        mon = basis_string(basis)
        if mon == "1":
            parts.append(str(coeff))
        else:
            parts.append(f"{coeff}*{mon}")
    return " + ".join(parts)


def substitute_x_plus_y_in_e6(red_e6: dict[tuple[int, int, int], int]) -> int:
    """Use x+y=-54u^2 on the specific reduced E6 form."""
    ux = red_e6.get((1, 1, 0), 0)
    uy = red_e6.get((1, 0, 1), 0)
    pure = red_e6.get((3, 0, 0), 0)
    assert ux == uy, red_e6
    assert set(red_e6) <= {(1, 1, 0), (1, 0, 1), (3, 0, 0)}, red_e6
    return pure - 54 * ux


def tournament_analysis() -> str:
    vertices = [
        ("symbolic E6 residual", 5),
        ("resultant E4 rung", 4),
        ("case split A/B1", 4),
        ("numerical B2 probe", 2),
        ("random real nullcone search", 1),
    ]
    ordered = " > ".join(name for name, _ in sorted(vertices, key=lambda t: -t[1]))
    hist = {}
    for _, score in vertices:
        hist[score] = hist.get(score, 0) + 1
    return "\n".join(
        [
            "Tournament Analysis:",
            "  vertices: proof obligations/certificates, not runners or charges",
            "  pairwise observable: exactness + assumptions removed + formalization readiness",
            "  switch/gauge: certificate A beats B if it removes a weaker/non-exact premise",
            f"  score histogram: {dict(sorted(hist.items()))}",
            f"  Hamiltonian path: {ordered}",
            "  directed cycles: none in this gauge; the ranking is transitive",
            "  tie path: if exactness ties, prefer the lower moment degree, then fewer side conditions",
            "  assumption challenge: charges, moment degrees, residual variables, case branches,",
            "    and proof obligations were considered as vertices.  The chosen quotient preserves",
            "    the GMC2 predicate 'two-sided nullcone is impossible' and intentionally destroys",
            "    coefficient phase information after replacing it by nonzero monomial certificates.",
        ]
    )


def main() -> None:
    e2 = moment_terms(2)
    e4 = moment_terms(4)
    e6 = moment_terms(6)

    print("Exact span-6 GMC(2) symbolic residual certificate")
    print("P = a Z^3 + b Z + c Zbar + d Zbar^3")
    print()
    print(f"E[P^2] = {polynomial_string(e2)}")
    print(f"E[P^4] = {polynomial_string(e4)}")
    print(f"E[P^6] = {polynomial_string(e6)}")
    print()

    red_e4 = reduce_b2(e4)
    red_e6 = reduce_b2(e6)
    residual = substitute_x_plus_y_in_e6(red_e6)

    print("B2 all-nonzero branch with u=ad, x=ac^3, y=b^3d and bc=-6ad:")
    print(f"  E[P^4] -> {reduced_string(red_e4)}")
    print("            = 24*(x + y + 54*u^2)")
    print(f"  E[P^6] -> {reduced_string(red_e6)}")
    print(f"  E[P^6] under E[P^4]=0 -> {residual}*u^3")
    print()

    assert red_e4 == {(0, 1, 0): 24, (0, 0, 1): 24, (2, 0, 0): 1296}
    assert red_e6 == {(1, 1, 0): 38880, (1, 0, 1): 38880, (3, 0, 0): 2566080}
    assert residual == 466560

    print("Case certificates:")
    print("  a=0: E2=2bc, so two-sided b!=0 forces c=0; then E4=24*b^3*d, forcing d=0.")
    print("  d=0,a!=0: E2=2bc, so two-sided c!=0 forces b=0; then E4=24*a*c^3, forcing c=0.")
    print("  a*d!=0: E2 and E4 force the B2 relations above; E6=466560*(ad)^3 != 0.")
    print()
    print("Conclusion: no two-sided nullcone member exists in the constant-coefficient")
    print("{+3,+1,-1,-3} span-6 stratum; moments 2,4,6 suffice.")
    print()
    print(tournament_analysis())


if __name__ == "__main__":
    main()
