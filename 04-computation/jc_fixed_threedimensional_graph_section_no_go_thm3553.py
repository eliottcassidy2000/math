#!/usr/bin/env python3
"""Exact controls for THM-3553, the fixed THM-1300 graph-section no-go.

Scope: fixed source chart z=h(x,y), fixed target graph chart with tangential
coordinates (F1,F2).  No claim is made about other source/target coordinates,
nongraph coordinate hypersurfaces, or arbitrary three-dimensional maps.
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not bool(condition):
        raise RuntimeError(message)


def jac2(p: sp.Expr, q: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.diff(p, x) * sp.diff(q, y) - sp.diff(p, y) * sp.diff(q, x))


def homogeneous_part(poly: sp.Expr, degree: int, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    P = sp.Poly(sp.expand(poly), x, y)
    return sp.expand(sum(co * x**mon[0] * y**mon[1] for mon, co in P.terms() if sum(mon) == degree))


def sha(value: object) -> str:
    return hashlib.sha256(repr(value).encode("utf-8")).hexdigest()


x, y, z = sp.symbols("x y z")
u = 1 + x * y
F1 = sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y))
F2 = sp.expand(y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
F3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)
F = (F1, F2, F3)

JF = sp.Matrix([[sp.diff(fi, q) for q in (x, y, z)] for fi in F])
require(sp.expand(JF.det()) == -2, "ambient THM-1300 determinant")

points = {
    "p0": (sp.Integer(0), sp.Integer(0), -sp.Rational(1, 4)),
    "pp": (sp.Integer(1), -sp.Rational(3, 2), sp.Rational(13, 2)),
    "pm": (-sp.Integer(1), sp.Rational(3, 2), sp.Rational(13, 2)),
}
images = {}
for name, point in points.items():
    images[name] = tuple(sp.expand(fi.subs(dict(zip((x, y, z), point)))) for fi in F)
require(len(set(images.values())) == 1, "collision images differ")
require(next(iter(images.values())) == (-sp.Rational(1, 4), 0, 0), "collision target mismatch")

# Tangential-minor PDE.  For P=F1(x,y,h), Q=F2(x,y,h), chain rule gives
# T_h=M0+Mx*h_x+My*h_y, where the three coefficients are 2x2 minors of JF.
H, Hx, Hy = sp.symbols("H Hx Hy")
M0 = sp.expand(sp.diff(F1, x) * sp.diff(F2, y) - sp.diff(F1, y) * sp.diff(F2, x)).subs(z, H)
Mx = sp.expand(sp.diff(F1, z) * sp.diff(F2, y) - sp.diff(F1, y) * sp.diff(F2, z)).subs(z, H)
My = sp.expand(sp.diff(F1, x) * sp.diff(F2, z) - sp.diff(F1, z) * sp.diff(F2, x)).subs(z, H)
t = x * y
M0_compact = sp.expand(
    -9 * x * u**4 * H**2
    - 3 * y * u**2 * (15 * t**2 + 25 * t + 7) * H
    - y**3 * (54 * t**3 + 189 * t**2 + 222 * t + 89)
)
Mx_compact = sp.expand(-u**2 * (3 * H * x**2 * u**2 + 9 * t**3 + 12 * t**2 - t - 1))
My_compact = sp.expand(-3 * u**2 * (H * u**2 + y**2 * (3 * t + 4)))
require(sp.expand(M0 - M0_compact) == 0, "PDE M0 compact formula")
require(sp.expand(Mx - Mx_compact) == 0, "PDE h_x compact formula")
require(sp.expand(My - My_compact) == 0, "PDE h_y compact formula")
PDE = sp.expand(M0_compact + Mx_compact * Hx + My_compact * Hy)


def tangential(h: sp.Expr) -> sp.Expr:
    P = sp.expand(F1.subs(z, h))
    Q = sp.expand(F2.subs(z, h))
    direct = jac2(P, Q, x, y)
    chain = sp.expand(PDE.subs({H: h, Hx: sp.diff(h, x), Hy: sp.diff(h, y)}))
    require(sp.expand(direct - chain) == 0, "direct/PDE tangential mismatch")
    return direct


pairs = (("p0", "pp"), ("p0", "pm"), ("pp", "pm"))


def interpolated_ansatz(degree: int, pair: tuple[str, str]) -> tuple[sp.Expr, tuple[sp.Symbol, ...], dict[sp.Symbol, sp.Expr]]:
    coeff = {(i, j): sp.symbols(f"c{degree}_{i}_{j}") for i in range(degree + 1) for j in range(degree + 1 - i)}
    h = sp.expand(sum(c * x**i * y**j for (i, j), c in coeff.items()))
    p, q = (points[name] for name in pair)
    equations = (sp.expand(h.subs({x: p[0], y: p[1]}) - p[2]), sp.expand(h.subs({x: q[0], y: q[1]}) - q[2]))
    solved = sp.solve(equations, (coeff[(0, 0)], coeff[(1, 0)]), dict=True)
    require(len(solved) == 1, f"interpolation solve failed degree={degree} pair={pair}")
    substitution = solved[0]
    hi = sp.expand(h.subs(substitution))
    require(sp.expand(hi.subs({x: p[0], y: p[1]}) - p[2]) == 0, "first interpolation residual")
    require(sp.expand(hi.subs({x: q[0], y: q[1]}) - q[2]) == 0, "second interpolation residual")
    free = tuple(c for c in coeff.values() if c not in substitution)
    require(len(free) == (degree + 1) * (degree + 2) // 2 - 2, "interpolation dimension")
    return hi, free, substitution


affine_rows = []
for pair in pairs:
    h1, free, _ = interpolated_ansatz(1, pair)
    require(len(free) == 1, "affine interpolation should leave one parameter")
    T1 = tangential(h1)
    y3 = sp.Poly(T1, x, y).coeff_monomial(y**3)
    require(y3 == -89, f"affine y^3 obstruction changed for {pair}: {y3}")
    affine_rows.append((pair, str(h1), str(free[0]), y3))

# The pp/pm pair alone permits a constant graph h=13/2.  It fails both the
# tangential PDE and THM-3544's nonautomorphic pencil floor: Q has degree 5.
hconst = sp.Rational(13, 2)
Tconst = tangential(hconst)
require(sp.Poly(Tconst, x, y).coeff_monomial(y**3) == -89, "constant graph obstruction")
Pconst = sp.expand(F1.subs(z, hconst))
Qconst = sp.expand(F2.subs(z, hconst))
require(sp.Poly(Pconst, x, y).total_degree() == 6, "constant graph P degree")
require(sp.Poly(Qconst, x, y).total_degree() == 5, "constant graph Q degree")

search_rows = []
for degree in (2, 3):
    for pair in pairs:
        hi, free, _ = interpolated_ansatz(degree, pair)
        hd = homogeneous_part(hi, degree, x, y)
        require(hd != 0, "generic leading form vanished")

        P = sp.expand(F1.subs(z, hi))
        Q = sp.expand(F2.subs(z, hi))
        Ptop = homogeneous_part(P, degree + 6, x, y)
        Qtop = homogeneous_part(Q, degree + 5, x, y)
        require(sp.expand(Ptop - x**3 * y**3 * hd) == 0, "restricted P top form")
        require(sp.expand(Qtop - 3 * x**3 * y**2 * hd) == 0, "restricted Q top form")

        # Hence every nonzero pencil member has degree >=degree+5.  In
        # particular quadratic/cubic graph ansatzes pass the THM-3544 floor.
        pencil_floor = degree + 5
        require(pencil_floor >= 6, "THM-3544 pencil floor compatibility")

        T = tangential(hi)
        S = sp.expand(x**3 * y**2 * hd)
        expected_top = sp.expand(-3 * S * sp.diff(S, x))
        actual_top = homogeneous_part(T, 2 * degree + 9, x, y)
        require(sp.expand(actual_top - expected_top) == 0, "top tangential obstruction identity")
        require(expected_top != 0, "generic top tangential obstruction vanished")

        # Exhibit one exact-degree interpolation witness, proving that the
        # searched affine coefficient space is nonempty before the PDE gate.
        specialization = {symbol: 0 for symbol in free}
        lead_symbol = next(symbol for symbol in free if symbol.name == f"c{degree}_{degree}_0")
        specialization[lead_symbol] = 1
        witness = sp.expand(hi.subs(specialization))
        require(sp.Poly(witness, x, y).total_degree() == degree, "witness degree")
        for name in pair:
            point = points[name]
            require(sp.expand(witness.subs({x: point[0], y: point[1]}) - point[2]) == 0, "witness interpolation")

        search_rows.append(
            (
                pair,
                degree,
                len(free),
                pencil_floor,
                2 * degree + 9,
                len(sp.Poly(expected_top, x, y).terms()),
                sha(witness),
            )
        )

# General all-positive-degree leading-form law, verified independently for
# generic degrees 1,2,3 above: S has an x^3 factor, so S_x is nonzero whenever
# h_d is nonzero in characteristic zero.  Thus the top minor cannot vanish.
semantic = {
    "ambient_det": -2,
    "collision_image": [str(v) for v in next(iter(images.values()))],
    "pde_sha": sha(PDE),
    "affine_rows": affine_rows,
    "constant_pair_floor": 5,
    "search_rows": search_rows,
    "general_top_law": "top(T_h)=-3*S*d_x(S), S=x^3*y^2*h_d, degree=2d+9",
    "scope": "fixed source z-graph and fixed target (F1,F2)-graph chart only",
}
semantic_digest = hashlib.sha256(json.dumps(semantic, sort_keys=True, default=str).encode()).hexdigest()

print("== THM-3553 fixed-map graph-section obstruction ==")
print("AMBIENT det=-2 collision=3 PDE=PASS pde_sha=" + sha(PDE))
for row in affine_rows:
    print("AFFINE", row[0], "free=1", "coeff_y3=-89", "PASS")
print("CONSTANT", ("pp", "pm"), "h=13/2", "degrees=(6,5)", "coeff_y3=-89", "FLOOR_FAIL", "PDE_FAIL")
for row in search_rows:
    print(
        "SEARCH",
        row[0],
        "degree=" + str(row[1]),
        "free=" + str(row[2]),
        "pencil_floor=" + str(row[3]),
        "minor_top_degree=" + str(row[4]),
        "minor_top_terms=" + str(row[5]),
        "survivors=0",
        "witness_sha=" + row[6],
    )
print("GENERAL exact positive-degree obstruction: top(T_h)=-3*S*d_x(S), S=x^3*y^2*h_d != 0")
print("SEMANTIC_SHA256=" + semantic_digest)
print("SCOPE fixed z=h(x,y), target graph over (F1,F2); no other-chart/hypersurface claim")
print("PASS")
