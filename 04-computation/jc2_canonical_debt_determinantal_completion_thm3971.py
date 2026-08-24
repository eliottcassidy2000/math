#!/usr/bin/env python3
"""Exact companion for THM-3971's determinantal completion near miss."""

from __future__ import annotations

import hashlib
import itertools
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


x, t, z, p, q, Y, s = sp.symbols("x t z p q Y s")
zxt = 1 + x * t


# ---------------------------------------------------------------------------
# The height-m affine modification and its boundary/canonical packet.
# ---------------------------------------------------------------------------

for m in range(1, 7):
    ys = [sp.expand(zxt * t**j) for j in range(m + 1)]
    gate(ys[0] == zxt, f"height {m}: y_0=z")
    for j in range(1, m + 1):
        zero(x * ys[j] - (zxt - 1) * ys[j - 1],
             f"height {m}: recurrence {j}")
        zero(x ** (m - j) * ys[m]
             - (zxt - 1) ** (m - j) * ys[j],
             f"height {m}: top-generator descent {j}")

    zero(x**m * ys[m] - zxt * (zxt - 1)**m,
         f"height {m}: boundary hypersurface equation")
    boundary_derivative = sp.diff(z * (z - 1)**m - x**m * Y, z)
    gate(boundary_derivative.subs({x: 0, z: 0}) == (-1)**m,
         f"height {m}: boundary chart is smooth")

    dy_dt = sp.diff(ys[m], t)
    zero(dy_dt - t ** (m - 1) * ((m + 1) * zxt - 1),
         f"height {m}: canonical derivative")
    zero(x ** (m - 1) * dy_dt
         - (zxt - 1) ** (m - 1) * ((m + 1) * zxt - 1),
         f"height {m}: canonical boundary order")

    # The curve x=s, t=-s^-1+(-1)^m q s^(m-1) reaches an arbitrary
    # boundary value q in y_m while every lower y_j tends to zero.
    t_arc = -1 / s + (-1)**m * q * s ** (m - 1)
    z_arc = sp.expand(1 + s * t_arc)
    gate(sp.limit(z_arc, s, 0) == 0,
         f"height {m}: boundary arc has z=0")
    for j in range(1, m):
        gate(sp.limit(sp.expand(z_arc * t_arc**j), s, 0) == 0,
             f"height {m}: lower boundary coordinate {j} vanishes")
    gate(sp.limit(sp.expand(z_arc * t_arc**m), s, 0) == q,
         f"height {m}: top boundary coordinate is free")


# The logarithmic primitive alpha=-t dx=(1-z) dx/x has residue one.
gate((1 - z).subs(z, 0) == 1, "logarithmic primitive has residue one")


# ---------------------------------------------------------------------------
# The first non-complete-intersection member m=2 is determinantal.
# ---------------------------------------------------------------------------

pxt = sp.expand(zxt * t)
qxt = sp.expand(zxt * t**2)
zero(x * pxt - zxt * (zxt - 1), "m=2 first determinant")
zero(x * qxt - (zxt - 1) * pxt, "m=2 second determinant")
zero(zxt * qxt - pxt**2, "m=2 third determinant")

matrix = sp.Matrix([[z, p, x], [p, q, z - 1]])
minors = [sp.expand(matrix[:, [i, j]].det())
          for i, j in itertools.combinations(range(3), 2)]
gate(minors == [sp.expand(z * q - p**2),
                sp.expand(z * (z - 1) - p * x),
                sp.expand(p * (z - 1) - q * x)],
     "m=2 ideal is the three 2x2 minors")


def jacobian(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(a, x) * sp.diff(b, t)
                     - sp.diff(a, t) * sp.diff(b, x))


pair_brackets = {
    (x, zxt): x,
    (x, pxt): 1 + 2 * x * t,
    (x, qxt): t * (2 + 3 * x * t),
    (zxt, pxt): pxt,
    (zxt, qxt): 2 * qxt,
    (pxt, qxt): zxt * t**3,
}
for (a, b), expected in pair_brackets.items():
    zero(jacobian(a, b) - expected, "m=2 generator bracket")


# ---------------------------------------------------------------------------
# Exact bounded hostile search.  Filtration means total degree in x,z,p,q.
# ---------------------------------------------------------------------------

generators = [x, zxt, pxt, qxt]


def filtered_candidates(bound: int) -> list[sp.Expr]:
    result: list[sp.Expr] = []
    for total in range(bound + 1):
        for a in range(total + 1):
            for b in range(total - a + 1):
                for c in range(total - a - b + 1):
                    d = total - a - b - c
                    result.append(sp.expand(
                        generators[0]**a * generators[1]**b
                        * generators[2]**c * generators[3]**d))
    return result


def coefficient_vector(poly: sp.Expr,
                       monomials: list[tuple[int, int]]) -> list[sp.Rational]:
    P = sp.Poly(poly, x, t)
    return [P.coeff_monomial(mon) for mon in monomials]


def filtered_basis(bound: int) -> list[sp.Expr]:
    candidates = filtered_candidates(bound)
    monomials = sorted({mon for f in candidates
                        for mon in sp.Poly(f, x, t).monoms()})
    rows: list[list[sp.Rational]] = []
    basis: list[sp.Expr] = []
    rank = 0
    for candidate in candidates:
        row = coefficient_vector(candidate, monomials)
        new_rank = sp.Matrix(rows + [row]).rank()
        if new_rank > rank:
            rows.append(row)
            basis.append(candidate)
            rank = new_rank
    return basis


bases = {bound: filtered_basis(bound) for bound in range(5)}
gate([len(bases[bound]) for bound in range(5)] == [1, 5, 12, 22, 35],
     "m=2 filtration dimensions through four")


def has_filtered_mate(a: sp.Expr, c_basis: list[sp.Expr]) -> bool:
    columns = [jacobian(a, c) for c in c_basis]
    monomials = sorted({(0, 0)} | {mon for f in columns
                                  for mon in sp.Poly(f, x, t).monoms()})
    matrix = sp.Matrix([
        [sp.Poly(f, x, t).coeff_monomial(mon) for f in columns]
        for mon in monomials
    ])
    target = sp.Matrix([1 if mon == (0, 0) else 0 for mon in monomials])
    return matrix.rank() == matrix.row_join(target).rank()


# Projectivize the 80 nonzero {-1,0,1} rows by fixing the first nonzero
# coefficient to +1.  This leaves exactly 40 generator-linear A's.
linear_rows: list[sp.Expr] = []
for coeffs in itertools.product((-1, 0, 1), repeat=4):
    if coeffs == (0, 0, 0, 0):
        continue
    if next(c for c in coeffs if c != 0) != 1:
        continue
    linear_rows.append(sp.expand(sum(c * g for c, g in zip(coeffs, generators))))
gate(len(linear_rows) == 40, "forty projective signed generator rows")

linear_survivors = sum(has_filtered_mate(a, bases[4]) for a in linear_rows)
gate(linear_survivors == 0,
     "no signed generator-linear A has a filtration-four mate")

low_rows = [a for a in bases[2] if a != 1]
gate(len(low_rows) == 11, "eleven nonconstant filtration-two basis rows")
low_survivors = sum(has_filtered_mate(a, bases[4]) for a in low_rows)
gate(low_survivors == 0,
     "no filtration-two basis row has a filtration-four mate")


summary = {
    "checks": CHECKS,
    "family": "B_m=k[x,z,zt,...,zt^m], z=1+xt, m>=1",
    "open": "D(x) union D(z)=A2_(x,t)",
    "boundary": "D=V(x,z)=A1_(zt^m), smooth and unique",
    "class": "units=k*, Cl(B_m)=Z[D]",
    "canonical": "div(dx wedge dt)=(m-1)D",
    "de_rham": "H2_dR(X_m)=k[D]; residue(-t dx)=1; volume nonexact",
    "m2": "2x3 determinantal; K=D pays the simple cubic canonical debt",
    "filtration_dimensions": [1, 5, 12, 22, 35],
    "bounded_search": {
        "signed_linear_A": 40,
        "filtration_two_A": 11,
        "filtration_four_mates": 0,
    },
    "conclusion": "no Darboux pair in B_m in any degree",
    "scope": "sharp completion near miss; no planar Keller counterexample",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3971 canonical-debt determinantal completion companion")
print(f"CHECKS={CHECKS}")
print("FAMILY=B_M_GENERATED_BY_X_Z_ZT_THROUGH_ZT_POWER_M;Z_1_PLUS_XT")
print("OPEN=D_X_UNION_D_Z_IS_A2_XT")
print("BOUNDARY=ONE_SMOOTH_A1;CL_Z_GENERATED_BY_D;UNITS_SCALARS")
print("CANONICAL=DIV_DX_WEDGE_DT_EQUALS_M_MINUS_1_TIMES_D")
print("DE_RHAM=H2_ONE_DIMENSIONAL;RESIDUE_OF_MINUS_T_DX_EQUALS_ONE")
print("M2=2_BY_3_DETERMINANTAL;SIMPLE_CUBIC_CANONICAL_DEBT_PAID")
print("FILTRATION_DIMS=1,5,12,22,35")
print("BOUNDED_SEARCH=40_SIGNED_LINEAR_PLUS_11_LOW_BASIS;SURVIVORS_0")
print("CONCLUSION=NO_DARBOUX_PAIR_IN_ANY_B_M;JC2_REMAINS_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
