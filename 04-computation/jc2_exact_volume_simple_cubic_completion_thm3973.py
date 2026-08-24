#!/usr/bin/env python3
"""Exact companion for THM-3973's exact-volume determinantal family."""

from __future__ import annotations

import hashlib
import itertools
import json
from math import ceil

from flint import fmpz_mat
import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(expression) == 0, message)


def integer_rank(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    return fmpz_mat(rows).rank()


x, t, z, p, y, u = sp.symbols("x t z p y u")


def jacobian(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(a, x) * sp.diff(b, t)
                     - sp.diff(a, t) * sp.diff(b, x))


# ---------------------------------------------------------------------------
# The one-parameter determinantal family B_n, n>=2.
# ---------------------------------------------------------------------------

for n in range(2, 7):
    zn = 1 + x**n * t
    pn = sp.expand(zn * t)
    yn = sp.expand(x ** (n - 1) * zn * t**2)

    zero(x**n * pn - zn * (zn - 1),
         f"n={n}: first determinantal relation")
    zero(pn * (zn - 1) - x * yn,
         f"n={n}: second determinantal relation")
    zero(zn * yn - x ** (n - 1) * pn**2,
         f"n={n}: third determinantal relation")

    matrix = sp.Matrix([[z, p, x], [x ** (n - 1) * p, y, z - 1]])
    minors = [sp.expand(matrix[:, [i, j]].det())
              for i, j in itertools.combinations(range(3), 2)]
    gate(minors == [sp.expand(z * y - x ** (n - 1) * p**2),
                    sp.expand(z * (z - 1) - x**n * p),
                    sp.expand(p * (z - 1) - x * y)],
         f"n={n}: exact 2-by-3 minor list")

    # The boundary chart is z(z-1)^2=x^(n+1)y and is smooth at x=z=0.
    F = z * (z - 1)**2 - x ** (n + 1) * y
    gate(sp.diff(F, z).subs({x: 0, z: 0}) == 1,
         f"n={n}: smooth boundary chart")

    dy_dt = sp.diff(yn, t)
    zero(x * dy_dt - (zn - 1) * (3 * zn - 1),
         f"n={n}: canonical order one")

    # beta=-t dx-(n-1)^-1 d(xt)=-dz/((n-1)x^(n-1)).
    beta_dx = -sp.Rational(n, n - 1) * t
    beta_dt = -sp.Rational(1, n - 1) * x
    zero(beta_dx + t + sp.Rational(1, n - 1) * t,
         f"n={n}: primitive dx coefficient")
    zero(beta_dt + sp.Rational(1, n - 1) * x,
         f"n={n}: primitive dt coefficient")
    zero(sp.diff(beta_dt, x) - sp.diff(beta_dx, t) - 1,
         f"n={n}: primitive differentiates to source volume")

    Fz = sp.diff(F, z)
    zero(sp.diff(F, x) + (n + 1) * x**n * y,
         f"n={n}: differentiated boundary dx coefficient")
    zero(sp.diff(F, y) + x ** (n + 1),
         f"n={n}: differentiated boundary dy coefficient")
    zero(Fz - (z - 1) * (3 * z - 1),
         f"n={n}: differentiated boundary dz coefficient")

    # The rational two-divisor compression has exact Jacobian one.
    un = x**n * t
    wn = 1 + 2 * un
    gn = un * (1 + un)
    P_rat = gn / ((n - 1) * wn**2)
    Q_rat = wn**3 / x ** (n - 1)
    zero(sp.together(jacobian(P_rat, Q_rat) - 1),
         f"n={n}: rational one-bracket compression")

    A_left_n = -wn * pn
    A_right_n = sp.Rational(6, n - 1) * x * pn
    zero(jacobian(A_left_n, x) - (1 + 6 * gn),
         f"n={n}: first bracket in uniform length-two identity")
    zero(jacobian(A_right_n, zn) + 6 * gn,
         f"n={n}: second bracket in uniform length-two identity")
    zero(jacobian(A_left_n, x) + jacobian(A_right_n, zn) - 1,
         f"n={n}: uniform length-two identity")


# ---------------------------------------------------------------------------
# Higher top generators give no new rings.
# ---------------------------------------------------------------------------

for m in range(2, 9):
    bezout_a, bezout_b, gcd = sp.gcdex(u + 1, u ** (m - 2), u)
    gate(gcd == 1, f"m={m}: coprime top-generator factors")
    zero(bezout_a * (u + 1) + bezout_b * u ** (m - 2) - 1,
         f"m={m}: top-generator Bezout identity")


# The exact negative-weight component exponents are independently checked.
for n in range(2, 8):
    for k in range(1, 22):
        feasible: list[tuple[int, int]] = []
        for a in range(0, ceil(k / n) + 3):
            for b in range(0, ceil(k / (n + 1)) + 3):
                if n * a + (n + 1) * b >= k:
                    feasible.append((a + 2 * b, a + b))
        gate(min(v[0] for v in feasible) == ceil(k / n),
             f"n={n},k={k}: u valuation floor")
        gate(min(v[1] for v in feasible) == ceil(k / (n + 1)),
             f"n={n},k={k}: u+1 valuation floor")


# The terminal scalar factor in the uniform two-by-two weight-support proof.
# These exact polynomial controls exercise every endpoint k=1 as well as
# several genuinely unequal ladders; the theorem gives the UFD derivation for
# arbitrary n,k.
Afac, Bfac, Lfac, Mfac = sp.symbols("Afac Bfac Lfac Mfac")
h_test = u * (u + 1) * (u**2 + u + 2)
K_test = u**2 + 2 * u + 3
for n in range(2, 7):
    for k in range(1, 6):
        upper = n * (k - 1) + 1
        f_test = Afac * h_test
        g_test = Bfac * h_test**k
        F_test = Lfac * K_test**upper
        G_test = Mfac * K_test
        scalar_row = (
            sp.diff(f_test, u) * G_test
            + n * f_test * sp.diff(G_test, u)
            - n * k * sp.diff(F_test, u) * g_test
            - upper * F_test * sp.diff(g_test, u)
        )
        terminal_factor = (
            K_test * sp.diff(h_test, u)
            + n * h_test * sp.diff(K_test, u)
        ) * (
            Afac * Mfac
            - k * upper * Lfac * Bfac * K_test ** (upper - 1)
            * h_test ** (k - 1)
        )
        zero(scalar_row - terminal_factor,
             f"n={n},k={k}: two-weight terminal factorization")
        gate(sp.degree(K_test * sp.diff(h_test, u)
                       + n * h_test * sp.diff(K_test, u), u) >= 1,
             f"n={n},k={k}: two-weight first factor nonconstant")


# ---------------------------------------------------------------------------
# Minimal member n=2: exact primitive and two-bracket collapse.
# ---------------------------------------------------------------------------

z2 = 1 + x**2 * t
p2 = sp.expand(z2 * t)
y2 = sp.expand(x * z2 * t**2)
w2 = 2 * z2 - 1

beta_a = -2 * p2 * w2
beta_b = 4 * x * p2
beta_c = -x * w2
zero(beta_a + beta_b * sp.diff(z2, x) + beta_c * sp.diff(p2, x)
     + 2 * t, "minimal primitive dx coefficient")
zero(beta_b * sp.diff(z2, t) + beta_c * sp.diff(p2, t)
     + x, "minimal primitive dt coefficient")
zero(jacobian(beta_a, x) + jacobian(beta_b, z2)
     + jacobian(beta_c, p2) - 1,
     "minimal three-term primitive bracket identity")

A_left = -w2 * p2
A_right = 6 * x * p2
zero(jacobian(A_left, x) + jacobian(A_right, z2) - 1,
     "minimal two-bracket collapse")


# ---------------------------------------------------------------------------
# Natural finite cubic control A=p, C=x+y on B_2.
# ---------------------------------------------------------------------------

Xtar, Ptar, Ctar, Wtar = sp.symbols("Xtar Ptar Ctar Wtar")
C2 = sp.expand(x + y2)
w_source = sp.expand(z2 - 1)
finite_cubic = (
    Xtar**3 - 2 * Ctar * Xtar**2
    + (Ctar**2 - Ptar - Ptar**3) * Xtar + Ptar * Ctar
)
zero(finite_cubic.subs({Xtar: x, Ptar: p2, Ctar: C2}),
     "natural finite cubic equation")
zero(p2 * y2 + x * y2**2 - x * p2**3,
     "natural finite cubic precursor identity")
zero(z2**2 - z2 - x**2 * p2,
     "natural finite cubic integral z equation")
zero(sp.together(z2 - (1 + x * (C2 - x) / p2)),
     "natural finite cubic generic recovery of z")

# The p=1 specialization has no constant or affine polynomial root.  The
# displayed coefficient rows are the exact exhaustive degree-at-most-one
# comparison used in the proof.
aroot, broot = sp.symbols("aroot broot")
specialized_root = sp.Poly(
    finite_cubic.subs({Ptar: 1, Xtar: aroot * Ctar + broot}), Ctar
)
zero(specialized_root.coeff_monomial(Ctar**3)
     - aroot * (aroot - 1)**2,
     "finite cubic affine-root top coefficient")
zero(specialized_root.as_expr().subs(aroot, 1)
     - (Ctar * (broot**2 - 1) + broot**3 - 2 * broot),
     "finite cubic affine-root endpoint")
gate(sp.Poly(finite_cubic.subs({Ptar: 1, Xtar: broot}), Ctar)
     .coeff_monomial(Ctar**2) == broot,
     "finite cubic constant-root C2 coefficient")

# Exact free normal basis {1,x,w}, w=z-1, and its discriminant/index ledger.
zero(x**2 - (C2 * x - p2 * w_source),
     "finite cubic free basis x2 row")
zero(x * w_source
     - (C2 * w_source + C2 - (1 + p2**2) * x),
     "finite cubic free basis xw row")
zero(w_source**2
     - (p2 * C2 * x - (1 + p2**2) * w_source),
     "finite cubic free basis w2 row")

Mx = sp.Matrix([
    [0, 0, Ctar],
    [1, Ctar, -(1 + Ptar**2)],
    [0, -Ptar, Ctar],
])
Mw = sp.Matrix([
    [0, Ctar, 0],
    [0, -(1 + Ptar**2), Ptar * Ctar],
    [1, Ctar, -(1 + Ptar**2)],
])
trace_basis = [sp.eye(3), Mx, Mw]
normal_discriminant = sp.factor(sp.Matrix([
    [sp.trace(left * right) for right in trace_basis]
    for left in trace_basis
]).det())
Delta2 = (
    4 * Ctar**4 * Ptar - 8 * Ctar**2 * Ptar**4
    + 20 * Ctar**2 * Ptar**2 + Ctar**2
    + 4 * Ptar**7 + 12 * Ptar**5 + 12 * Ptar**3 + 4 * Ptar
)
zero(normal_discriminant - Delta2,
     "finite cubic normal-basis discriminant")
zero(sp.discriminant(finite_cubic, Xtar) - Ptar**2 * Delta2,
     "finite cubic monogenic index-square discriminant")
gate(Delta2.subs(Ptar, 0) == Ctar**2,
     "finite cubic p=0 is generically unramified")

# All three p=0,C!=0 normalization addresses satisfy the specialized table.
fiber_rows = [
    Xtar**2 - Ctar * Xtar + Ptar * Wtar,
    Xtar * Wtar - Ctar * Wtar - Ctar + (1 + Ptar**2) * Xtar,
    Wtar**2 - Ptar * Ctar * Xtar + (1 + Ptar**2) * Wtar,
]
for address in ((0, -1), (Ctar, -1), (Ctar, 0)):
    for row in fiber_rows:
        zero(row.subs({Ptar: 0, Xtar: address[0], Wtar: address[1]}),
             "finite cubic p=0 normalization address")

R2 = sp.expand(x * (2 * z2 - 1) + y2)
J2 = sp.factor(jacobian(p2, C2))
zero(x * J2 + R2, "finite cubic bracket numerator")
zero(J2 + t * (t**2 + 2) * x**2 + t**2 + 1,
     "finite cubic irreducible source ramification equation")
gate(sp.rem(t**2 + 1, t, domain=sp.QQ) == 1,
     "finite cubic ramification nonsquare odd t valuation")
zero(R2 - x * ((2 * z2 - 1) + p2**2 / z2),
     "finite cubic L1 residual factor")
gate(sp.factor(J2.subs(x, 0)) == -(1 + t**2),
     "finite cubic residual ramification meets L1")


# ---------------------------------------------------------------------------
# All-degree no-mate controls for the four canonical generators.
# ---------------------------------------------------------------------------

f_coeffs = sp.symbols("f0:8")
Fpoly = sum(f_coeffs[i] * u**i for i in range(len(f_coeffs)))

for n in range(2, 7):
    un = x**n * t
    zn = 1 + un
    pn = t * zn
    yn = x ** (n - 1) * zn * t**2

    # z: J(A,z)=x^(n-1)(xA_x-ntA_t).
    Atest = sum(sp.symbols(f"a{n}_0:{5}")[i] * x**i * t**i
                for i in range(5))
    zero(jacobian(Atest, zn)
         - x ** (n - 1) * (x * sp.diff(Atest, x)
                            - n * t * sp.diff(Atest, t)),
         f"n={n}: z Euler derivation")

    # p: only the weight-one row xF(u) can contribute a constant.
    Ap = x * Fpoly.subs(u, un)
    p_ode = (n * u * (1 + u) * sp.diff(Fpoly, u)
             + (1 + 2 * u) * Fpoly)
    zero(jacobian(Ap, pn) - p_ode.subs(u, un),
         f"n={n}: p mate ODE")

    # y: only weight two x^2F(u), and its bracket is divisible by u^(m-1)
    # with m=2 in the canonical presentation.
    Ay = x**2 * Fpoly.subs(u, un)
    y_ode = u * ((n + 1) * u * (1 + u) * sp.diff(Fpoly, u)
                 + 2 * (2 + 3 * u) * Fpoly)
    zero(jacobian(Ay, yn) - y_ode.subs(u, un),
         f"n={n}: y mate ODE")

    lam = sp.symbols(f"lambda_{n}", nonzero=True)
    corrected_y = sp.expand(yn + lam * pn)
    critical_address = {x: 1 / lam, t: -lam**n}
    zero(sp.together(sp.diff(corrected_y, x).subs(critical_address)),
         f"n={n}: y+lambda p critical x derivative")
    zero(sp.together(sp.diff(corrected_y, t).subs(critical_address)),
         f"n={n}: y+lambda p critical t derivative")


# ---------------------------------------------------------------------------
# Minimal exact generator filtration and sparse Darboux search.
# ---------------------------------------------------------------------------

generators = [x, z2, p2, y2]


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


def filtered_basis(bound: int) -> list[sp.Expr]:
    candidates = filtered_candidates(bound)
    monomials = sorted({mon for f in candidates
                        for mon in sp.Poly(f, x, t).monoms()})
    rows: list[list[int]] = []
    basis: list[sp.Expr] = []
    rank = 0
    for candidate in candidates:
        poly = sp.Poly(candidate, x, t)
        row = [int(poly.coeff_monomial(mon)) for mon in monomials]
        new_rank = integer_rank(rows + [row])
        if new_rank > rank:
            rows.append(row)
            basis.append(candidate)
            rank = new_rank
    return basis


bases = {bound: filtered_basis(bound) for bound in range(5)}
gate([len(bases[bound]) for bound in range(5)] == [1, 5, 14, 28, 47],
     "minimal filtration dimensions through four")

# There is no pair when both entries are generator-linear, for arbitrary
# coefficients rather than only a finite coefficient box.
a_coeffs = sp.symbols("a0:4")
c_coeffs = sp.symbols("c0:4")
A_linear = sum(a * g for a, g in zip(a_coeffs, generators))
C_linear = sum(c * g for c, g in zip(c_coeffs, generators))
linear_rows = sp.Poly(jacobian(A_linear, C_linear) - 1, x, t).coeffs()
linear_gb = sp.groebner(linear_rows, *a_coeffs, *c_coeffs,
                        order="grevlex")
gate([poly.as_expr() for poly in linear_gb.polys] == [1],
     "no arbitrary generator-linear Darboux pair")

# Every generator-linear function which is nonconstant on the boundary has
# an affine critical point.  The resultant records the complete parameter
# split; its nonzero roots are genuine because the t-leading rows stay
# nonzero whenever x is nonzero.
Apar, Bpar, Cpar = sp.symbols("Apar Bpar Cpar")
boundary_linear = y2 + Apar * x + Bpar * z2 + Cpar * p2
critical_resultant = sp.factor(sp.resultant(
    sp.diff(boundary_linear, x), sp.diff(boundary_linear, t), t))
Hcritical = (
    3 * Bpar**3 * x**7
    + (9 * Apar**2 - 6 * Apar * Bpar * Cpar) * x**6
    - Bpar**2 * Cpar * x**5
    + (2 * Apar * Cpar**2 + 3 * Bpar**2) * x**4
    + (4 * Apar * Cpar - 3 * Bpar * Cpar**2) * x**3
    + (-4 * Apar + 2 * Bpar * Cpar) * x**2
    + Cpar**3 * x - Cpar**2
)
zero(critical_resultant - 3 * x**3 * Hcritical,
     "boundary-linear critical resultant")
gate(sp.Poly(sp.diff(boundary_linear, x), t).LC() == 3 * x**2,
     "boundary-linear x-derivative has stable t-leading row")
gate(sp.Poly(sp.diff(boundary_linear, t), t).LC() == 3 * x**3,
     "boundary-linear t-derivative has stable t-leading row")

# Precompute every J(A_i,C_j) coefficient vector.  This makes the 7020-row
# exact integer-rank census fast enough to replay twice.
low_basis = bases[2][1:]
mate_basis = bases[4]
bracket_polys = [[sp.Poly(jacobian(a, c), x, t) for c in mate_basis]
                 for a in low_basis]
search_monomials = sorted({(0, 0)} | {
    mon for row in bracket_polys for poly in row for mon in poly.monoms()
})
tensor = [[[
    int(bracket_polys[i][j].coeff_monomial(mon))
    for j in range(len(mate_basis))
] for mon in search_monomials] for i in range(len(low_basis))]
target = [1 if mon == (0, 0) else 0 for mon in search_monomials]


def has_mate(indices: tuple[int, ...], signs: tuple[int, ...]) -> bool:
    data = [[
        sum(sign * tensor[index][row][column]
            for index, sign in zip(indices, signs))
        for column in range(len(mate_basis))
    ] for row in range(len(search_monomials))]
    rank = integer_rank(data)
    augmented = [row + [entry] for row, entry in zip(data, target)]
    return rank == integer_rank(augmented)


gate(len(low_basis) == 13, "thirteen nonconstant filtration-two basis rows")
for i in range(len(low_basis)):
    gate(not has_mate((i,), (1,)),
         f"filtration-two basis row {i} has no filtration-four mate")

sparse_rows = 0
for support in (2, 3, 4):
    for indices in itertools.combinations(range(len(low_basis)), support):
        for tail_signs in itertools.product((-1, 1), repeat=support - 1):
            signs = (1,) + tail_signs
            sparse_rows += 1
            gate(not has_mate(indices, signs),
                 "signed sparse filtration-two row has no mate")
gate(sparse_rows == 7020, "complete signed support-two-through-four count")


summary = {
    "checks": CHECKS,
    "family": "B_n=k[x,z,p,y], z=1+x^n t, p=zt, y=x^(n-1)zt^2, n>=2",
    "presentation": "2x2 minors of [[z,p,x],[x^(n-1)p,y,z-1]]",
    "open_boundary": "U=A2_(x,t); D=V(x,z)=A1_y",
    "passport": "units=k*; Cl=Z[D]; div(dx wedge dt)=D",
    "exact_volume": "beta=-dz/((n-1)x^(n-1)) is regular and d beta=dx wedge dt",
    "grading": "B_-k=x^-k u^ceil(k/n)(u+1)^ceil(k/(n+1)) k[u]",
    "collapse": "all higher top generators give the same B_n",
    "brackets": "constant has bracket length at most two; every two-by-two weight support impossible",
    "rational_pair": "P=u(u+1)/((n-1)(1+2u)^2); Q=(1+2u)^3/x^(n-1)",
    "finite_cubic": "B_2 finite free rank3 over k[p,x+y], normal basis 1,x,z-1; residual interior ramification",
    "minimal_search": "dims 1,5,14,28,47; linear ideal [1]; 7033 sparse rows have no F4 mate",
    "scope": "positive completion passport; natural finite cubic non-Keller; finite Keller pair open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3973 exact-volume simple-cubic completion companion")
print(f"CHECKS={CHECKS}")
print("FAMILY=B_N_X_Z_P_Y;Z_1_PLUS_XN_T;P_ZT;Y_XN_MINUS1_ZT2;N_GE_2")
print("PRESENTATION=MINORS_OF_Z_P_X__XN_MINUS1P_Y_Z_MINUS1")
print("OPEN_BOUNDARY=U_A2_XT;D_A1_Y;UNITS_SCALARS;CL_ZD")
print("CANONICAL=DIV_SOURCE_VOLUME_D;KAPPA_1")
print("EXACT_VOLUME=BETA_MINUS_DZ_OVER_N_MINUS1_XN_MINUS1;RHO_0")
print("TOP_GENERATORS=ALL_M_GE_2_COLLAPSE_TO_B_N")
print("GRADING=EXACT_NEGATIVE_WEIGHT_U_AND_U_PLUS1_VALUATION_FLOORS")
print("BRACKETS=LENGTH_AT_MOST_2;NO_TWO_BY_TWO_WEIGHT_SUPPORT_OR_CANONICAL_GENERATOR_PAIR")
print("RATIONAL_COMPRESSION=J_P_Q_1;DENOMINATORS_X_AND_2Z_MINUS1")
print("FINITE_CUBIC=B2_FREE_RANK3_OVER_K_P_X_PLUS_Y;NORMAL_NONMONOGENIC;INTERIOR_RAMIFICATION")
print("MINIMAL_FILTRATION_DIMS=1,5,14,28,47")
print("MINIMAL_SEARCH=ARBITRARY_LINEAR_EMPTY;7033_LOW_SPARSE_ROWS_NO_F4_MATE")
print("CONCLUSION=POSITIVE_COMPLETION_PASSPORT;NATURAL_FINITE_CUBIC_NONKELLER;FINITE_KELLER_MAP_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
