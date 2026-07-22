#!/usr/bin/env python3
"""Exact referee for THM-2084 (source-fiber cubic Gauss--Manin gate).

Universe: Q[p,q,z,dp,dq] and the four Hermite representatives of reduced
degrees 4, 5, 7, and 8.  This is an identity check, not a bounded search.

Reproduce:
    python3 04-computation/jc2_cubic_fiber_gauss_manin_gate_codex_20260722.py
    python3 -O 04-computation/jc2_cubic_fiber_gauss_manin_gate_codex_20260722.py
"""

from __future__ import annotations

import sympy as sp


z, p, q, dp, dq, a = sp.symbols("z p q dp dq a")


def dx(expr: sp.Expr) -> sp.Expr:
    """Formal x derivative for an expression in p,q."""
    return sp.expand(sp.diff(expr, p) * dp + sp.diff(expr, q) * dq)


def L(expr: sp.Expr) -> sp.Expr:
    """Hamiltonian field of z^3+pz+q, with z held in partial_x."""
    return sp.expand(
        (dp * z + dq) * sp.diff(expr, z)
        - (3 * z**2 + p) * (sp.diff(expr, p) * dp + sp.diff(expr, q) * dq)
    )


E = {
    1: z,
    2: z**2 + sp.Rational(2, 3) * p,
    4: (
        z**4
        + sp.Rational(4, 3) * p * z**2
        + sp.Rational(4, 3) * q * z
        + sp.Rational(2, 9) * p**2
    ),
    5: (
        z**5
        + sp.Rational(5, 3) * p * z**3
        + sp.Rational(5, 3) * q * z**2
        + sp.Rational(5, 9) * p**2 * z
        + sp.Rational(10, 9) * p * q
    ),
    7: (
        z**7
        + sp.Rational(7, 3) * p * z**5
        + sp.Rational(7, 3) * q * z**4
        + sp.Rational(14, 9) * p**2 * z**3
        + sp.Rational(28, 9) * p * q * z**2
        + (sp.Rational(14, 81) * p**3 + sp.Rational(14, 9) * q**2) * z
        + sp.Rational(14, 27) * p**2 * q
    ),
    8: (
        z**8
        + sp.Rational(8, 3) * p * z**6
        + sp.Rational(8, 3) * q * z**5
        + sp.Rational(20, 9) * p**2 * z**4
        + sp.Rational(40, 9) * p * q * z**3
        + (sp.Rational(40, 81) * p**3 + sp.Rational(20, 9) * q**2) * z**2
        + sp.Rational(40, 27) * p**2 * q * z
        - sp.Rational(10, 243) * p**4
        + sp.Rational(40, 27) * p * q**2
    ),
}

Phi = {
    1: p,
    2: 2 * q,
    4: sp.Rational(4, 3) * p * q,
    5: -sp.Rational(5, 27) * p**3 + sp.Rational(5, 3) * q**2,
    7: -sp.Rational(7, 81) * p**4 + sp.Rational(14, 9) * p * q**2,
    8: -sp.Rational(40, 81) * p**3 * q + sp.Rational(40, 27) * q**3,
}

R = {
    1: q,
    2: -sp.Rational(1, 3) * p**2,
    4: -sp.Rational(4, 27) * p**3 + sp.Rational(2, 3) * q**2,
    5: -sp.Rational(5, 9) * p**2 * q,
    7: -sp.Rational(28, 81) * p**3 * q + sp.Rational(14, 27) * q**3,
    8: sp.Rational(8, 243) * p**5 - sp.Rational(20, 27) * p**2 * q**2,
}


print("Hermite/Gauss--Manin identities")
for n in E:
    residual = sp.expand(L(E[n]) - z * dx(Phi[n]) - dx(R[n]))
    assert residual == 0
    print(f"  E_{n}: PASS")


# Verify the rank-two connection in Q[p,q,t][z]/(z^3+pz+q-t).
t = sp.symbols("t")
Ax, Bx, Cx, A0, B0, C0 = sp.symbols("Ax Bx Cx A0 B0 C0")
connection = {
    0: -p * Ax + 3 * (q - t) * Bx + dq * B0,
    1: 2 * p * Bx + 3 * (q - t) * Cx + dp * B0 + 2 * dq * C0,
    2: -3 * Ax + 2 * p * Cx + 2 * dp * C0,
}

# Direct reduction of L(A+Bz+Cz^2), treating A,B,C as functions of (x,t).
raw = (
    -(3 * z**2 + p) * (Ax + Bx * z + Cx * z**2)
    + (dp * z + dq) * (B0 + 2 * C0 * z)
)
raw = sp.Poly(sp.expand(raw), z)
reduced = 0
for (power,), coeff in raw.terms():
    term = coeff * z**power
    while sp.Poly(term, z).degree() >= 3:
        poly = sp.Poly(term, z)
        top = poly.degree()
        lc = poly.LC()
        term = sp.expand(term - lc * z ** (top - 3) * (z**3 + p * z + q - t))
    reduced += term
reduced = sp.Poly(sp.expand(reduced), z)
for k in range(3):
    assert sp.expand(reduced.coeff_monomial(z**k) - connection[k]) == 0
print("rank-two connection: PASS")


# Boundary polynomials for a possible pole h: z=h, p=a h^2,
# q=-(1+a)h^3+regular.  Their exact nonzero resultants are the pole gate.
boundary_expected = {
    4: 48,
    5: -11025,
    7: 13070456784,
    8: 21902400000000,
}

print("boundary resultant gates")
for n in (4, 5, 7, 8):
    e_num = sp.together(E[n].subs({z: 1, p: a, q: -(1 + a)})).as_numer_denom()[0]
    phi_num = sp.together(Phi[n].subs({p: a, q: -(1 + a)})).as_numer_denom()[0]
    resultant = int(sp.resultant(e_num, phi_num, a))
    assert resultant == boundary_expected[n]
    # SymPy keeps the harmless common integer content 5 in the n=8 row.
    assert sp.Poly(sp.gcd(e_num, phi_num), a).degree() == 0
    print(f"  n={n}: resultant={resultant} PASS")


# Degree-at-infinity noncancellation on the only balanced slopes.
P0, Q0 = sp.symbols("P0 Q0", nonzero=True)

# n=5: q^2=p^3/9, and R_5 has the unique p^2 q leading term.
assert sp.expand(R[5].subs(q**2, p**3 / 9)) != 0

# n=7: q^2=p^3/18 leaves the weight-nine R_7 term nonzero.
R7_on_branch = sp.rem(sp.Poly(R[7], q), sp.Poly(q**2 - p**3 / 18, q)).as_expr()
assert sp.factor(R7_on_branch) == -sp.Rational(77, 243) * p**3 * q

# n=8: q^2=p^3/3 leaves the weight-ten R_8 term nonzero.
assert sp.factor(R[8].subs(q**2, p**3 / 3)) == -sp.Rational(52, 243) * p**5

# The extra n=8 Newton edge b/a=1 is killed by the unique p^5 term of R_8.
assert sp.Poly(R[8], p, q).coeff_monomial(p**5) == sp.Rational(8, 243)
print("degree-at-infinity gates: PASS")


# First still-open degree: derive E_10 from the all-n polynomial-part formula
# and verify the two covariants recorded as the next frontier.
def faber(n: int) -> sp.Expr:
    out = 0
    for i in range(n // 2 + 1):
        for j in range(n // 3 + 1):
            if 2 * i + 3 * j <= n:
                out += (
                    sp.binomial(sp.Rational(n, 3), i + j)
                    * sp.binomial(i + j, i)
                    * p**i
                    * q**j
                    * z ** (n - 2 * i - 3 * j)
                )
    return sp.expand(out)


Phi10 = sp.Rational(70, 243) * p * q * (-p**3 + 6 * q**2)
R10 = sp.Rational(35, 2187) * (p**6 - 36 * p**3 * q**2 + 27 * q**4)
assert sp.expand(L(faber(10)) - z * dx(Phi10) - dx(R10)) == 0
print("n=10 frontier covariant: PASS")


# Positive tame control: P=y^3+x, Q=y has constant Jacobian 1 and reduced
# complementary fiber degree one, exactly the allowed endpoint.
x, y = sp.symbols("x y")
P_control = y**3 + x
Q_control = y
J_control = sp.diff(P_control, x) * sp.diff(Q_control, y) - sp.diff(P_control, y) * sp.diff(Q_control, x)
assert J_control == 1
print("tame cubic/linear control: PASS")

print("RESULT: PASS")
