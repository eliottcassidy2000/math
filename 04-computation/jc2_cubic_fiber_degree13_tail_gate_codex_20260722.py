#!/usr/bin/env python3
"""Exact referee for THM-2110 (cubic source-fiber degree 13 gate).

Universe: exact rational polynomial arithmetic in Q[p,q,z,a,Y].  The script
checks the degree-thirteen Faber identity, the finite centering-pole resultant,
the upper-Newton support separation from every lower reduced Faber term, and
the explicit Bezout certificate for the degree-at-infinity primitive.  The
all-n adjacent-tail lemma is proved symbolically in the theorem from its
coefficient recurrence; the finite sweep here is only a hostile control.

Reproduce:
    python3 04-computation/jc2_cubic_fiber_degree13_tail_gate_codex_20260722.py
    python3 -O 04-computation/jc2_cubic_fiber_degree13_tail_gate_codex_20260722.py
"""

from __future__ import annotations

import sympy as sp


z, p, q, dp, dq, a, Y = sp.symbols("z p q dp dq a Y")


def faber(n: int) -> sp.Expr:
    """Polynomial part at z=infinity of (z^3+pz+q)^(n/3)."""
    return sp.expand(
        sum(
            sp.binomial(sp.Rational(n, 3), i + j)
            * sp.binomial(i + j, i)
            * p**i
            * q**j
            * z ** (n - 2 * i - 3 * j)
            for i in range(n // 2 + 1)
            for j in range(n // 3 + 1)
            if 2 * i + 3 * j <= n
        )
    )


def tail(n: int, offset: int) -> sp.Expr:
    """Three times [z^(-offset)] (z^3+pz+q)^(n/3)."""
    weight = n + offset
    return sp.expand(
        sum(
            3
            * sp.binomial(sp.Rational(n, 3), i + j)
            * sp.binomial(i + j, i)
            * p**i
            * q**j
            for i in range(weight // 2 + 1)
            for j in range(weight // 3 + 1)
            if 2 * i + 3 * j == weight
        )
    )


def dx(expr: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(expr, p) * dp + sp.diff(expr, q) * dq)


def L(expr: sp.Expr) -> sp.Expr:
    return sp.expand(
        (dp * z + dq) * sp.diff(expr, z)
        - (3 * z**2 + p)
        * (sp.diff(expr, p) * dp + sp.diff(expr, q) * dq)
    )


E13 = faber(13)
Phi13 = tail(13, 1)
R13 = tail(13, 2)

expected_phi13 = (
    sp.Rational(65, 6561) * p * (p**6 - 63 * p**3 * q**2 + 189 * q**4)
)
expected_r13 = (
    sp.Rational(91, 6561) * q * (5 * p**6 - 60 * p**3 * q**2 + 27 * q**4)
)

assert sp.expand(Phi13 - expected_phi13) == 0
assert sp.expand(R13 - expected_r13) == 0
assert sp.expand(L(E13) - z * dx(Phi13) - dx(R13)) == 0
print("degree-13 Faber/Gauss--Manin identity: PASS")


# If h has a centering pole and p ~ a h^2, q ~ -(1+a)h^3, the
# degree-thirteen leading coefficients are the following primitive polynomials.
e_boundary = sp.factor(E13.subs({z: 1, p: a, q: -(1 + a)}))
phi_boundary = sp.factor(Phi13.subs({p: a, q: -(1 + a)}))
e_primitive = 13 * a**6 - 78 * a**5 - 351 * a**4 - 468 * a**3 - 234 * a**2 + 27
phi_primitive = a * (
    a**6
    - 63 * a**5
    + 63 * a**4
    + 693 * a**3
    + 1134 * a**2
    + 756 * a
    + 189
)
assert sp.expand(e_boundary - sp.Rational(35, 6561) * e_primitive) == 0
assert sp.expand(phi_boundary - sp.Rational(65, 6561) * phi_primitive) == 0
boundary_resultant = int(sp.resultant(e_primitive, phi_primitive, a))
assert boundary_resultant == -(3**21) * (17**5) * (23**2)
assert sp.Poly(sp.gcd(e_primitive, phi_primitive), a).degree() == 0
assert E13.subs({z: 1, p: 0, q: -1}) == sp.Rational(35, 243)
print(f"degree-13 centering-pole resultant: {boundary_resultant} PASS")


# Every lower Phi support monomial is strictly componentwise dominated by a
# monomial of Phi_13 when deg(p),deg(q)>0.  Hence lower Faber representatives
# cannot join the leading degree cancellation.
reduced_degrees = (1, 2, 4, 5, 7, 8, 10, 11, 13)


def support(expr: sp.Expr) -> set[tuple[int, int]]:
    return {
        monomial
        for monomial, coefficient in sp.Poly(expr, p, q).terms()
        if coefficient != 0
    }


top_support = support(Phi13)
for m in reduced_degrees[:-1]:
    for i, j in support(tail(m, 1)):
        assert any(I >= i and J >= j and (I, J) != (i, j) for I, J in top_support)
print("lower-Phi upper-Newton separation at n=13: PASS")


# On the only internal top-edge tie, deg(p):deg(q)=2:3.  With
# Y=q_lead^2/p_lead^3, Phi_13 and R_13 have the two displayed factors.
# This exact Bezout identity is the degree-thirteen instance of the theorem's
# all-n adjacent-tail noncollision lemma.
F = 1 - 63 * Y + 189 * Y**2
G = 5 - 60 * Y + 27 * Y**2
assert sp.expand((63 * Y - 134) * F + (105 - 441 * Y) * G) == 391
assert sp.resultant(F, G, Y) == -(3**5) * (17**2) * 23
print("degree-13 infinity primitive Bezout certificate: PASS")


# Hostile finite controls for the all-n coefficient-recurrence lemma.  This
# does not replace its UFD proof: it checks both residue classes through n=61,
# including the repeated-root discriminant points p=-3,q=2 and its scaling.
for n in range(1, 62):
    if sp.gcd(n, 3) != 1:
        continue
    phi_n = tail(n, 1)
    r_n = tail(n, 2)
    if n % 2:
        # Remove the forced p/q monomial prefactors and reduce to Y=q^2/p^3.
        phi_y = sp.Poly(phi_n.subs({p: 1, q**2: Y}), q).as_expr()
        r_y = sp.Poly(r_n.subs({p: 1, q**2: Y}), q).as_expr()
        if sp.Poly(phi_y, q).degree() == 1:
            phi_y = sp.cancel(phi_y / q)
        if sp.Poly(r_y, q).degree() == 1:
            r_y = sp.cancel(r_y / q)
        assert sp.Poly(sp.gcd(sp.together(phi_y), sp.together(r_y)), Y).degree() == 0
    assert not (phi_n.subs({p: -3, q: 2}) == 0 and r_n.subs({p: -3, q: 2}) == 0)
print("odd-n projective-gcd and all-n discriminant controls through n=61: PASS")


# The excluded common-zero exception is real: at the trivial cubic z^3 both
# adjacent tails vanish.  The theorem explicitly excludes (p,q)=(0,0).
assert Phi13.subs({p: 0, q: 0}) == 0
assert R13.subs({p: 0, q: 0}) == 0
print("trivial-cubic boundary control: PASS")

print("RESULT: PASS")
