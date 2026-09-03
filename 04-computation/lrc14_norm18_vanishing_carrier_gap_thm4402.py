#!/usr/bin/env python3
"""Exact primary certificate for THM-4402.

The certificate treats the infinite family

    w=(1,m,16m-1), c=C=(1,-16,1), m=6t+5, t>=2.

It checks the algebraic relation and an explicit carrier lift, exhausts every
possible relation of l1 norm at most 17 uniformly in m, and verifies from the
six-term carrier formula that the component length is 3/[7(16m-1)].
"""

from __future__ import annotations

from itertools import product

import sympy as sp


def check(condition: bool, label: str) -> None:
    """Optimization-safe certificate check."""
    if not bool(condition):
        raise RuntimeError(f"FAIL: {label}")


def l1(vector: tuple[int, int, int]) -> int:
    return sum(abs(entry) for entry in vector)


def cross(left, right):
    return (
        sp.expand(left[1] * right[2] - left[2] * right[1]),
        sp.expand(left[2] * right[0] - left[0] * right[2]),
        sp.expand(left[0] * right[1] - left[1] * right[0]),
    )


def dot(left, right):
    return sp.expand(sum(a * b for a, b in zip(left, right)))


def nonnegative_rational_for_m_ge_17(expression: sp.Expr, m: sp.Symbol) -> bool:
    """Certify a rational expression by coefficient positivity after m=u+17."""
    u = sp.symbols("u", nonnegative=True, integer=True)
    numerator, denominator = sp.fraction(sp.cancel(expression.subs(m, u + 17)))
    numerator_poly = sp.Poly(numerator, u, domain=sp.QQ)
    denominator_poly = sp.Poly(denominator, u, domain=sp.QQ)
    numerator_coefficients = numerator_poly.all_coeffs()
    denominator_coefficients = denominator_poly.all_coeffs()
    return (
        all(coefficient >= 0 for coefficient in numerator_coefficients)
        and any(coefficient > 0 for coefficient in numerator_coefficients)
        and all(coefficient >= 0 for coefficient in denominator_coefficients)
        and any(coefficient > 0 for coefficient in denominator_coefficients)
    )


m = sp.symbols("m", integer=True, positive=True)
t = sp.symbols("t", integer=True, nonnegative=True)
W = 16 * m - 1
w = (sp.Integer(1), m, W)
c = (sp.Integer(1), sp.Integer(-16), sp.Integer(1))
n = (sp.Integer(0), sp.Integer(1), sp.Integer(16))
C = cross(w, n)
r = sp.Rational(3, 14)

m_family = 6 * (t + 2) + 5
W_family = sp.expand(16 * m_family - 1)
check(m_family.subs(t, 0) == 17, "family begins at m=17")
check(sp.Poly(m_family, t).all_coeffs() == [6, 17], "m positivity")
check(sp.Poly(W_family - m_family, t).all_coeffs() == [90, 254],
      "strictly increasing distinct speeds")
check(sp.expand(m_family - (2 * (3 * t + 8) + 1)) == 0,
      "m is identically odd")
check(sp.expand(W_family - (2 * (48 * t + 135) + 1)) == 0,
      "W is identically odd")
check(sp.expand(m_family - (3 * (2 * t + 5) + 2)) == 0,
      "m is identically 2 modulo 3")
check(sp.expand(W_family - (3 * (32 * t + 90) + 1)) == 0,
      "W is identically 1 modulo 3")

check(dot(c, w) == 0, "relation c dot w")
check(C == c, "explicit raw-carrier lift w cross n=c")
check(dot(c, n) == 0, "zero defect c dot n")
check(l1(tuple(int(entry) for entry in c)) == 18, "relation norm")
check(sp.gcd_list(c) in (1, -1), "relation primitive")

# When m=6t+5, the speed residues are (1,2,1), the lift residues are
# (0,1,1), and -w_i^{-1} n_i gives owners (0,1,2).
speed_residues = (1, 2, 1)
lift_residues = (0, 1, 1)
owner_residues = tuple(
    (-pow(speed, -1, 3) * lift) % 3
    for speed, lift in zip(speed_residues, lift_residues)
)
check(owner_residues == (0, 1, 2), "distinct owner residues")
check(tuple(int(entry % 3) for entry in c) == (1, 2, 1),
      "ternary-unit full support")

# A relation a+b*m+d*(16m-1)=0 is
# (a-d)+(b+16d)m=0.  Enumerating the finite coefficient l1-ball through 17
# therefore audits all m at once, not a finite speed window.
coefficient_vectors = 0
nonidentity_integer_roots = set()
for a, b, d in product(range(-17, 18), repeat=3):
    if (a, b, d) == (0, 0, 0) or l1((a, b, d)) > 17:
        continue
    coefficient_vectors += 1
    intercept = a - d
    slope = b + 16 * d
    if slope == 0:
        check(intercept != 0, "no identity relation below norm 18")
        continue
    root = sp.Rational(-intercept, slope)
    if root.q == 1:
        nonidentity_integer_roots.add(int(root))
        check(root < 17, "no relation of norm at most 17 for m>=17")

# Odd speeds force every integer relation to have even l1 norm, so the
# exhaustion through 17 and the displayed norm-18 relation prove minimality.
check(all((a + b + d) % 2 == l1((a, b, d)) % 2
          for a, b, d in product(range(-2, 3), repeat=3)),
      "absolute-value parity identity")

components = (
    2 * r,
    2 * r / m,
    2 * r / W,
    r + r / m - abs(int(C[2])) / m,
    r + r / W - abs(int(C[1])) / W,
    r / m + r / W - abs(int(C[0])) / (m * W),
)
target = sp.factor(sp.Rational(3, 7) / W)
differences = tuple(sp.factor(component - target) for component in components)
check(components[2] == target, "third interval gives target length")
check(all(nonnegative_rational_for_m_ge_17(difference, m)
          for difference in differences if difference != 0),
      "target is the minimum of all six carrier terms for m>=17")
check(sp.limit(target, m, sp.oo) == 0, "component length tends to zero")

print("THM-4402 LRC14 NORM-18 VANISHING LIVE-CARRIER GAP")
print("family: m=6t+5, t>=2; w=(1,m,16m-1)")
print("primitive distinct positive odd ternary-unit speed triple: PASS")
print("relation/carrier: c=C=(1,-16,1), ||c||_1=18")
print("explicit lift: n=(0,1,16), w cross n=C, c dot n=0")
print("speed residues mod 3:", speed_residues)
print("owner residues mod 3:", owner_residues)
print("coefficient vectors exhausted through l1=17:", coefficient_vectors)
print("largest nonidentity integer root m:", max(nonidentity_integer_roots))
print("all relation norms are even for odd speeds: PASS")
print("carrier formula terms:")
for index, component in enumerate(components, start=1):
    print(f"  q{index} =", sp.factor(component))
print("differences q_i-3/[7(16m-1)]:")
for index, difference in enumerate(differences, start=1):
    print(f"  d{index} =", difference)
print("L_w(C)=", target)
print("limit as m tends to infinity: 0")
print("scope: no uniform positive lower gap for a live component in the minimal norm-18 shell")
print("NO_CLAIM=LRC14_or_any_scale-dependent_owner-conditioned_bound")
print("RESULT=PASS")
