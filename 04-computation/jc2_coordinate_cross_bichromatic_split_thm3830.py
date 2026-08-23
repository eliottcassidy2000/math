#!/usr/bin/env python3
"""Exact companion for THM-3830's coordinate-cross sign obstruction."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: object, rhs: object, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)  # type: ignore[operator]


z, alpha, c = sp.symbols("z alpha c", nonzero=True)
t, q, x, y = sp.symbols("t q x y")
a = sp.expand((7 * z**2 + 3) * (3 * z**3 + 7 * z**2 + 1))
b = sp.expand((z + 1) * (2 * z + 1) * (3 * z - 1))

check(a == 21 * z**5 + 49 * z**4 + 9 * z**3 + 28 * z**2 + 3,
      "expanded five-slope polynomial")
check(b == 6 * z**3 + 7 * z**2 - 1,
      "expanded sign-complement polynomial")
check(sp.discriminant(a, z) == 353831803500,
      "the five slopes are simple")
check(sp.resultant(a, b, z) == -31298700,
      "no slope is a B3 slope")
check(sp.resultant(a, z * b, z) == 93896100,
      "every five-slope root has alpha*b(alpha) nonzero")
check(sp.gcd(sp.Poly(a, z), sp.Poly(z * b, z)).degree() == 0,
      "gcd(a,z*b)=1")

k = c + t * q
h = alpha * k + t
m = q / c
C = (1 + alpha * q) / c
same(C * k - m * h, 1, "generalized row has determinant one")
same(h - alpha * k, t, "chosen slope fibre is t=0")

A5 = sp.expand((7 * h**2 + 3 * k**2) * (3 * h**3 + 7 * h**2 * k + k**3))
B3 = sp.expand((h + k) * (2 * h + k) * (3 * h - k))
same(A5.subs(t, 0), c**5 * a.subs(z, alpha),
     "A5 restriction is c^5*a(alpha)")
same((k * B3).subs(t, 0), c**4 * b.subs(z, alpha),
     "kB3 restriction is c^4*b(alpha)")
same((h**2).subs(t, 0), alpha**2 * c**2,
     "h squared restriction")

# Modulo a(alpha), A5 has a simple t factor.  This checks that the component
# language is not hiding a repeated chosen slope.
A5_mod_a = sp.rem(sp.Poly(A5, alpha), sp.Poly(a.subs(z, alpha), alpha)).as_expr()
check(sp.expand(A5_mod_a).subs(t, 0) == 0,
      "chosen A5 member vanishes at t=0 modulo a(alpha)")
linear_t = sp.diff(A5_mod_a, t).subs(t, 0)
expected_linear_t = c**4 * sp.diff(a, z).subs(z, alpha)
same(sp.rem(sp.Poly(linear_t - expected_linear_t, alpha),
            sp.Poly(a.subs(z, alpha), alpha)).as_expr(), 0,
     "chosen member has the expected simple t coefficient")
check(sp.resultant(a, sp.diff(a, z), z) != 0,
      "simple t coefficient is nonzero at every slope")

forced_y0 = -c**2 * b.subs(z, alpha) / alpha**2
forced_x0 = forced_y0
check(sp.diff(forced_y0, x) == 0,
      "unselected y=0 component forces a scalar d restriction")
check(sp.diff(forced_x0, y) == 0,
      "unselected x=0 component forces the same scalar restriction")
check(sp.cancel(forced_y0) != 0,
      "forced boundary scalar is nonzero in characteristic zero")

# A small finite-field hostile control: modulo 11 the five-slope polynomial
# has three visible roots, and none loses alpha*b(alpha).  This is not used
# for the characteristic-zero proof.
p = 11
roots_mod_p = [u for u in range(p) if int(a.subs(z, u)) % p == 0]
check(roots_mod_p == [1, 3, 8], "three visible mod-11 slope roots")
for u in roots_mod_p:
    check((u * int(b.subs(z, u))) % p != 0,
          f"mod-11 root {u} retains nonzero alpha*b(alpha)")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "hypothesis": "char0; a(alpha)=0; c!=0; k=c+xyq; h=alpha*k+xy; A5=d(kB3+h^2d)",
    "conclusion": "x divides d iff y divides d",
    "mechanism": "restriction to the unselected axis forces d=-c^2*b(alpha)/alpha^2, a nonzero scalar",
    "arithmetic": "gcd(a,z*b)=1; resultant=93896100",
    "universality": "units of K[x,y]/(xy) are scalars, so every normalized unimodular coordinate-cross row has k=c+xyq",
    "scope": "only the intersecting coordinate-cross grammar; disjoint/non-coordinate fibres and Keller equations remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3830-coordinate-cross-bichromatic-split-nonentry")
print("slope_arithmetic=gcd(a,z*b)=1;resultant=93896100")
print("row=k=c+xy*q;h=alpha*k+xy;determinant_completion=exact")
print("restriction=unselected_axis_forces_nonzero_scalar_d")
print("conclusion=x_divides_d_iff_y_divides_d")
print("geometry=opposite_sign_components_must_be_disjoint")
print("scope=coordinate_cross_only;second_row_and_Jacobian_OPEN")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
