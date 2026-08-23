#!/usr/bin/env python3
"""Exact companion for THM-3827's generic-fibre genus floor."""

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


T, Z, s, x = sp.symbols("T Z s x")
H = (
    84 * T**7 + 36 * T**6 * Z**2 + 196 * T**6 * Z
    + 84 * T**5 * Z**3 + 36 * T**5 * Z**2 + 49 * T**4 * Z**4
    + 112 * T**4 * Z**3 - 12 * T**3 * Z**5 - 14 * T**2 * Z**6
    + 12 * T**2 * Z**5 + Z**8
)

H_poly = sp.Poly(H, Z)
check(H_poly.degree() == 8, "sidecar has degree eight")
check(H_poly.LC() == 1, "sidecar is monic")

P = Z**8 - 14 * Z**6 + 24 * Z**5 + 49 * Z**4 + 28 * Z**3 + 196 * Z - 84
same(H.subs(T, -1), P, "T=-1 hostile fibre")
P_poly = sp.Poly(P, Z)
check(sp.gcd(P_poly, P_poly.diff()).degree() == 0,
      "hostile fibre is squarefree")

Delta = sp.factor(sp.discriminant(H, Z))
Delta_factor = (
    19144454963200 * T**44
    * (
        466560000 * T**8 - 96264205056 * T**7
        - 180219652272 * T**6 - 35458022604096 * T**5
        - 85186647741177 * T**4 - 616061164068 * T**3
        - 497726825264 * T**2 + 2712579840 * T - 678144960
    )
)
same(Delta, Delta_factor, "exact generic discriminant factorization")
check(Delta != 0, "generic sidecar is squarefree")
same(Delta.subs(T, -1), sp.discriminant(P, Z),
     "discriminant specializes at the hostile fibre")

# A nonconstant composition cannot annihilate a nonzero one-variable
# polynomial at a transcendental element.  This explicit control is not the
# human proof, but catches accidental use of T rather than p(s).
p_control = s**3 - 2 * s + 5
check(sp.Poly(Delta.subs(T, p_control), s).degree() > 0,
      "nonconstant Stein composition retains nonzero discriminant")

# A monic squarefree polynomial of degree 2g+2 has hyperelliptic genus g.
genus_target = (H_poly.degree() - 2) // 2
check(genus_target == 3, "sidecar genus is three")
for genus_source in (0, 1, 2):
    # Riemann--Hurwitz would require 2g_X-2 >= 2g_Y-2 for degree >=1.
    check(2 * genus_source - 2 < 2 * genus_target - 2,
          f"Riemann--Hurwitz excludes source genus {genus_source}")

# A second completed square turns the atlas pencil into a five-slope
# reducibility gate.
h, k = sp.symbols("h k")
H_hk = H.subs({T: h, Z: k})
A_pencil = (7 * h**2 + 3 * k**2) * (3 * h**3 + 7 * h**2 * k + k**3)
B_pencil = (h + k) * (2 * h + k) * (3 * h - k)
same(H_hk, (k * B_pencil) ** 2 + 4 * h**2 * A_pencil,
     "five-slope completed square")

z = sp.symbols("z")
A_affine = sp.expand(A_pencil.subs({h: z, k: 1}))
check(sp.degree(A_affine, z) == 5, "pencil packet has five finite slopes")
check(sp.discriminant(A_affine, z) == 353831803500,
      "five pencil slopes are distinct")
check(A_pencil.subs({h: 1, k: 0}) == 21,
      "pencil packet has no root at the k-zero point")
check(sp.resultant(A_affine, sp.expand((k * B_pencil).subs({h: z, k: 1})), z)
      == -31298700,
      "size-one subset cannot have the kB complement")
check(sp.gcd(sp.Poly(h**2, h, k), sp.Poly(k * B_pencil, h, k)).total_degree() == 0,
      "size-two cancellation cannot make h squared divide kB")

# If d has homogeneous degree r in A=d(kB+h^2d), all r except 1 and 2
# have an unmatched extreme degree before coefficient arithmetic begins.
possible_subset_sizes = []
for r in range(6):
    degrees = (r + 4, 2 * r + 2)
    if max(degrees) == 5 or degrees[0] == degrees[1]:
        possible_subset_sizes.append(r)
check(possible_subset_sizes == [1, 2],
      "only subset sizes one and two survive degree screening")

# Exact squarefree/non-squarefree special-fibre controls for h=p(g).
p_squarefree = (x - 1) * (x + 2)
p_repeated = (x - 1) ** 2 * (x + 2)
check(sp.gcd(sp.Poly(p_squarefree, x), sp.Poly(p_squarefree, x).diff()).degree() == 0,
      "squarefree arm control")
check(sp.gcd(sp.Poly(p_repeated, x), sp.Poly(p_repeated, x).diff()).degree() > 0,
      "repeated-root arm hostile control")


source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "hypothesis": "h=p(g), p nonconstant, K(g) relatively algebraically closed in K(x,y)",
    "sidecar": "W^2=H(p(g),Z), monic squarefree degree eight, geometric genus three",
    "mechanism": "a nonconstant k gives a curve map; Riemann-Hurwitz excludes generic-fibre genus <=2",
    "intersection": "K[x,y] intersect K(g)=K[g] by a denominator-root valuation",
    "arm": "etale pullback makes p(g) squarefree and k nonconstant on every component",
    "pencil": "one of h=0 or the five fibres h-alpha_i*k=0 must be reducible",
    "scope": "generative closed factor existence is CITED, not computational; no JC counterexample",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases")
print("sidecar=monic_squarefree_degree8_hyperelliptic_genus3")
print("hostile_fibre=T=-1;gcd(P,Pprime)=1")
print("closed_polynomial=K(g)_relatively_algebraically_closed_in_K(x,y)")
print("obstruction=generic_fibre_genus_0_1_2_impossible")
print("intersection=K[x,y]_intersect_K(g)=K[g]")
print("boundary=p(g)_squarefree;each_arm_component_requires_nonconstant_k")
print("pencil=reducible_fibre_among_h_and_five_fixed_slopes")
print("stein_factor=existence_CITED_not_computational")
print("scope=no_JC_counterexample")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
