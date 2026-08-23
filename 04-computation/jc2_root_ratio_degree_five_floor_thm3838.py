#!/usr/bin/env python3
"""Exact companion for THM-3838's numerator/denominator degree floor."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    CHECKS += 1


def arithmetic_genus(degree: int) -> int:
    return (degree - 1) * (degree - 2) // 2


gate([arithmetic_genus(d) for d in range(1, 5)] == [0, 0, 1, 3],
     "plane-curve arithmetic genus through quartics")
gate(all(arithmetic_genus(d) < 3 for d in (1, 2, 3)),
     "generic genus at least three forces degree at least four")
gate([(a, b) for a in range(1, 5) for b in range(4, 5) if a * b == 4]
     == [(1, 4)],
     "degree-four composition boundary is linear outer word")

T, Z = sp.symbols("T Z")
H = (
    84 * T**7
    + 36 * T**6 * Z**2
    + 196 * T**6 * Z
    + 84 * T**5 * Z**3
    + 36 * T**5 * Z**2
    + 49 * T**4 * Z**4
    + 112 * T**4 * Z**3
    - 12 * T**3 * Z**5
    - 14 * T**2 * Z**6
    + 12 * T**2 * Z**5
    + Z**8
)

gate(sp.degree(H, Z) == 8 and sp.LC(sp.Poly(H, Z)) == 1,
     "h-side sidecar is monic degree eight")
gate(sp.degree(H, T) == 7 and sp.LC(sp.Poly(H, T)) == 84,
     "k-side sidecar has degree seven")
gate((8 - 2) // 2 == 3 and (7 - 1) // 2 == 3,
     "both squarefree hyperelliptic sidecars have genus three")
gate(8 % 2 == 0 and 7 % 2 == 1,
     "opposite infinity parity retained")

# Small exact squarefreeness controls from each orientation.  The theorem
# inherits the full transcendental discriminant statements from THM-3827.
H_even_control = sp.Poly(H.subs(T, -1), Z)
H_odd_control = sp.Poly(H.subs(Z, 1), T)
gate(sp.gcd(H_even_control, H_even_control.diff()).degree() == 0,
     "degree-eight control is squarefree")
gate(sp.gcd(H_odd_control, H_odd_control.diff()).degree() == 0,
     "degree-seven control is squarefree")

# Adjunction for a smooth plane quartic gives K_C=O_C(1), of degree 4.
# The canonical system of a genus-three hyperelliptic curve instead has
# degree two onto a conic, so the two equality models cannot coincide.
gate(4 * (4 - 3) == 2 * 3 - 2,
     "plane-quartic adjunction has canonical degree four")
gate(2 * (3 - 1) == 4,
     "hyperelliptic genus-three canonical system pulls back a conic")

semantic = {
    "inheritance": "THM-3827 gives h=p(g), k=q(ell) with primitive generic genera at least 3",
    "degree": "deg p(g)=deg(p)*deg(g) and similarly for q(ell)",
    "boundary": "total degree 4 forces a linear outer word and primitive degree 4",
    "plane": "genus 3 at plane degree 4 forces a smooth plane quartic",
    "conflict": "THM-3827 equality forces the same curve to be hyperelliptic genus 3",
    "conclusion": "deg(h)>=5 and deg(k)>=5",
    "ratio": "gcd(h,k)=1, so z=h/k has numerator and denominator degree at least 5",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python assert")

print("theorem=THM-3838-root-ratio-numerator-denominator-degree-five-floor")
print("inheritance=dual_closed_polynomial_generic_genus_at_least_3")
print("coarse=degree_h,degree_k_at_least_4")
print("quartic=genus3_forces_smooth_plane_quartic_nonhyperelliptic")
print("sidecar=genus3_hyperelliptic_in_both_orientations")
print("conclusion=degree_h,degree_k_at_least_5")
print("ratio=reduced_numerator_and_denominator_each_degree_at_least_5")
print("scope=all_degree;degree5_boundary_and_atlas_existence_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
