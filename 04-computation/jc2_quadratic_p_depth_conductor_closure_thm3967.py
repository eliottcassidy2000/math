#!/usr/bin/env python3
"""Exact companion for THM-3967 (quadratic-P-depth conductor closure)."""

import hashlib
import json
import sympy as sp


CHECKS = 0


def gate(condition, label):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise AssertionError(label)


def zero(expr, label):
    gate(sp.expand(expr) == 0, label)


h, P, T = sp.symbols("h P T")
a, b, d, r, s, p = sp.symbols("a b d r s p", nonzero=True)

# Hidden polynomial and forced h^2 row.
q = a * P**2 + b * P + d
K = sp.expand(q.subs(P, h**2) - 2 * h**3)
zero(K - (a * h**4 - 2 * h**3 + b * h**2 + d), "hidden quartic")
zero(K.subs(d, 0) - h**2 * (a * h**2 - 2 * h + b), "forced h2 factor")

# A quadratic repeated factor would square-fill the quartic. Its h^3 row
# fixes s=-1/a and its absent h row then forces p=0.
M = h**2 + s * h + p
M2 = sp.Poly(sp.expand(a * M**2), h)
gate(M2.coeff_monomial(h**3) == 2 * a * s, "quadratic square h3 row")
gate(M2.coeff_monomial(h) == 2 * a * s * p, "quadratic square h1 row")
zero((2 * a * s).subs(s, -1 / a) + 2, "h3 fixes s")
zero(M.subs({s: -1 / a, p: 0}) - h * (h - 1 / a),
     "quadratic factor becomes reducible")

# Repeated-root equations and exact graph normal form.
b_root = 3 * r - 2 * a * r**2
d_root = a * r**4 - r**3
zero(sp.diff(K, h).subs({h: r, b: b_root}), "repeated root derivative")
zero(K.subs({h: r, b: b_root, d: d_root}), "repeated root value")
q_graph = 3 * r * P - r**3 + a * (P - r**2) ** 2
zero(q_graph - q.subs({b: b_root, d: d_root}), "graph coefficient identity")
residual = a * h**2 + (2 * a * r - 2) * h + r * (a * r - 1)
zero(K.subs({b: b_root, d: d_root}) - (h - r) ** 2 * residual,
     "full repeated-root factorization")
zero(sp.discriminant(residual, h) - 4 * (1 - a * r),
     "residual discriminant")
zero(residual.subs(h, r) - r * (4 * a * r - 3),
     "triple-root seam")

# The full hidden discriminant freezes the coefficient hypersurface.
Delta = (16 * a**3 * d**2 - 8 * a**2 * b**2 * d + a * b**4
         + 36 * a * b * d - b**3 - 27 * d)
zero(sp.discriminant(K, h) - 16 * d * Delta,
     "hidden quartic discriminant")

# The a=0 endpoint is reducible.
q_linear = 3 * r * P - r**3
F_linear = T**3 - 3 * P * T - q_linear
zero(F_linear - (T + r) * (T**2 - r * T + r**2 - 3 * P),
     "a0 factorization")

# Exact cleared-denominator rows for r=S/R. The proof uses gcd(R,S)=1 to
# infer R|a, then the two reductions 2*a1*S=3 and a1*S=1 modulo R.
R, S, a1 = sp.symbols("R S a1", nonzero=True)
den_b = sp.expand((b_root - b) * R**2).subs(r, S / R)
den_d = sp.expand((d_root - d) * R**4).subs(r, S / R)
zero(den_b - (3 * S * R - 2 * a * S**2 - b * R**2),
     "cleared b denominator row")
zero(den_d - (a * S**4 - S**3 * R - d * R**4),
     "cleared d denominator row")
zero((2 * a1 * S - 3) - 2 * (a1 * S - 1) + 1,
     "incompatible denominator residues leave one")

# In the forced h^2 row, repetition is exactly ab=1 and is constant over
# k[t].
L_forced = a * h**2 - 2 * h + b
gate(sp.discriminant(L_forced, h) == 4 * (1 - a * b),
     "forced-row repetition iff ab1")

# A positive graph control.
aa, rr = sp.symbols("aa rr", nonzero=True)
bb = 3 * rr - 2 * aa * rr**2
dd = aa * rr**4 - rr**3
Kpos = aa * h**4 - 2 * h**3 + bb * h**2 + dd
gate(sp.rem(sp.Poly(Kpos, h), sp.Poly(h - rr, h)) == 0,
     "graph control root")
gate(sp.rem(sp.Poly(sp.diff(Kpos, h), h), sp.Poly(h - rr, h)) == 0,
     "graph control repeated")

# Depth-three hostile: a genuinely quadratic repeated hidden factor, but a
# polynomial T-root makes the cubic reducible.
t = sp.symbols("t")
s3 = P + t
q3 = 3 * s3 * P - s3**3
K3 = sp.expand(q3.subs(P, h**2) - 2 * h**3)
s3h = h**2 + t
zero(K3 + (h - s3h) ** 2 * (2 * h + s3h), "depth3 nongraph square")
M3 = h**2 - h + t
gate(sp.discriminant(M3, h) == 1 - 4 * t, "depth3 factor discriminant")
F3 = sp.expand(T**3 - 3 * P * T - q3)
zero(F3 - (T + s3) * (T**2 - s3 * T + s3**2 - 3 * P),
     "depth3 domain failure")

summary = {
    "checks": CHECKS,
    "universe": "irreducible T^3-3PT-q with deg_P(q)<=2",
    "hidden": "K=a*h^4-2*h^3+b*h^2+d",
    "quadratic_repeat": "square-fill forces p=0 and reducibility",
    "linear_repeat": "b=3r-2ar^2; d=ar^4-r^3",
    "denominator": "coprime denominator congruences force r in k[t]",
    "classification": "squarefree normal / moving P2 / THM3964 graph",
    "conclusion": "quadratic coefficient depth has no Keller chart",
    "sharpness": "depth3 nongraph repeated factor exists but cubic splits",
    "scope": "depth>=3 affine-P graph multipliers and JC2 remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3967 quadratic-P-depth conductor closure companion")
print(f"CHECKS={CHECKS}")
print("UNIVERSE=IRREDUCIBLE_NATURAL_CUBIC;DEG_P_Q_LE_2")
print("HIDDEN=AH4_MINUS_2H3_PLUS_BH2_PLUS_D")
print("QUADRATIC_REPEAT=SQUARE_FILL_FORCES_REDUCIBLE_M")
print("LINEAR_REPEAT=EXACT_POLYNOMIAL_GRAPH_AFTER_UFD_DENOMINATOR_GATE")
print("TRICHOTOMY=SQUAREFREE_NORMAL;MOVING_P2;THM3964_GRAPH")
print("CONCLUSION=NO_SAME_FIELD_PLANAR_KELLER_CHART")
print("SHARP_DEPTH3=NONGRAPH_REPEAT_EXISTS_BUT_DOMAIN_SPLITS")
print("SCOPE=AFFINE_P_MULTIPLIER;HIGHER_DEPTH;NONMONOGENIC;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
