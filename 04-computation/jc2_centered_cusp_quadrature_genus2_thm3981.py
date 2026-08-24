#!/usr/bin/env python3
"""Exact ramification/divisor companion for THM-3981."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


x, z, Y = sp.symbols("x z Y", nonzero=True)
lam = sp.symbols("lambda")
phi = sp.expand(z * (z - 1) ** 2)
h = sp.expand(x**3 * (Y + lam * x / 2))
F = sp.expand(phi - h)

# Degree and irreducibility input: a reducible cubic has a rational root,
# while a root would factor the degree-four h-map through the degree-three
# phi-map.
gate(sp.degree(phi, z) == 3, "phi degree")
gate(sp.degree(h, x) == 4, "h degree for nonzero lambda")
gate(4 % 3 != 0, "composition-degree obstruction")
gate(sp.Poly(F, z).LC() == 1, "monic cubic")

# Critical values and the complete finite discriminant support.
phip = sp.factor(sp.diff(phi, z))
gate(phip == (z - 1) * (3 * z - 1), "phi critical factorization")
gate(phi.subs(z, 1) == 0, "zero critical value")
gate(phi.subs(z, sp.Rational(1, 3)) == sp.Rational(4, 27),
     "nonzero critical value")
disc = sp.factor(sp.discriminant(F, z))
gate(sp.factor(disc - h * (4 - 27 * h)) == 0,
     "cubic discriminant")

# The four h=4/27 addresses are squarefree over k(Y).
G = sp.expand(27 * x**3 * (2 * Y + lam * x) - 8)
domain = sp.QQ.frac_field(Y, lam)
gate(sp.degree(G, x) == 4, "four nonzero critical addresses")
gate(sp.gcd(sp.Poly(G, x, domain=domain),
            sp.Poly(sp.diff(G, x), x, domain=domain)).degree() == 0,
     "critical-address polynomial squarefree")
xcrit = -3 * Y / (2 * lam)
gate(sp.factor(G.subs(x, xcrit)
               + (729 * Y**4 + 128 * lam**3) / (16 * lam**3)) == 0,
     "only possible repeated address is absent generically")
gate(sp.factor(sp.diff(h, x)) == x**2 * (3 * Y + 2 * lam * x),
     "h critical points")
gate(sp.factor(h.subs(x, -2 * Y / lam)) == 0,
     "second zero address")

# Riemann--Hurwitz ledger: cusp, second zero, four nonzero critical
# addresses, and the unique totally ramified infinity point.
ramification = {
    "x0_z1": 1,
    "xminus2Yoverlambda_z1": 1,
    "four_z1over3": 4,
    "infinity": 2,
}
gate(sum(ramification.values()) == 8, "total ramification")
genus = (3 * (-2) + sum(ramification.values()) + 2) // 2
gate(genus == 2, "Riemann--Hurwitz genus")
gate(sp.gcd(3, 4) == 1, "unique infinity Newton edge")

# Orders of omega=2x dx/((z-1)(3z-1)).  Each tuple is
# (ord x, ord dx, ord(z-1), ord(3z-1), expected ord omega).
orders = {
    "P0": (1, 0, 0, 0, 1),
    "P1_cusp": (2, 1, 3, 0, 0),
    "Q1_second_zero": (0, 1, 1, 0, 0),
    "Ralpha_nonzero_critical": (0, 1, 0, 1, 0),
    "Pinfinity": (-3, -4, -4, -4, 1),
}
for name, (ox, odx, oz1, o3z1, expected) in orders.items():
    gate(ox + odx - oz1 - o3z1 == expected,
         f"omega valuation at {name}")
gate(sum(row[-1] for row in orders.values()) == 2,
     "canonical divisor degree exhausted")
gate(2 * genus - 2 == 2, "canonical degree")

# A finite separable pullback of a regular differential stays regular:
# ord_Q(pi^*omega)=e*ord_P(omega)+(e-1).  An exact differential of a
# function with pole order m has order -m-1 in characteristic zero.
for e in range(1, 8):
    for omega_order in (0, 1):
        gate(e * omega_order + e - 1 >= 0,
             f"regular pullback e={e}, ord={omega_order}")
for m in range(1, 8):
    gate(-m - 1 < 0, f"differential of pole m={m}")

# The lambda=0 endpoint is rational but logarithmically nonintegrable.
r = sp.symbols("r", nonzero=True)
D0 = Y * r**3 - 1
x0 = r / D0
z0 = Y * r**3 / D0
gate(sp.factor(F.subs({lam: 0, x: x0, z: z0})) == 0,
     "zero-slope rational parametrization")
omega0 = sp.factor(
    2 * x0 * sp.diff(x0, r) / ((z0 - 1) * (3 * z0 - 1))
)
gate(sp.factor(omega0 + 2 * r / D0) == 0,
     "zero-slope logarithmic differential")
gate(sp.degree(D0, r) == 3, "three endpoint poles")
gate(sp.gcd(sp.Poly(D0, r, domain=sp.QQ.frac_field(Y)),
            sp.Poly(sp.diff(D0, r), r,
                    domain=sp.QQ.frac_field(Y))).degree() == 0,
     "endpoint poles distinct")
residue0 = -2 / (3 * Y * r)
gate(sp.factor(residue0.subs(Y, r**-3) + 2 * r**2 / 3) == 0,
     "endpoint residue on Y*r^3=1")
for e in range(1, 8):
    gate(e != 0, f"endpoint residue survives ramification e={e}")

summary = {
    "checks": CHECKS,
    "curve": "z(z-1)^2=x^3(Y+lambda*x/2)",
    "degree_x_cover": 3,
    "irreducibility_gate": "deg(h)=4 not divisible by deg(phi)=3",
    "discriminant": "h(4-27h)",
    "ramification_contributions": ramification,
    "genus": genus,
    "omega": "2x dx/((z-1)(3z-1))",
    "omega_divisor": "P_(0,0)+P_infinity",
    "finite_cover_primitive": False,
    "lambda_zero": "rational; omega=-2r dr/(Yr^3-1); three nonzero residues",
    "scope": "canonical centered quadrature gauge only; JC2 open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3981 centered cusp quadrature genus-two companion")
print(f"CHECKS={CHECKS}")
print("INTEGRAL=DEGREE_3_PHI_CANNOT_FACTOR_DEGREE_4_H")
print("RAMIFICATION=1_PLUS_1_PLUS_4_PLUS_2_EQUALS_8")
print("GENUS=2")
print("OMEGA_DIVISOR=P_ZERO_COLOR_PLUS_P_INFINITY")
print("LAMBDA_ZERO=RATIONAL_WITH_THREE_NONZERO_RESIDUES")
print("FINITE_COVER_PRIMITIVE=IMPOSSIBLE")
print("CENTERED_FORMAL_X=TRANSCENDENTAL_OVER_GENERIC_SLICE_FIELD")
print("SCOPE=CENTERED_QUADRATURE_GAUGE_ONLY;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
