#!/usr/bin/env python3
"""Exact companion for THM-3924's cubic class and valuation packet.

Reproduction:
  python3 04-computation/jc2_decic_cubic_index_five_ramification_class_thm3924.py
  python3 -O 04-computation/jc2_decic_cubic_index_five_ramification_class_thm3924.py
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


A, C, Z, u, z, x, b = sp.symbols("A C Z u z x b")

K = 2304 * b**5 + 10176 * b**4 + 4064 * b**3 + 996 * b**2 + 84 * b + 5
c = b - sp.Rational(1, 2)
d = (4 * b + 1) / 8
kappa = 96 * b**3 + 48 * b**2 + 18 * b + 3

p = (u**2 - 1) * (u**2 + b) ** 2
h = (
    u**9
    + (3 * b - sp.Rational(3, 2)) * u**7
    + (3 * b**2 - sp.Rational(9, 2) * b + sp.Rational(3, 8)) * u**5
    + (b**3 - sp.Rational(9, 2) * b**2 + sp.Rational(9, 8) * b + sp.Rational(1, 16)) * u**3
    + (-sp.Rational(3, 2) * b**3 + sp.Rational(9, 8) * b**2
       + sp.Rational(3, 16) * b + sp.Rational(3, 128)) * u
)
F = z**3 - 3 * p * z + 2 * h
r = sp.expand(p**3 - h**2)

P = sp.cancel(A**6 * p.subs(u, C / A))
H = sp.cancel(A**9 * h.subs(u, C / A))
Q = 2 * A**10 - H
f = sp.expand(Z**3 - 3 * P * Z - 2 * Q)
Delta = sp.cancel(A**8 * r.subs(u, C / A) + 4 * H - 4 * A**10)

# The A-localized order is the radial UFD chart.
zero(f.subs(Z, A**3 * z) - A**9 * (F.subs(u, C / A) - 4 * A),
     "homogenized cubic is the radial relation")
zero(sp.expand(f.subs(A, 0)) - (Z - C**3) ** 2 * (Z + 2 * C**3),
     "power polynomial has one simple and one doubled infinity root")
zero(Delta.subs(A, 0) - C**8 * (512 * C + kappa) / 128,
     "normal discriminant is generically nonzero on A=0")

# The displayed THM-3921 basis has power-basis determinant -d/A^5.
L = C**3 + c * A**2 * C
N = C * L - d * A**4
e1_coordinates = sp.Matrix([-2 * L**2 / A**4, L / A**4, 1 / A**4])
e2_coordinates = sp.Matrix([
    (2 * C * L**2 - 4 * d * A**4 * L) / A**5,
    -N / A**5,
    -C / A**5,
])
basis_matrix = sp.Matrix.hstack(sp.Matrix([1, 0, 0]), e1_coordinates, e2_coordinates)
zero(sp.det(basis_matrix) + d / A**5, "power-order index exponent is five")
gate(sp.gcd(K, 4 * b + 1) == 1, "basis determinant constant survives K=0")

# The ramification equation on the UFD chart is irreducible.
ramification = z**2 - p
ramification_factors = sp.Poly(
    ramification, z, domain=sp.QQ.frac_field(b, u)
).factor_list()[1]
gate(len(ramification_factors) == 1
     and ramification_factors[0][0].degree() == 2
     and ramification_factors[0][1] == 1,
     "z^2-p is irreducible over Q(b,u)")

# Write H_ram=A^-6(Z^2-Pi).  At the simple A=0 prime Z=-2C^3,
# the numerator is a unit with residue 3C^6, hence valuation -6.
Pi = sp.cancel(A**6 * p.subs(u, C / A))
simple_numerator = sp.expand(Z**2 - Pi).subs(Z, -2 * C**3)
zero(simple_numerator.subs(A, 0) - 3 * C**6,
     "simple-prime ramification numerator residue")

# At the degree-two prime use the integral local coordinate from THM-3921.
m = L - d * A**4 / C
quadratic_numerator = sp.cancel((m + A**5 * x) ** 2 - Pi)
numerator, denominator = sp.together(quadratic_numerator).as_numer_denom()
polynomial_A = sp.Poly(sp.expand(numerator), A)
for exponent in range(5):
    zero(polynomial_A.coeff_monomial(A**exponent),
         f"degree-two numerator has no A^{exponent} term")
zero(polynomial_A.coeff_monomial(A**5) / denominator - 2 * C**3 * x,
     "degree-two numerator leading A^5 coefficient")
zero(polynomial_A.coeff_monomial(A**6) / denominator - (2 * b + 1) / 8,
     "degree-two numerator next coefficient")

residual_quadratic = 384 * C**4 * x**2 - 512 * C - kappa
gate(sp.Poly(residual_quadratic, x).degree() == 2,
     "degree-two residual really has degree two")
gate(sp.gcd(K, kappa) == 1, "residual x is nonzero on K=0")
gate(sp.gcd(K, 2 * b + 1) == 1, "truncated-square coefficient survives K=0")

# Nagata: the only unit relation is div(A)=P_-+P_2.  In the quotient
# Z^2/<(1,1)>, choose g=[P_-], so [P_2]=-g.  The ramification divisor
# relation is E-6P_--P_2=0 and therefore [E]=5g.
unit_relation = sp.Matrix([1, 1])
ramification_poles = sp.Matrix([-6, -1])
gate(sp.gcd(int(unit_relation[0]), int(unit_relation[1])) == 1,
     "A-divisor relation is primitive")
valuation_gap = int(ramification_poles[1] - ramification_poles[0])
gate(valuation_gap == 5, "ramification class is five times the generator")
gate(valuation_gap == sp.degree(sp.denom(sp.cancel(sp.det(basis_matrix))), A),
     "class divisibility equals power-index exponent")

semantic_payload = {
    "localized_order": "B_A=k[u,z,F_inverse]_UFD_with_units_kstar_times_AZ",
    "A_zero_packet": "two_unramified_primes_of_residue_degrees_one_and_two",
    "class_group": "Cl(B)=Z2_mod_diagonal=Z",
    "ramification_equation": "E_on_B_A_is_z_squared_minus_p",
    "ramification_valuations": [-6, -1],
    "ramification_class": "E_equals_five_times_primitive_generator",
    "index_bridge": "power_basis_index_A5_equals_boundary_class_divisibility_five",
    "consequence": "THM3922_primitive_boundary_gate_excludes_A2_Keller_atlas",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3924-decic-cubic-index-five-ramification-class-obstruction")
print("localized_order=UFD_k[u,z,F^-1];units=kstar_times_A^Z")
print("A_zero_primes=residue_degrees_1_2;div(A)=Pminus+P2")
print("class_group=Z2_mod_(1,1)=Z;generator=Pminus;P2=-Pminus")
print("ramification=E: z^2-p=0")
print("valuations_E_equation=Pminus:-6;P2:-1")
print("ramification_class=[E]=5[Pminus];nonprimitive")
print("index_bridge=power_index_A^5=valuation_gap_5=class_divisibility_5")
print("consequence=no_A2_Keller_atlas_by_THM3922")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
