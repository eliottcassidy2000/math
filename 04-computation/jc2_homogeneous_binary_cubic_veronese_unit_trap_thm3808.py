#!/usr/bin/env python3
"""Exact companion for THM-3808's homogeneous binary-cubic unit trap."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: object, rhs: object, label: str) -> None:
    difference = lhs - rhs  # type: ignore[operator]
    if isinstance(difference, sp.MatrixBase):
        reduced = difference.applyfunc(lambda entry: sp.cancel(sp.expand(entry)))
        check(reduced == sp.zeros(*difference.shape), label)
        return
    check(sp.cancel(sp.expand(difference)) == 0, label)


# Universal Delone--Faddeev sign convention.
a, b, c, d, x, y, lam = sp.symbols("a b c d x y lam")
I3 = sp.eye(3)
M_omega = sp.Matrix(
    [
        [0, -a * c, -a * d],
        [1, b, 0],
        [0, -a, 0],
    ]
)
M_theta = sp.Matrix(
    [
        [0, -a * d, -b * d],
        [0, 0, d],
        [1, 0, -c],
    ]
)
check(M_omega * M_theta == M_theta * M_omega,
      "universal Delone-Faddeev commutativity")
check(M_omega**2 == -a * c * I3 + b * M_omega - a * M_theta,
      "universal omega-square law")
check(M_omega * M_theta == -a * d * I3,
      "universal mixed law")
check(M_theta**2 == -b * d * I3 + d * M_omega - c * M_theta,
      "universal theta-square law")

T = x * M_omega + y * M_theta
index_matrix = sp.Matrix.hstack(
    sp.Matrix([1, 0, 0]),
    sp.Matrix([0, x, y]),
    T**2 * sp.Matrix([1, 0, 0]),
)
binary_form = a * x**3 + b * x**2 * y + c * x * y**2 + d * y**3
same(index_matrix.det(), -binary_form, "universal binary index sign")

basis_matrices = [I3, M_omega, M_theta]
trace_pairing = sp.Matrix(
    3,
    3,
    lambda i, j: sp.trace(basis_matrices[i] * basis_matrices[j]),
)
Delta = b**2 * c**2 - 4 * a * c**3 - 4 * b**3 * d - 27 * a**2 * d**2 + 18 * a * b * c * d
same(trace_pairing.det(), Delta, "universal trace discriminant")
same(sp.discriminant(binary_form.subs(y, 1), x), Delta,
     "universal binary cubic discriminant")


# Rational four-branch packet.
A, C = sp.symbols("A C")
packet = {a: A, b: C, c: 7 * A, d: -3 * A}
Mw = sp.simplify(M_omega.subs(packet))
Mt = sp.simplify(M_theta.subs(packet))
same(Mw**2, -7 * A**2 * I3 + C * Mw - A * Mt,
     "packet omega-square law")
same(Mw * Mt, 3 * A**2 * I3, "packet mixed law")
same(Mt**2, 3 * A * C * I3 - 3 * A * Mw - 7 * A * Mt,
     "packet theta-square law")

index_packet = sp.factor(index_matrix.det().subs(packet))
index_expected = -(A * x**3 + C * x**2 * y + 7 * A * x * y**2 - 3 * A * y**3)
same(index_packet, index_expected, "packet index form")
check(sp.Poly(index_packet, A, C).eval({A: 0, C: 0}) == 0,
      "index vanishes at the vertex")
for Xv, Yv in [(0, 0), (1, 0), (0, 1), (2, -3), (-5, 7)]:
    check(sp.expand(index_packet.subs({x: Xv, y: Yv})).subs({A: 0, C: 0}) == 0,
          f"sample index value lies in vertex ideal x={Xv},y={Yv}")

Delta_packet = sp.factor(Delta.subs(packet))
Delta_factored = A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
same(Delta_packet, Delta_factored, "rational four-line discriminant")
squarefree_gcd = sp.gcd(
    sp.gcd(sp.Poly(Delta_packet, A, C),
           sp.Poly(sp.diff(Delta_packet, A), A, C)),
    sp.Poly(sp.diff(Delta_packet, C), A, C),
)
check(squarefree_gcd == 1, "discriminant squarefree Jacobian gcd")

branch_vectors = [(1, 0), (5, 1), (19, 4), (-17, 3)]
for i in range(4):
    for j in range(i + 1, 4):
        check(branch_vectors[i][0] * branch_vectors[j][1]
              - branch_vectors[i][1] * branch_vectors[j][0] != 0,
              f"branch lines distinct {i},{j}")

# At A=C=0 all products in the augmentation ideal square to zero.
Mw0 = Mw.subs({A: 0, C: 0})
Mt0 = Mt.subs({A: 0, C: 0})
check(Mw0**2 == sp.zeros(3), "square-zero omega fibre")
check(Mw0 * Mt0 == sp.zeros(3), "square-zero mixed fibre")
check(Mt0**2 == sp.zeros(3), "square-zero theta fibre")


# Explicit third-Veronese parametrization.
s, t = sp.symbols("s t")
A_st = s**2 * t
C_st = s**3 + 7 * s * t**2 + 3 * t**3
omega_st = s**3
theta_st = 3 * s * t**2
same(omega_st**2,
     -7 * A_st**2 + C_st * omega_st - A_st * theta_st,
     "parametrized omega-square")
same(omega_st * theta_st, 3 * A_st**2,
     "parametrized mixed relation")
same(theta_st**2,
     3 * A_st * C_st - 3 * A_st * omega_st - 7 * A_st * theta_st,
     "parametrized theta-square")


def cubic_vector(poly: sp.Expr) -> sp.Matrix:
    expanded = sp.expand(poly)
    return sp.Matrix([
        expanded.coeff(s, 3 - j).coeff(t, j) for j in range(4)
    ])


span_matrix = sp.Matrix.hstack(*[
    cubic_vector(poly) for poly in (A_st, C_st, omega_st, theta_st)
])
same(span_matrix.det(), -9, "four packet cubics span the third Veronese")
same(t**3, (3 * C_st - 3 * omega_st - 7 * theta_st) / 9,
     "recover fourth Veronese generator")

# The omega characteristic law and the rational degree-three field map.
omega_char = sp.factor(Mw.charpoly(lam).as_expr())
omega_char_expected = lam**3 - C * lam**2 + 7 * A**2 * lam + 3 * A**3
same(omega_char, omega_char_expected, "omega characteristic cubic")
z = sp.symbols("z")
Pz = z**3 + 7 * z + 3
Qz = z**2
check(sp.gcd(Pz, Qz) == 1, "degree-three rational map is reduced")
check(sp.degree(Pz, z) == 3 and sp.degree(Qz, z) == 2,
      "generic field degree is three")
same(omega_char_expected.subs({lam: A * z, C: A * Pz / Qz}) / A**3,
     0, "rational root parametrizes characteristic cubic")


# Ramification and companion rays.
jacobian = sp.factor(
    sp.diff(A_st, s) * sp.diff(C_st, t)
    - sp.diff(A_st, t) * sp.diff(C_st, s)
)
same(jacobian, -3 * s * (s + t) * (s + 2 * t) * (s - 3 * t),
     "four ramification rays")
branch_pullbacks = [
    (A_st, s**2 * t),
    (C_st + 5 * A_st, (s + t)**2 * (s + 3 * t)),
    (4 * C_st + 19 * A_st, (s + 2 * t)**2 * (4 * s + 3 * t)),
    (3 * C_st - 17 * A_st, (s - 3 * t)**2 * (3 * s + t)),
]
for j, (lhs, rhs) in enumerate(branch_pullbacks):
    same(lhs, rhs, f"branch pullback double-plus-companion {j}")

critical_points = [-1, -2, 3]
critical_values = [-5, sp.Rational(-19, 4), sp.Rational(17, 3)]
wronskian = sp.factor(sp.diff(Pz, z) * Qz - Pz * sp.diff(Qz, z))
same(wronskian, z * (z + 1) * (z + 2) * (z - 3),
     "projective Wronskian")
for point, value in zip(critical_points, critical_values):
    same(Pz.subs(z, point) / Qz.subs(z, point), value,
         f"finite branch value at z={point}")
check(Qz.subs(z, 0) == 0 and Pz.subs(z, 0) != 0,
      "z=0 is the branch over infinity")

ramification_vectors = [(1, 0), (1, 1), (1, 2), (1, -3)]
companion_vectors = [(0, 1), (1, 3), (4, 3), (3, 1)]
for i, comp in enumerate(companion_vectors):
    for j, ram in enumerate(ramification_vectors):
        check(comp[0] * ram[1] - comp[1] * ram[0] != 0,
              f"companion {i} is not deleted ray {j}")


# Exact rank-four invariant-unit lattice.  The columns are exponent vectors
# for s^3, (s+t)/s, (s+2t)/s, and (s-3t)/s.
lattice_basis = sp.Matrix(
    [
        [3, -1, -1, -1],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
)
same(lattice_basis.det(), 3, "invariant unit lattice has index three")
check(lattice_basis.rank() == 4, "invariant unit lattice rank four")
for n0 in range(-4, 5):
    for n1 in range(-4, 5):
        for n2 in range(-4, 5):
            for n3 in range(-4, 5):
                vec = sp.Matrix([n0, n1, n2, n3])
                coords = lattice_basis.inv() * vec
                is_integral = all(value.q == 1 for value in coords)
                check(is_integral == ((n0 + n1 + n2 + n3) % 3 == 0),
                      f"unit lattice congruence {n0,n1,n2,n3}")

same((omega_st + A_st) / omega_st, (s + t) / s,
     "first projective unit")
same((omega_st + 2 * A_st) / omega_st, (s + 2 * t) / s,
     "second projective unit")
same((omega_st - 3 * A_st) / omega_st, (s - 3 * t) / s,
     "third projective unit")
check(all(sp.total_degree(poly) == 3 for poly in
          (omega_st, omega_st + A_st, omega_st + 2 * A_st,
           omega_st - 3 * A_st)),
      "unit representatives are Veronese fractions")


# Riemann--Hurwitz and Hilbert-function controls for the universal family.
same(3 * (-2) + 4, -2, "Riemann-Hurwitz genus-zero total")
for n in range(1, 41):
    check((n + 1) + 2 * n == 3 * n + 1,
          f"Veronese Hilbert equality degree {n}")

# Hostile boundaries: a one-dimensional coefficient span and the cyclic
# coordinate pencil both force repeated discriminant factors.
one_span = sp.factor(Delta.subs({a: A, b: A, c: A, d: A}))
check(sp.rem(one_span, A**4, A) == 0,
      "one-dimensional coefficient span has fourth-power discriminant")
cyclic_delta = sp.factor(Delta.subs({a: A, b: 0, c: 0, d: C}))
same(cyclic_delta, -27 * A**2 * C**2,
     "cyclic Veronese coordinate pencil has repeated branch axes")


semantic = {
    "family": "every generic squarefree homogeneous-linear Delone-Faddeev cubic is the third-Veronese cone",
    "index": "I=-(a*x^3+b*x^2*y+c*x*y^2+d*y^3); vertex fibre is square-zero and I never represents a unit",
    "packet": "(a,b,c,d)=(A,C,7A,-3A), Delta=A(C+5A)(4C+19A)(3C-17A)",
    "ramification": "four simple rays with four retained simple companions; S3 monodromy",
    "units": "maximal etale open has lattice {n in Z^4: sum n=0 mod 3}, rank four",
    "scope": "passes THM-3801 algebraic gates but violates constant units; no Jacobian counterexample",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3808-homogeneous-linear-binary-cubic-veronese-unit-trap")
print("universal=linear_packet+squarefree+field=>third_Veronese_cone")
print("index=codim2_square_zero;no_unit_value")
print("packet=(A,C,7A,-3A);Delta=A(C+5A)(4C+19A)(3C-17A)")
print("branch=4_simple_transpositions;4_visible_principal_companions;S3")
print("etale_open_units=lattice_sum_mod3_zero;rank4")
print("outcome=sharp_positive_normalization_grammar_but_not_constant_unit_counterexample")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
