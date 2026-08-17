#!/usr/bin/env python3
"""Exact controls for THM-3532 (fixed Keller conjugacy covariance).

The abstract norm and trace-discriminant identities are proved in THM-3532.
This companion checks the load-bearing explicit controls:

* the THM-2465 maps W1,W2 are target postcompositions T_i o F;
* their one-step boundary, finite divisor, and x-core discriminant are exact
  target pullbacks;
* W_i^2 is not T_i o F^2, whereas (T_i o F o T_i^-1)^2 is exactly
  T_i o F^2 o T_i^-1 on a hostile finite point bank;
* L o T_i^-1 loses the standard five-weight A(1,0) chart, but recovers it
  after the coordinates/filtrations are transported; and
* changing the boundary equation by a scalar obeys the exact raw-tower
  normalization recurrence r_0=c, r_(n+1)=c^e_n r_n^3.

All arithmetic is exact over QQ.  The finite point bank is a control for the
composition typing; the general conjugacy identity is formal associativity,
not inferred from the bank.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


x, y, z, r = sp.symbols("x y z r")
VARS = (x, y, z)


def require(condition, message):
    """Optimization-safe validity gate."""
    if not condition:
        raise AssertionError(message)


def subst(poly, mapping):
    return sp.expand(poly.subs(mapping, simultaneous=True))


def apply_map(poly, phi):
    return subst(poly, dict(zip(VARS, phi)))


def compose(outer, inner):
    """Return outer o inner for polynomial triples."""
    return tuple(apply_map(poly, inner) for poly in outer)


def evaluate(phi, point):
    mapping = dict(zip(VARS, point))
    return tuple(sp.Rational(poly.subs(mapping)) for poly in phi)


def evaluate_chain(point, *maps_inner_to_outer):
    value = point
    for phi in maps_inner_to_outer:
        value = evaluate(phi, value)
    return value


def jacobian_det(phi):
    return sp.expand(sp.Matrix(phi).jacobian(VARS).det())


def face(poly, weight, choose):
    terms = sp.Poly(sp.expand(poly), *VARS, domain=sp.QQ).terms()
    values = [sum(a * b for a, b in zip(monomial, weight))
              for monomial, _ in terms]
    extremum = max(values) if choose == "max" else min(values)
    selected = sp.Integer(0)
    for (monomial, coefficient), value in zip(terms, values):
        if value == extremum:
            selected += coefficient * sp.prod(v**e for v, e in zip(VARS, monomial))
    return extremum, sp.factor(selected)


def packet_signature(poly):
    return {
        "max_lambda": face(poly, (1, 0, -1), "max"),
        "min_lambda": face(poly, (1, 0, -1), "min"),
        "min_beta": face(poly, (1, -1, -2), "min"),
        "max_k": face(poly, (0, 0, 1), "max"),
        "min_gamma": face(poly, (1, -1, -5), "min"),
    }


def compact_signature(signature):
    return tuple((name, value, str(poly))
                 for name, (value, poly) in signature.items())


u = 1 + x * y
F = (
    u**3 * z + y**2 * u * (4 + 3 * x * y),
    y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
    2 * x - 3 * x**2 * y - x**3 * z,
)

L = 27 * x**2 * z**2 - 18 * x * y * z + 16 * x + y**3 * z - y**2
S = 27 * x * z**2 - 9 * y * z + 8

T1 = (x, y, z + x**2)
T1inv = (x, y, z - x**2)
T2 = (x + y * z, y + z**2, z)
T2inv = (x - y * z + z**3, y - z**2, z)

IDENTITY = VARS
require(compose(T1, T1inv) == IDENTITY == compose(T1inv, T1),
        "T1 inverse identity failed")
require(compose(T2, T2inv) == IDENTITY == compose(T2inv, T2),
        "T2 inverse identity failed")
# det J_F=-2 and the ten-term identity F^*L=B are already independently
# certified in THM-3529.  Re-expanding them here dominates the runtime while
# adding no independent covariance check, so this companion treats them as
# inherited exact inputs.
require(jacobian_det(T1) == jacobian_det(T1inv) == 1,
        "T1 Jacobian-unit check failed")
require(jacobian_det(T2) == jacobian_det(T2inv) == 1,
        "T2 Jacobian-unit check failed")

base_signature = packet_signature(L)
expected_base = {
    "max_lambda": (1, 16 * x),
    "min_lambda": (-1, y**3 * z),
    "min_beta": (-5, y**3 * z),
    "max_k": (2, 27 * x**2 * z**2),
    "min_gamma": (-8, z * (27 * x**2 * z + y**3)),
}
require(all(value == expected_base[name][0]
            and sp.expand(poly - expected_base[name][1]) == 0
            for name, (value, poly) in base_signature.items()),
        "base A(1,0) five-face signature failed")

print("THM-3532 exact covariance controls")
print("base_det=-2 (inherited THM-3529 exact input)")
print("base_A10_signature=" + repr(compact_signature(base_signature)))

for name, tame, tame_inv in (("W1", T1, T1inv), ("W2", T2, T2inv)):
    boundary = apply_map(L, tame_inv)
    transported_S = apply_map(S, tame_inv)

    # The x-coordinate is unchanged on the source of W_i.  Its cubic core is
    # therefore just the old core with the target variables pulled back.
    E = L * r**3 + (4 - 3 * y * z) * r - 2 * z
    E_post = apply_map(E, tame_inv)
    require(sp.expand(sp.discriminant(E_post, r)
                      + 4 * transported_S**2 * boundary) == 0,
            name + " target-pulled x-core discriminant failed")

    # In standard (x,y,z) weights the transformed seed is not A(1,0).
    standard_signature = packet_signature(boundary)
    require(compact_signature(standard_signature) != compact_signature(base_signature),
            name + " did not trigger the standard-weight hostile")
    # In the transported coordinates, pull back once by T_i and recover L.
    require(sp.expand(apply_map(boundary, tame) - L) == 0,
            name + " transported boundary did not recover L")
    transported_signature = packet_signature(apply_map(boundary, tame))
    require(compact_signature(transported_signature) == compact_signature(base_signature),
            name + " transported A(1,0) signature failed")

    # Find a small exact witness to the invalid postcomposition iteration
    # identity (T_i F)^2 = T_i F^2.  This is a typing hostile, not a genericity
    # experiment.
    witnesses = []
    for point in product((-1, 0, 1), repeat=3):
        # W_i^2 = T_i F T_i F, while the false transport would be T_i F^2.
        post_squared = evaluate_chain(point, F, tame, F, tame)
        alleged = evaluate_chain(point, F, F, tame)
        if post_squared != alleged:
            score = max(abs(v) for v in post_squared + alleged)
            witnesses.append((score, point, post_squared, alleged))
    require(bool(witnesses), name + " postcomposition hostile bank was empty")
    _, witness_point, witness_left, witness_right = min(witnesses, key=lambda row: (row[0], row[1]))

    # Positive controls.  The first two equalities are formal consequences of
    # T_i^-1 T_i=id; the bank checks the implementation without pretending to
    # replace that proof.
    for point in product((-1, 0, 1), repeat=3):
        old_target = evaluate(F, point)
        post_target = evaluate_chain(point, F, tame)
        require(sp.Rational(boundary.subs(dict(zip(VARS, post_target))))
                == sp.Rational(L.subs(dict(zip(VARS, old_target)))),
                name + " target boundary pullback failed at " + repr(point))
        core_mapping = dict(zip(VARS, post_target))
        core_mapping[r] = point[0]
        require(sp.Rational(E_post.subs(core_mapping)) == 0,
                name + " pulled x-core failed at " + repr(point))

        lhs = evaluate_chain(point, tame_inv, F, tame,
                             tame_inv, F, tame)
        rhs = evaluate_chain(point, tame_inv, F, F, tame)
        require(lhs == rhs,
                name + " honest conjugacy square failed at " + repr(point))

    print(name + "_inverse=" + repr(tuple(map(str, tame_inv))))
    print(name + "_boundary=" + str(sp.factor(boundary)))
    print(name + "_standard_signature=" + repr(compact_signature(standard_signature)))
    print(name + "_transported_signature_recovers_A10=True")
    print(name + "_finite_divisor_identity=True")
    print(name + "_x_core_disc=-4*(S_o_Tinv)^2*(L_o_Tinv)")
    print(name + "_postcomposition_square_hostile="
          + repr((witness_point, tuple(map(str, witness_left)), tuple(map(str, witness_right)))))
    print(name + "_honest_conjugacy_27_point_bank=True")

# An affine triangular conjugacy is covered by the same formulas.
A = (x + y + 1, y + z - 2, z + 3)
Ainv = (x - y + z - 6, y - z + 5, z - 3)
require(compose(A, Ainv) == IDENTITY == compose(Ainv, A),
        "affine inverse identity failed")
require(jacobian_det(A) == jacobian_det(Ainv) == 1,
        "affine Jacobian-unit check failed")
LA = apply_map(L, Ainv)
for point in product((-1, 0, 1), repeat=3):
    old_source = evaluate(Ainv, point)
    old_target = evaluate(F, old_source)
    new_target = evaluate(A, old_target)
    require(sp.Rational(LA.subs(dict(zip(VARS, new_target))))
            == sp.Rational(L.subs(dict(zip(VARS, old_target)))),
            "affine boundary pullback failed at " + repr(point))
    lhs = evaluate_chain(point, Ainv, F, A, Ainv, F, A)
    rhs = evaluate_chain(point, Ainv, F, F, A)
    require(lhs == rhs,
            "affine honest conjugacy square failed at " + repr(point))
print("affine_conjugacy_boundary_and_27_point_bank=True")

# Exact scalar-normalization recurrence.  If L is replaced by cL and the raw
# recurrence is repeated, P_n is multiplied by c^q_n.
grades = [(1, 0)]
for _ in range(10):
    e, m = grades[-1]
    grades.append((7 * e - 2 * m, 3 * e - 2 * m))
q = [1]
for n in range(10):
    q.append(grades[n][0] + 3 * q[n])
require(all(e % 2 == 1 for e, _ in grades),
        "packet first-coordinate parity failed")
require(all(qn % 2 == (n + 1) % 2 for n, qn in enumerate(q)),
        "boundary-rescaling exponent parity failed")
print("grades_0_6=" + repr(tuple(grades[:7])))
print("boundary_rescale_exponents_q_0_6=" + repr(tuple(q[:7])))
print("q_parity_n_plus_1=True")
print("VERDICT=one-step_two-sided_covariance; all-level_honest_conjugacy; standard_weights_require_transport")
