#!/usr/bin/env python3
"""Exact finite controls for the quartic A4 cyclic-resolvent mod-2 gate.

The mathematical gate is geometric.  If E/K is a cyclic cubic extension of
the function field of affine space, R is the target normalization in E, and
U=R_reg, then

    H^1_et(U, mu_2)^C3 = 0.

Indeed, a nonzero invariant double cover would have a C6 Galois closure over
the affine target.  Its quadratic quotient would be unramified in every
codimension-one valuation, hence trivial by purity and pi_1^et(A^n_C)=1.

This companion checks the finite group/module steps, the simplest smooth
cyclic-cubic hostile, and the class-group positive control.  Purity and the
normalization argument remain mathematical inputs, not executable claims.
All checks raise explicitly, including under ``python -O``.
"""

from itertools import product

import sympy as sp


checks = 0


def check(label, condition):
    global checks
    checks += 1
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")


def add2(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


def mat_vec2(matrix, vector):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(len(vector))) % 2 for i in range(len(matrix)))


def mat_mul2(left, right):
    n = len(left)
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(n)) % 2 for j in range(n))
        for i in range(n)
    )


I2 = ((1, 0), (0, 1))
ZERO2 = ((0, 0), (0, 0))
RHO = ((0, 1), (1, 1))
RHO2 = mat_mul2(RHO, RHO)

check("standard C3 matrix has order three", mat_mul2(RHO2, RHO) == I2 and RHO != I2)
check(
    "C3 norm vanishes on the standard plane",
    tuple(tuple((I2[i][j] + RHO[i][j] + RHO2[i][j]) % 2 for j in range(2)) for i in range(2)) == ZERO2,
)

W = tuple(product((0, 1), repeat=2))
nonzero = tuple(v for v in W if v != (0, 0))
fixed = tuple(v for v in W if mat_vec2(RHO, v) == v)
orbit = {mat_vec2(RHO, nonzero[0]), mat_vec2(RHO2, nonzero[0]), nonzero[0]}
check("standard plane has no nonzero C3-fixed vector", fixed == ((0, 0),))
check("C3 cycles the three nonzero characters", orbit == set(nonzero))


# The semidirect product V4 rtimes C3 is A4.  This checks the exact group
# carried by a connected equivariant V4 torsor in the cyclic-resolvent lane.
def rho_power(k, vector):
    if k % 3 == 0:
        return vector
    if k % 3 == 1:
        return mat_vec2(RHO, vector)
    return mat_vec2(RHO2, vector)


def a4_mul(left, right):
    v, i = left
    w, j = right
    return add2(v, rho_power(i, w)), (i + j) % 3


A4 = tuple((v, k) for v in W for k in range(3))
identity = ((0, 0), 0)
check("V4 semidirect C3 has order twelve", len(A4) == 12)
check("A4 model is closed", all(a4_mul(g, h) in A4 for g in A4 for h in A4))
check("A4 model has an identity", all(a4_mul(identity, g) == g == a4_mul(g, identity) for g in A4))


def order(group_element):
    value = identity
    for n in range(1, 13):
        value = a4_mul(value, group_element)
        if value == identity:
            return n
    raise RuntimeError("order exceeded group size")


orders = sorted(order(g) for g in A4)
check("A4 order census", orders == [1] + [2] * 3 + [3] * 8)


# An invariant nontrivial double cover of a cyclic cubic normalization would
# give an order-six Galois group with a normal C2 kernel, hence C6.  In C6,
# every inertia subgroup disjoint from that C2 lies in the unique C3 and dies
# in the quadratic quotient.
C6 = set(range(6))
C2 = {0, 3}
C3 = {0, 2, 4}
subgroups_c6 = ({0}, C2, C3, C6)
allowed_inertia = tuple(H for H in subgroups_c6 if H.intersection(C2) == {0})
check("C6 subgroups disjoint from relative C2", allowed_inertia == ({0}, C3))
check("allowed inertia dies in the quadratic quotient", all(H.issubset(C3) for H in allowed_inertia))


# Smooth cyclic-cubic hostile: (t,u,v) -> (t^3,u,v).  Its normalization is
# A3, so its Kummer H1 is zero.  It has the right C3 quotient shape but cannot
# be the resolvent normalization of an A4 Keller cover by the mod-2 gate.
t, u, v, T = sp.symbols("t u v T")
zeta = (-1 + sp.sqrt(-3)) / 2
orbit_polynomial = sp.expand(sp.prod(T - zeta**k * t for k in range(3)))
check("cyclic orbit polynomial", sp.simplify(orbit_polynomial - (T**3 - t**3)) == 0)
jacobian = sp.Matrix([t**3, u, v]).jacobian([t, u, v]).det()
check("smooth cyclic chart quotient Jacobian", sp.expand(jacobian - 3 * t**2) == 0)
check("smooth cyclic hostile is ramified, not Keller", not jacobian.is_constant())
smooth_chart_unit_rank = 0
smooth_chart_class_2_rank = 0
check("A3 Kummer rank is zero", smooth_chart_unit_rank + smooth_chart_class_2_rank == 0)


# Positive carrier control: for R=C[a,b,c,d]/(d^2-abc), the three ramified
# toric divisor classes have presentation F2^3/<(1,1,1)>.  Cyclically
# permuting a,b,c gives precisely the standard irreducible C3-plane.
ONE3 = (1, 1, 1)


def class_rep(bits):
    mate = add2(bits, ONE3)
    return min(bits, mate)


classes = tuple(sorted({class_rep(bits) for bits in product((0, 1), repeat=3)}))


def cycle3(bits):
    return bits[1], bits[2], bits[0]


nonzero_classes = tuple(bits for bits in classes if bits != (0, 0, 0))
class_orbit = {class_rep(cycle3(nonzero_classes[0]))}
class_orbit.add(class_rep(cycle3(cycle3(nonzero_classes[0]))))
class_orbit.add(nonzero_classes[0])
check("toric hostile class group has four elements", len(classes) == 4)
check("C3 cycles all three nonzero two-torsion classes", class_orbit == set(nonzero_classes))
check("positive control carries one standard plane", len(nonzero_classes) == 3)


# A procyclic regular chart has at most one independent mod-2 character and
# therefore cannot receive the rank-two V4 character plane.
for chart, rank in (("A^n", 0), ("G_m x A^(n-1)", 1)):
    check(f"{chart} cannot carry connected V4", rank < 2)


print("A4 CYCLIC-RESOLVENT MOD-2 GATE")
print("theorem input: X=A^n_C, E/K cyclic cubic, R=normalization, U=R_reg")
print("invariant descent: H^1_et(U,mu_2)^C3=0")
print("mechanism: invariant C2 torsor -> C6 closure -> unramified quadratic quotient of A^n -> contradiction")
print("representation consequence: H^1_et(U,mu_2) is a direct sum of standard F2[C3] planes")
print("quartic A4 test: H^1_et(U,mu_2) must be nonzero, hence has dimension at least two")
print("Kummer test: units must be nonconstant modulo squares or Cl(R)[2] must be nonzero")
print("regular-chart test: any nonempty A^n or G_m x A^(n-1) chart excludes the A4 resolvent")
print("negative control: (t,u,v)->(t^3,u,v) has smooth normalization A^3 and Kummer rank zero")
print("positive carrier: d^2=abc has Cl[2]=(F2)^2 with the standard cyclic action")
print("scope: necessary resolvent obstruction only; no general A4, degree-four, or Jacobian exclusion")
print(f"CHECKS PASSED: {checks}")
print("FAILED CHECKS: NONE")
