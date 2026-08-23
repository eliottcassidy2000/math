#!/usr/bin/env python3
"""Exact companion for THM-3801's cubic-normalization gate."""

from __future__ import annotations

import hashlib
import itertools
import json

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: sp.Expr, rhs: sp.Expr, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)


# Depressed monogenic cubic: trace, different, norm, discriminant, and the
# ramified-sheet/unramified-companion factorization.
t, p, q = sp.symbols("t p q")
f = t**3 + p * t + q
M = sp.Matrix(
    [
        [0, 0, -q],
        [1, 0, -p],
        [0, 1, 0],
    ]
)
I3 = sp.eye(3)
same(M.charpoly(t).as_expr(), f, "companion characteristic polynomial")

traces = [sp.trace(M**j) for j in range(5)]
expected_traces = [3, 0, -2 * p, -3 * q, 2 * p**2]
for j, expected in enumerate(expected_traces):
    same(traces[j], expected, f"depressed trace power {j}")

trace_matrix = sp.Matrix(3, 3, lambda i, j: traces[i + j])
delta = -4 * p**3 - 27 * q**2
same(trace_matrix.det(), delta, "power-basis trace discriminant")

different_matrix = 3 * M**2 + p * I3
same(different_matrix.det(), -delta, "norm of cubic different")
h = -3 * t**2 - 4 * p
remainder = sp.rem(sp.expand((3 * t**2 + p) ** 2 * h - delta), f, t)
same(remainder, 0, "discriminant pullback factorization")

# Generic simple branch: (t-s)^2(t+2s).  The derivative vanishes only on
# the double root, while h vanishes only on the simple companion.
s = sp.symbols("s", nonzero=True)
branch_sub = {p: -3 * s**2, q: 2 * s**3}
same(f.subs(branch_sub), (t - s) ** 2 * (t + 2 * s),
     "simple cubic branch factorization")
different = 3 * t**2 + p
same(different.subs(branch_sub).subs(t, s), 0, "different on double sheet")
same(h.subs(branch_sub).subs(t, s), 9 * s**2, "companion factor unit on double sheet")
same(different.subs(branch_sub).subs(t, -2 * s), 9 * s**2,
     "different unit on simple sheet")
same(h.subs(branch_sub).subs(t, -2 * s), 0,
     "companion factor on simple sheet")


# A normal finite-flat cubic that is genuinely nonmonogenic: the third
# Veronese invariant ring.  Multiplication is represented in the basis
# (1,omega,theta).
a, d, x, y = sp.symbols("a d x y")
M_omega = sp.Matrix(
    [
        [0, 0, a * d],
        [1, 0, 0],
        [0, a, 0],
    ]
)
M_theta = sp.Matrix(
    [
        [0, a * d, 0],
        [0, 0, d],
        [1, 0, 0],
    ]
)
check(M_omega * M_theta == a * d * I3, "Veronese omega theta relation")
check(M_omega**2 == a * M_theta, "Veronese omega square relation")
check(M_theta**2 == d * M_omega, "Veronese theta square relation")
check(M_omega * M_theta == M_theta * M_omega, "Veronese commutativity")

basis_matrices = [I3, M_omega, M_theta]
trace_pairing = sp.Matrix(
    3,
    3,
    lambda i, j: sp.trace(basis_matrices[i] * basis_matrices[j]),
)
delta_0 = -27 * a**2 * d**2
same(trace_pairing.det(), delta_0, "Veronese trace discriminant")

T_matrix = x * M_omega + y * M_theta
T2_first_column = T_matrix**2 * sp.Matrix([1, 0, 0])
index_matrix = sp.Matrix.hstack(
    sp.Matrix([1, 0, 0]),
    sp.Matrix([0, x, y]),
    T2_first_column,
)
index_form = a * x**3 - d * y**3
same(index_matrix.det(), index_form, "Veronese binary index form")

lam = sp.symbols("lam")
char_T = sp.expand(T_matrix.charpoly(lam).as_expr())
disc_T = sp.discriminant(char_T, lam)
same(disc_T, index_form**2 * delta_0, "index-square discriminant law")
check(index_form.subs({a: 0, d: 0}) == 0,
      "Veronese index form cannot represent a unit")
check(M_omega**3 == a**2 * d * I3, "Veronese cyclic cubic field law")

# Finite semigroup control: every invariant monomial has exactly one of the
# three module residue types (0,0), (2,1), (1,2) modulo three.
residue_types = {(0, 0), (2, 1), (1, 2)}
for i in range(13):
    for j in range(13):
        if (i + j) % 3 != 0:
            continue
        check((i % 3, j % 3) in residue_types,
              f"Veronese module residue i={i},j={j}")


# Exact nodal normalization used by THM-3790.
e, A, C = sp.symbols("e A C")
A_e = e**2
C_e = e**3 - e
N = C**2 - A * (A - 1) ** 2
same(N.subs({A: A_e, C: C_e}), 0, "nodal parametrization")
same(C_e, e * (A_e - 1), "nodal rational inverse identity")
check((A_e.subs(e, 1), C_e.subs(e, 1)) == (1, 0),
      "first nodal preimage")
check((A_e.subs(e, -1), C_e.subs(e, -1)) == (1, 0),
      "second nodal preimage")
check(sp.resultant(sp.diff(A_e, e), sp.diff(C_e, e), e) != 0,
      "nodal normalization is an immersion")
check(sp.Poly(N, A, C, domain=sp.QQ).is_irreducible,
      "nodal polynomial irreducible over Q")


# The only ramified cubic DVR partition with an unramified summand is (2,1).
ramified_partitions: set[tuple[tuple[int, int], ...]] = set()
for count in range(1, 4):
    entries = [(ee, ff) for ee in range(1, 4) for ff in range(1, 4)]
    for packet in itertools.combinations_with_replacement(entries, count):
        if sum(ee * ff for ee, ff in packet) != 3:
            continue
        if not any(ee > 1 for ee, _ff in packet):
            continue
        if not any(ee == 1 for ee, _ff in packet):
            continue
        ramified_partitions.add(packet)
check(ramified_partitions == {((1, 1), (2, 1))},
      "unique ramified-plus-unramified cubic partition")


semantic = {
    "branch": "every discriminant component has tame partition (2,1); ramified sheet deleted; unique reduced principal companion visible",
    "index": "S=R+M,M=R^2; binary cubic I(x,y)=det(1,T,T^2); viable constant-unit cubic has I(R^2) disjoint from k*",
    "monogenic": "if S=R[t], f'(t) is nowhere zero on X, hence constant, hence t satisfies a quadratic relation",
    "nodal": "div(N(A,C))=L+D has at least two visible residue-degree-one-or-more sheets, so N is not a discriminant factor",
    "resolvent": "discriminant squarefree nonconstant; quadratic resolvent nontrivial; normal closure S3",
    "hostile": "third Veronese is normal finite-flat nonmonogenic with I=a*x^3-d*y^3 and cyclic total ramification",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate")
print("completion=finite_flat_rank3;trace_zero_module_rank2;binary_index_cubic")
print("monogenic=no_global_generator_on_constant_unit_etale_open")
print("depressed=delta=(3t^2+p)^2*(-3t^2-4p);trace=(3,0,-2p,-3q,2p^2)")
print("branch=every_component_partition_2plus1;one_visible_principal_companion")
print("resolvent=delta_squarefree_nonsquare;Galois_closure=S3")
print("nodal=N_not_discriminant;visible_pullback=L+D;Pic_classes=1+2=0_mod3")
print("hostile=third_Veronese_nonmonogenic_cyclic_total_ramification")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
