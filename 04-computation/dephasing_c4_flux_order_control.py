#!/usr/bin/env python3
"""Exact C4 control for the first phase-sensitive dephasing correction."""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


n = 4
z = sp.symbols("z", nonzero=True)
ii = sp.I

# Unit hoppings around 0-1-2-3-0, with the full Wilson phase on edge 3->0.
h = sp.zeros(n)
for a, b, value in (
    (0, 1, sp.Integer(1)),
    (1, 2, sp.Integer(1)),
    (2, 3, sp.Integer(1)),
    (3, 0, z),
):
    h[a, b] = value
    h[b, a] = 1 / value

basis = [(a, b) for a in range(n) for b in range(n)]
index = {pair: k for k, pair in enumerate(basis)}
k_super = sp.zeros(n * n)

# K=-i[H,-] on the matrix-unit basis E_ab.
for a, b in basis:
    column = index[a, b]
    for m in range(n):
        k_super[index[m, b], column] += -ii * h[m, a]
        k_super[index[a, m], column] += ii * h[b, m]

p_indices = [index[i, i] for i in range(n)]
q_indices = [index[a, b] for a, b in basis if a != b]
b_block = k_super.extract(p_indices, q_indices)
c_block = k_super.extract(q_indices, p_indices)
d_block = k_super.extract(q_indices, q_indices)

a3 = sp.simplify(b_block * d_block * c_block)
a4_raw = sp.simplify(b_block * d_block * d_block * c_block)
delta = sp.simplify(z + 1 / z - 2)
phase_matrix = sp.simplify((a4_raw - a4_raw.subs(z, 1)) / delta)

expected = sp.Matrix(
    [
        [2, -4, 6, -4],
        [-4, 2, -4, 6],
        [6, -4, 2, -4],
        [-4, 6, -4, 2],
    ]
)

require(a3 == sp.zeros(n), "the triangle/third-order term must vanish on C4")
require(phase_matrix == expected, "unexpected fourth-order Wilson matrix")
require(
    sp.simplify(a4_raw - a4_raw.subs(z, 1) - delta * expected) == sp.zeros(n),
    "fourth-order phase factorization failed",
)
row_sums = [sp.simplify(sum(expected.row(i))) for i in range(n)]
require(row_sums == [0, 0, 0, 0], "population conservation failed")

print("C4 STRONG-DEPHASING FLUX CONTROL")
print("convention=K=-i[H,-], unit C4 hoppings, Wilson phase z on edge 3->0")
print("A3_zero=True")
print("delta=z+z^-1-2=2(cos(Phi)-1)")
print(f"phase_matrix={expected.tolist()}")
print(f"row_sums={row_sums}")
print("VERDICT=phase first appears in the fourth Schur-walk coefficient")
