#!/usr/bin/env python3
"""Exact audit for THM-2201's F_13 root-fibre Hasse carrier."""

from itertools import combinations
from math import comb

P = 13
N = 13


def jet(vector):
    return tuple(
        sum(vector[k] * comb(k, j) for k in range(N)) % P
        for j in range(N)
    )


def translate(vector, amount):
    out = [0] * N
    for k, value in enumerate(vector):
        out[(k + amount) % N] = value
    return tuple(out)


def translated_jet(old_jet, amount):
    return tuple(
        sum(comb(amount, j - r) * old_jet[r] for r in range(j + 1)) % P
        for j in range(N)
    )


def determinant_mod(matrix):
    a = [row[:] for row in matrix]
    det = 1
    for col in range(len(a)):
        pivot = next(row for row in range(col, len(a)) if a[row][col] % P)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        pivot_value = a[col][col] % P
        det = det * pivot_value % P
        inv = pow(pivot_value, -1, P)
        a[col] = [(entry * inv) % P for entry in a[col]]
        for row in range(col + 1, len(a)):
            factor = a[row][col] % P
            if factor:
                a[row] = [
                    (a[row][j] - factor * a[col][j]) % P
                    for j in range(len(a))
                ]
    return det % P


pascal = [[comb(k, j) % P for k in range(N)] for j in range(N)]
pascal_det = determinant_mod(pascal)
assert pascal_det == 1

translation_checks = 0
for amount in range(N):
    for k in range(N):
        basis = [0] * N
        basis[k] = 1
        assert jet(translate(basis, amount)) == translated_jet(jet(basis), amount)
        translation_checks += 1

endpoint_checks = 0
zero = (0,) * N
for old in range(N):
    for new in range(N):
        moved = [0] * N
        moved[old] = -1
        moved[new] += 1
        expected = tuple(
            (-comb(old, j) + comb(new, j)) % P for j in range(N)
        )
        assert jet(moved) == expected
        endpoint_checks += 1

occupancy_checks = 0
for size in range(11):
    for subset in combinations(range(N), size):
        indicator = [0] * N
        for k in subset:
            indicator[k] = 1
        assert (jet(indicator)[0] == 0) == (size == 0)
        occupancy_checks += 1

full_fibre_jet = jet((1,) * N)
assert full_fibre_jet == (0,) * 12 + (1,)

x, y = 2, 7
left = ({x}, {x, y})
right = ({x}, {y})
assert set().union(*left) == set().union(*right)
assert set().union(*left[1:]) != set().union(*right[1:])

print("THM-2201 cyclic root-fibre Hasse-jet exact audit")
print("field=F_13; sheets=13; algebra=F_13[epsilon]/(epsilon^13)")
print(f"pascal_det_mod13={pascal_det}")
print(f"basis_translation_checks={translation_checks}")
print(f"endpoint_delete_insert_checks={endpoint_checks}")
print(f"occupancy_subsets_size_at_most_10={occupancy_checks}")
print(f"full_fibre_jet={full_fibre_jet}")
print("aggregate_deletion_hostile_pair=PASS")
print("PASS")
