#!/usr/bin/env python3
"""Finite controls for the classic lattice nonduplication proof.

The universal conclusion is proved in CLASSIC_LATTICE_NONDUPLICATION.md,
not inferred from this bounded census. No arbitrary infinite pieces are
enumerated. Python standard library only; failures remain active under -O.
"""
from itertools import product
from pathlib import Path
import json

checks = 0
def need(condition, name):
    global checks
    checks += 1
    if not condition:
        raise ValueError(name)

def apply(M, p):
    return tuple(sum(M[i][j]*p[j] for j in (0, 1)) for i in (0, 1))

def multiply(A, B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in (0, 1))
                       for j in (0, 1)) for i in (0, 1))

def inverse(M):
    a, b = M[0]
    c, d = M[1]
    determinant = a*d-b*c
    need(abs(determinant) == 1, 'integral inverse exists')
    return ((d//determinant, -b//determinant),
            (-c//determinant, a//determinant))

def rho(p):
    x, y = p
    return max(abs(x), abs(y), abs(x+y))

def hex_size(r):
    return 3*r*r+3*r+1

def translate(points, t):
    return {(p[0]+t[0], p[1]+t[1]) for p in points}

I = ((1, 0), (0, 1))
rotation = ((0, -1), (1, 1))
reflection = ((1, 1), (0, -1))
group = {I}
while True:
    enlarged = group | {multiply(M, g) for M in group for g in (rotation, reflection)}
    if enlarged == group:
        break
    group = enlarged
need(len(group) == 12, 'D6 means order-twelve dihedral group')
inverses = {M: inverse(M) for M in group}
for M in group:
    need(multiply(M, inverses[M]) == I, 'matrix inverse')
    need(all(rho(apply(M, p)) == rho(p) for p in product(range(-4, 5), repeat=2)),
         'hexagonal radius preserved')

square_cases = hex_cases = isometry_cases = partition_cases = 0
for R in range(9):
    box = set(product(range(-R, R+1), repeat=2))
    hexagon = {p for p in box if rho(p) <= R}
    L = 2*R+1
    need(len(box) == L*L, 'square cardinality')
    need(len(hexagon) == hex_size(R), 'hexagonal cardinality')
    for M in group:
        need({apply(M, p) for p in hexagon} == hexagon, 'hexagon invariant under D6')
    for t in product(range(-3, 4), repeat=2):
        minus_t = (-t[0], -t[1])
        shifted = translate(box, minus_t)
        exact_defect = L*L-max(L-abs(t[0]), 0)*max(L-abs(t[1]), 0)
        need(len(shifted-box) == exact_defect, 'exact rectangular overlap')
        need(len(shifted ^ box) == 2*exact_defect, 'square symmetric difference')
        need(exact_defect <= L*(abs(t[0])+abs(t[1])), 'square one-sided boundary bound')
        square_cases += 1
        shifted_hex = translate(hexagon, minus_t)
        d = rho(t)
        shell_bound = 6*R*d+3*d*(d+1)
        need(hex_size(R+d)-hex_size(R) == shell_bound, 'exact outer-shell size')
        need(len(shifted_hex-hexagon) <= shell_bound, 'hexagonal one-sided boundary bound')
        need(len(shifted_hex ^ hexagon) == 2*len(shifted_hex-hexagon), 'hexagonal symmetric difference')
        need(all(rho(p) <= R+d for p in shifted_hex), 'common enlarged hexagon containment')
        hex_cases += 1
        for M in group:
            preimage = {apply(inverses[M], p) for p in shifted_hex}
            expected = translate(hexagon, apply(inverses[M], minus_t))
            need(preimage == expected, 'affine isometry preimage')
            need(len(preimage-hexagon) <= shell_bound, 'isometry boundary bound')
            isometry_cases += 1

    # Infinite quadrant partition, tested only on the finite sets involved.
    # This checks the bookkeeping, not the universal arbitrary-piece theorem.
    label = lambda p: int(p[0] >= 0)+2*int(p[1] >= 0)
    motions = [(I, (3, -2)), (rotation, (-1, 3)),
               (reflection, (0, -3)), (multiply(rotation, reflection), (-2, 1))]
    sets = [{apply(inverses[M], p) for p in translate(hexagon, (-t[0], -t[1]))}
            for M, t in motions]
    left = sum(sum(label(p) == i for p in F) for i, F in enumerate(sets))
    right = len(hexagon)+sum(len(F-hexagon) for F in sets)
    need(left <= right, 'arbitrary-piece finite-intersection inequality control')
    partition_cases += 1

for D in range(21):
    square_gap = 2*(6*D+1)**2-(8*D+1)**2
    hex_gap = 2*hex_size(3*D)-hex_size(4*D)
    need(square_gap == 8*D*D+8*D+1 > 0, 'explicit square contradiction R=3D')
    need(hex_gap == 6*D*D+6*D+1 > 0, 'explicit hexagonal contradiction R=3D')

result = {'status': 'FINITE_EXACT_PASS', 'checks': checks,
          'scope': 'Finite controls for an elementary classic proof; no universal plane or three-dimensional theorem.',
          'square_translation_cases': square_cases,
          'hexagonal_translation_cases': hex_cases,
          'hexagonal_affine_isometry_cases': isometry_cases,
          'quadrant_partition_controls': partition_cases,
          'radii_tested': [0, 8], 'translation_coordinate_range': [-3, 3],
          'D6_order': len(group), 'explicit_cutoff_D_tested': [0, 20],
          'zero_displacement_case_included': True,
          'universal_conclusion_source': 'The finite-boundary proof, not this finite census.'}
text = json.dumps(result, indent=2)+'\n'
Path(__file__).with_name('classic_lattice_nonduplication_checks.json').write_text(text)
print(text, end='')
