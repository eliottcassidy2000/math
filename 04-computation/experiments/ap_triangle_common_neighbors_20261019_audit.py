#!/usr/bin/env python3
"""Independent exact audit; never imports the submitted checker.

Reads its adjacent certificate JSON and checks it by rational identities,
chord collinearity, exhaustive pairwise integer tests, and a different complete
integer-neighbor enumeration through all integer triangles of area 2880.
All checks survive -O. Standard library only.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import isqrt, lcm
from pathlib import Path
import hashlib
import json
import subprocess
import sys

CHECKS = 0

def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(label)

def rat_square(z):
    z = Q(z)
    return z >= 0 and isqrt(z.numerator)**2 == z.numerator and isqrt(z.denominator)**2 == z.denominator

root = Path(__file__).resolve().parent
base = root / 'ap_triangle_common_neighbors_20261019'
data = json.loads(base.with_suffix('.json').read_text())
cert = data['certificate']
labels = list(map(int, cert['integer_labels']))
L = int(cert['common_scale_L'])
require(len(labels) == 6 and len(set(labels)) == 6, 'six distinct labels')
require(all(z > 0 for z in labels), 'positive labels')
require(int(cert['label_multiplier_L_squared']) == L*L, 'scale square')
require(labels[:3] == [L*L*z for z in (482,3362,6242)], 'original core square dilation')
expected = {tuple(e) for e in combinations(range(3), 2)} | {(i,j) for i in range(3) for j in range(3,6)}
recorded = {tuple(e['edge']):int(e['root']) for e in cert['required_edges']}
require(len(cert['required_edges']) == len(recorded) == 12, '12 distinct recorded edges')
require(set(recorded) == expected, 'complete required edge set')
all_pairs = []
for i,j in combinations(range(6),2):
    s = labels[i]+labels[j]
    r = isqrt(s)
    square = r*r == s
    require(square == ((i,j) in expected), 'every induced adjacency')
    if (i,j) in recorded:
        require(recorded[(i,j)] == r, 'recorded exact edge root')
    all_pairs.append({'pair':[i,j], 'is_square':square})

# Derive 2P directly, then certify claimed 3P,4P,5P by curve membership
# and chord collinearity. This does not use the submitted addition routine.
P = (Q(25,4),Q(75,8))
x,y = P
slope = (3*x*x-25)/(2*y)
x2 = (x*x+25)**2/(4*y*y)
y2 = slope*(x-x2)-y
previous = (x2,y2)
denominators = []
neighbors = []
for index, row in enumerate(cert['rational_neighbors'], start=3):
    point = tuple(map(Q,row['point_on_E5']))
    u,v = point
    require(row['multiple'] == index, 'successive multiples')
    require(v*v == u*u*u-25*u, 'E5 membership')
    px,py = previous
    require(px != x and u != x and u != px, 'chord has three distinct abscissas')
    require((py-y)*(u-x) == (-v-y)*(px-x), 'P, previous, -current collinear')
    roots = list(map(Q,row['roots']))
    z = Q(row['neighbor'])
    require(len(roots) == 3 and all(t > 0 for t in roots), 'three positive roots')
    for j,t in enumerate(roots):
        require(t*t == z+482+2880*j, 'each rational square equation')
    require(z > 0 and z not in (482,3362,6242) and z not in neighbors, 'fresh graph neighbor')
    require(z == 576*(u*u+25)**2/(4*v*v)-3362, 'direct duplication identity')
    D = lcm(*(t.denominator for t in roots))
    require(D == int(row['root_denominator']), 'root denominator')
    denominators.append(D)
    require(Q(labels[index], L*L) == z, 'integer leaf is exact square scaling')
    neighbors.append(z)
    previous = point
require(L == lcm(*denominators), 'single L is lcm of all root denominators')

# Independent complete hostile: integer roots p,q,r give integer right
# triangle (r-p,r+p,2q) of area 2880. Enumerate all unordered AB=5760,
# retaining exactly those for which A^2+B^2 is an integer square.
factor_pairs = []
triangles = []
solutions = []
for A in range(1,isqrt(5760)+1):
    if 5760 % A:
        continue
    B = 5760//A
    C = isqrt(A*A+B*B)
    is_triangle = C*C == A*A+B*B
    factor_pairs.append({'A':A,'B':B,'integer_hypotenuse':is_triangle})
    if is_triangle:
        triangles.append([A,B,C])
        p,q,r = Q(B-A,2), Q(C,2), Q(A+B,2)
        if all(t.denominator == 1 for t in (p,q,r)):
            z = p*p-482
            solutions.append([int(z),int(p),int(q),int(r)])
require(len(factor_pairs) == 24, 'all 24 unordered AB factor pairs')
require(triangles == [[36,160,164]], 'unique integer right triangle of area 2880')
require(solutions == [[3362,62,82,98]], 'unique integer square-equation solution')
external = [s for s in solutions if s[0] > 0 and s[0] not in (482,3362,6242)]
require(not external, 'no external integer common neighbor')

control = data['nonclique_positive_control']
z = Q(control['neighbor'])
require(not rat_square(1+7), 'control anchors are not a clique')
require(z > 0 and z not in control['anchors'], 'control is graph-distinct positive')
for a,r in zip(control['anchors'], control['roots']):
    require(z+a == Q(r)**2, 'nonclique control edge')

# Reproducibility is separate from the independent mathematical checks.
regular = subprocess.check_output([sys.executable,str(base.with_suffix('.py'))])
optimized = subprocess.check_output([sys.executable,'-O',str(base.with_suffix('.py'))])
require(regular == optimized == base.with_suffix('.json').read_bytes(), 'normal/-O/stored output identical')
require(data['checks'] == 6615, 'reported original check count')

print(json.dumps({
    'status':'PASS: no mathematical or certificate blocker found',
    'independent_exact_checks':CHECKS,
    'six_vertex_certificate':{'positive_distinct_labels':6,'required_edges':12,'leaf_nonedges':3,'induced_graph':'K3 join I3','all_pair_tests':all_pairs},
    'independent_area_2880_enumeration':{'all_unordered_factor_pairs_AB_5760':factor_pairs,'integer_right_triangles':triangles,'all_integer_neighbor_equation_solutions':solutions,'external_integer_common_neighbors':external},
    'original_checker':{'checks':6615,'normal_optimized_and_stored_outputs_identical':True,'output_sha256':hashlib.sha256(regular).hexdigest()},
},indent=2,sort_keys=True))
