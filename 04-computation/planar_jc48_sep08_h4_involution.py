#!/usr/bin/env python3
"""Exact H4 reflection-set closure and the retained-sheet Euler bound.
All arithmetic is in Z[phi], phi^2=phi+1. No prior implementation is imported.
"""
from itertools import permutations, product
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
import json

GATES = 0

def need(value, message):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(message)

ZERO = (0, 0)

def add(x, y):
    return (x[0] + y[0], x[1] + y[1])

def neg(x):
    return (-x[0], -x[1])

def mul(x, y):
    return (x[0] * y[0] + x[1] * y[1],
            x[0] * y[1] + x[1] * y[0] + x[1] * y[1])

def dot(x, y):
    out = ZERO
    for a, b in zip(x, y):
        out = add(out, mul(a, b))
    return out

def negative(v):
    return tuple(neg(x) for x in v)

def canon(v):
    return min(v, negative(v))

def even(p):
    return sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4)) % 2 == 0

# Scaled unit H4 roots, all of norm squared4. The three families are disjoint.
roots = set()
for i in range(4):
    for sign in (-1, 1):
        row = [ZERO] * 4
        row[i] = (2 * sign, 0)
        roots.add(tuple(row))
need(len(roots) == 8, 'coordinate roots')
for signs in product((-1, 1), repeat=4):
    roots.add(tuple((s, 0) for s in signs))
need(len(roots) == 24, 'hypercube roots')
for p in permutations(range(4)):
    if not even(p):
        continue
    for signs in product((-1, 1), repeat=3):
        row = (ZERO, (signs[0], 0), (0, signs[1]), (-signs[2], signs[2]))
        roots.add(tuple(row[i] for i in p))
roots = sorted(roots)
need(len(roots) == 120, 'complete root count')
for r in roots:
    need(dot(r, r) == (4, 0), 'root norm')
    need(negative(r) in roots, 'opposite root')
lines = sorted({canon(r) for r in roots})
need(len(lines) == 60, 'root lines')
root_index = {r: i for i, r in enumerate(roots)}
line_index = {r: i for i, r in enumerate(lines)}

def reflect(r, v):
    q = dot(r, v)
    out = []
    for ri, vi in zip(r, v):
        h = mul(q, ri)
        need(h[0] % 2 == 0 and h[1] % 2 == 0, 'integral reflection coefficients')
        out.append(add(vi, (-h[0] // 2, -h[1] // 2)))
    result = tuple(out)
    need(result in root_index, 'root system reflection closure')
    return result

root_action = [tuple(root_index[reflect(r, v)] for v in roots) for r in lines]
line_action = [tuple(line_index[canon(roots[a[root_index[v]]])] for v in lines) for a in root_action]

def compose(a, b):
    return tuple(a[j] for j in b)

def power(a, n):
    out = tuple(range(len(a)))
    for _ in range(n):
        out = compose(out, a)
    return out

for a, b in zip(root_action, line_action):
    need(set(a) == set(range(120)) and power(a, 2) == tuple(range(120)), 'actual reflection involution')
    need(set(b) == set(range(60)) and power(b, 2) == tuple(range(60)), 'line reflection involution')
# This literal simple tuple has a chain3,3,5 Gram matrix after signs orient
# adjacent dots negatively. Its exact representation has order14400 below.
simple = (0, 1, 45, 55)
expected = {(0, 1): ((2, 0), (-2, 0)), (1, 2): ((2, 0), (-2, 0)),
            (2, 3): ((0, 2), (0, -2))}
for i in range(4):
    for j in range(i + 1, 4):
        need(dot(lines[simple[i]], lines[simple[j]]) in expected.get((i, j), (ZERO,)), 'simple Gram matrix')
        n = 5 if (i, j) == (2, 3) else 3 if j == i + 1 else 2
        ab = compose(root_action[simple[i]], root_action[simple[j]])
        need(power(ab, n) == tuple(range(120)) and ab != tuple(range(120)), 'Coxeter product relation')
generators = [root_action[i] for i in simple]
identity = tuple(range(120))
group = {identity}
todo = deque([identity])
while todo:
    a = todo.popleft()
    for b in generators:
        c = compose(a, b)
        if c not in group:
            group.add(c)
            todo.append(c)
need(len(group) == 14400, 'faithful root image order')
need(all(a in group for a in root_action), 'all sixty reflections in generated image')
# The root-line action is used only for conjugation/closed reflection sets.
# Its kernel is the scalar sign; the120-root action above avoids that loss.
line_group = {tuple(line_index[canon(roots[g[root_index[v]]])] for v in lines) for g in group}
need(len(line_group) == 7200, 'line action loses the scalar sign')
for i in range(60):
    for j in range(60):
        lhs = compose(root_action[i], compose(root_action[j], root_action[i]))
        need(lhs == root_action[line_action[i][j]], 'conjugation equals reflected root line')

FULL = (1 << 60) - 1

def members(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask -= bit

def close(mask):
    vals = list(members(mask))
    todo = deque(vals)
    while todo:
        i = todo.popleft()
        for j in vals[:]:
            for k in (line_action[i][j], line_action[j][i]):
                bit = 1 << k
                if not mask & bit:
                    mask |= bit
                    vals.append(k)
                    todo.append(k)
                    if mask == FULL:
                        return mask
    return mask

# No parabolic, stabilizer, degree, cycle-type or table filter is imposed.
# Enumerate every nonempty conjugation-closed subset containing fixed line0.
closed = {1}
todo = deque([1])
extensions = 0
while todo:
    old = todo.popleft()
    for i in range(60):
        if old & (1 << i):
            continue
        new = close(old | (1 << i))
        extensions += 1
        need(new & old == old and new & (1 << i), 'closure retains its seed')
        if new not in closed:
            closed.add(new)
            todo.append(new)
for mask in closed:
    need(mask & 1, 'distinguished reflection retained')
    values = list(members(mask))
    need(all(mask & (1 << line_action[i][j]) for i in values for j in values), 'final conjugation closure')
histogram = Counter(mask.bit_count() for mask in closed)
need(len(closed) == 221 and extensions == 11727, 'complete closure universe')
need(max(mask.bit_count() for mask in closed if mask != FULL) == 16, 'proper reflection maximum')

# Independent benchmark from Douglass--Pfeiffer--Roehrle, arXiv1101.5893v3,
# Table9, printedp15. Tuples(type,class_size,number_of_reflections).
# A class contributes class_size*reflection_count/60 containing a fixed root.
table = [('A1',60,1), ('A1^2',450,2), ('A2',200,3), ('I2(5)',72,5),
         ('A1A2',600,4), ('I2(5)A1',360,6), ('A3',300,6), ('H3',60,15),
         ('A1^3',300,3), ('H4',1,60), ('H3A1',60,16), ('I2(5)^2',36,10),
         ('A4',60,10), ('A2^2',100,6), ('D4',25,12), ('A1^4',75,4)]
benchmark = Counter()
for name, size, count in table:
    need(size * count % 60 == 0, 'table incidence integer')
    benchmark[count] += size * count // 60
need(histogram == benchmark, 'independent complete table fingerprint')

# Actual transitive conjugation action on60 reflections attains fixed ratio4/15.
orbit = {0}
todo = deque([0])
while todo:
    i = todo.popleft()
    for s in simple:
        j = line_action[s][i]
        if j not in orbit:
            orbit.add(j)
            todo.append(j)
need(len(orbit) == 60, 'reflection class transitive')
fix_counts = [sum(i == j for i, j in enumerate(a)) for a in line_action]
need(set(fix_counts) == {16}, 'sharp actual fixed-point ratio')
need(Fraction(16,60) == Fraction(4,15), 'fixed fraction reduction')
# Appending one fixed letter violates that fraction: transitivity is essential.
need(Fraction(17,61) > Fraction(4,15), 'nontransitive fixed-letter hostile')
# General scalar identity, and an explicit finite sanity universe.
# 3D-8k =13D/15+8(4D/15-k); no bounded scan supplies the all-D proof.
for d in range(2, 241):
    for k in range(0, (4*d)//15 + 1):
        need(3*d - 8*k == Fraction(13*d,15) + 8*(Fraction(4*d,15)-k), 'Euler decomposition')
        need(d - 2*k > 0 and 3*d - 8*k > 1, 'three-node contradiction')
need(3*60 - 8*16 == 52, 'sharp action Euler hostile')
report = {'roots':120,'root_lines':60,'simple':simple,'root_image_order':len(group),
          'line_image_order':len(line_group),'closed_sets_containing_one_reflection':len(closed),
          'extension_operations':extensions,'size_histogram':sorted(histogram.items()),
          'largest_proper_reflection_set':16,'sharp_fixed_ratio':'4/15',
          'closed_set_sha256':sha256(json.dumps(sorted(closed)).encode()).hexdigest(),
          'Euler_lower_bound':'13D/15','scalar_sanity_degree_range':[2,240]}
semantic = json.dumps(report,sort_keys=True,separators=(',',':')).encode()
print(json.dumps(report,sort_keys=True,indent=2))
print('semantic_sha256='+sha256(semantic).hexdigest())
print('always_active_exact_gates='+str(GATES))
