#!/usr/bin/env python3
"""Exact checks for AP square-sum triangles and the Berggren bridge.

No external packages; explicit exceptions remain active under python -O.
Theorems are proved in the companion note. Counts are finite checks only.
"""
from collections import defaultdict
from itertools import combinations
from math import gcd, isqrt, lcm, pi, sin, sqrt

CHECKS = 0

def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)

def ppt(m, n):
    return n*n-m*m, 2*m*n, n*n+m*m

def children(m, n):
    return (('L', n, 2*n-m), ('M', n, 2*n+m), ('R', m, n+2*m))

def admissible(a, b, c):
    return c*c > 4*a*b

def image(a, b, c, g=1):
    return tuple(2*g*g*x for x in (c*c-4*a*b, c*c, c*c+4*a*b))

def invert(vertices):
    x, y, z = vertices
    check(0 < x < y < z and x+z == 2*y, 'domain')
    sums = (x+y, x+z, y+z)
    roots = tuple(isqrt(s) for s in sums)
    check(tuple(r*r for r in roots) == sums, 'square edges')
    check(all(r % 2 == 0 for r in roots), 'even roots')
    u, v, w = (r//2 for r in roots)
    check((u+w) % 2 == 0, 'half-root parity')
    p, q = (w+u)//2, (w-u)//2
    g = gcd(p, q)
    check(v % g == 0, 'hypotenuse content')
    return p//g, q//g, v//g, g

def validate(a, b, c, g):
    vertices = image(a, b, c, g)
    A, B = max(a, b), min(a, b)
    check(invert(vertices) == (A, B, c, g), 'bijection inverse')
    check(all(gcd(x, y) == 2*g*g for x,y in combinations(vertices, 2)), 'vertex pairwise gcd')
    roots = (2*g*(A-B), 2*g*c, 2*g*(A+B))
    check(all(gcd(x, y) == 2*g for x,y in combinations(roots, 2)), 'root pairwise gcd')
    check((vertices[1]-vertices[0]) % (96*g*g) == 0, 'gap divisibility')
    check(isqrt(vertices[0]//2)**2 != vertices[0]//2, 'left endpoint type')
    check(isqrt(vertices[2]//2)**2 != vertices[2]//2, 'right endpoint type')
    x,y,z=vertices
    check([abs(y-z), abs(z-x)] == [abs(y-x), abs(x-z)], 'midpoint rotation edge differences')

# Complete positive primitive Euclid census by hypotenuse.
LIMIT_C = 1_000_000
primitive = []
for n in range(2, isqrt(LIMIT_C)+1):
    for m in range(1, n):
        if n*n+m*m > LIMIT_C:
            break
        if (m+n) % 2 and gcd(m,n) == 1:
            a,b,c=ppt(m,n)
            primitive.append((a,b,c,m,n))
admitted = [p for p in primitive if admissible(*p[:3])]
for a,b,c,m,n in admitted:
    validate(a,b,c,1)
    choices = [name for name,m1,n1 in children(m,n) if admissible(*ppt(m1,n1))]
    check(choices == (['L'] if b>a else ['R']), 'unique admitted child')
    check(not admissible(*ppt(n,2*n+m)), 'middle-child obstruction')
print(f'Primitive PPT census C<={LIMIT_C}: {len(primitive)}; admissible {len(admitted)}; ratio {len(admitted)/len(primitive):.9f}')

# Complete rooted ternary tree; exact predicates, no floating thresholds.
level = [('',1,2)]
seen = set()
tree_counts = []
for depth in range(11):
    count = 0
    nxt = []
    for word,m,n in level:
        check((m,n) not in seen, 'tree duplicate')
        seen.add((m,n))
        a,b,c=ppt(m,n)
        good=admissible(a,b,c)
        count += good
        if good:
            check(word.endswith('LLL') or word.endswith('RRR'), 'suffix necessity')
        if word.endswith('LLLL') or word.endswith('RRRR'):
            check(good, 'suffix sufficiency')
        # Four repeated specified letters work from every node, including root.
        ml,nl=m,n
        for _ in range(4):
            ml,nl=nl,2*nl-ml
        check(admissible(*ppt(ml,nl)), 'uniform L entry')
        check(admissible(*ppt(m,n+8*m)), 'uniform R entry')
        if word and word[-1] in 'LR':
            letter=word[-1]
            run=len(word)-len(word.rstrip(letter))
            if run <= 2:
                check(not good, 'entry cannot occur before third letter')
            if run >= 4:
                check(good, 'entry occurs by fourth letter')
        for name,m1,n1 in children(m,n):
            nxt.append((word+name,m1,n1))
    if depth >= 4:
        check(2*3**(depth-4) <= count <= 2*3**(depth-3), 'depth-density bounds')
    tree_counts.append((depth,len(level),count))
    level=nxt
print('Tree depth,total,admissible:', tree_counts)
print(f'All {len(seen)} nodes through depth 10 pass exact branch and suffix controls.')

# Complete direct enumeration of graph triangles, independent of Euclid parameters.
N = 1_000_000
brute = set()
for q in range(2, isqrt(2*N)+1, 2):
    y=q*q//2
    if y>=N:
        break
    for p in range(isqrt(y)+1,q):
        x=p*p-y
        z=2*y-x
        if z>N:
            continue
        r=isqrt(y+z)
        if r*r==y+z:
            brute.add((x,y,z))
generated=set()
for a,b,c,m,n in admitted:
    base=image(a,b,c)
    if base[2]>N:
        continue
    for g in range(1,isqrt(N//base[2])+1):
        v=tuple(g*g*t for t in base)
        generated.add(v)
        validate(a,b,c,g)
check(brute==generated, 'independent graph census mismatch')
ordered=sorted(brute,key=lambda v:(v[2],v))
check(ordered[0]==(482,3362,6242), 'minimum')
print(f'Graph-triangle census max<={N}: {len(brute)}; Euclid-square-scale set agrees exactly.')
print('First five by largest vertex:',ordered[:5])

# Primitive and all-content counting tables; finite evidence, not proof of asymptotics.
for bound in (10**6,10**8,10**10,10**12):
    bases=[image(a,b,c)[2] for a,b,c,m,n in admitted if image(a,b,c)[2]<=bound]
    all_count=sum(isqrt(bound//z) for z in bases)
    print(f'max<={bound}: base count {len(bases)}, all-content count {all_count}, base/sqrt(N)={len(bases)/sqrt(bound):.9f}')

# Earliest two primitive representations sharing a hypotenuse and center.
fibres=defaultdict(list)
for a,b,c,m,n in admitted:
    fibres[c].append((a,b,c))
first_shared=min(c for c,ps in fibres.items() if len(ps)>1)
check(first_shared==145, 'first shared primitive center')
shared=sorted(image(*p) for p in fibres[first_shared])
check(shared==[(14594,42050,69506),(22466,42050,61634)], 'shared-center example')
print('First shared primitive center:',first_shared,shared)

# Arbitrarily large windmills are proved in the note; test k=1,...,8 exactly.
for k in range(1,9):
    ps=[ppt(j,j+1) for j in range(4,4+k)]
    H=lcm(*(p[2] for p in ps))
    triples=[image(*p,H//p[2]) for p in ps]
    centers={t[1] for t in triples}
    vertices={x for t in triples for x in t}
    edges={tuple(sorted(e)) for t in triples for e in combinations(t,2)}
    check(centers=={2*H*H}, 'windmill center')
    check(len(vertices)==2*k+1 and len(edges)==3*k, 'windmill distinct triangles')
    check(all(isqrt(a+b)**2==a+b for a,b in edges), 'windmill square edges')
print('Common-center windmill constructions: k=1..8 verified (subgraphs, not claimed induced).')

# Explicit hostile controls.
check(not admissible(3,4,5), 'root must fail positivity')
check(image(16,63,65)==(386,8450,16514), 'not first-vertex minimality')
check(tuple(2*x for x in (482,3362,6242)) not in brute, 'nonsquare scaling must fail')
print('Hostiles: inadmissible root, first-vertex counterexample, nonsquare scaling, gcd-1 exclusion retained.')
print(f'PASS: {CHECKS} explicit checks.')
