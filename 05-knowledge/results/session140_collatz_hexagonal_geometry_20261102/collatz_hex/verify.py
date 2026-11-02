#!/usr/bin/env python3
"""Exact finite controls for Session140 Collatz geometry; no third-party packages."""
import json


def check(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def C(n):
    return n // 2 if n % 2 == 0 else 3*n+1


def T(n):
    return n // 2 if n % 2 == 0 else (3*n+1)//2


def U(n):
    m = 3*n+1
    while m % 2 == 0:
        m //= 2
    return m


def iterate(f, n, steps):
    for _ in range(steps):
        n = f(n)
    return n


def predC(n):
    return [2*n] + ([(n-1)//3] if n % 6 == 4 else [])


def predT(n):
    return [2*n] + ([(2*n-1)//3] if n % 3 == 2 else [])


def fibre(pred, n, depth):
    layer = {n}
    for _ in range(depth):
        layer = {p for m in layer for p in pred(m)}
    return layer


K = {1:2, 2:1, 4:2, 5:3, 7:4, 8:1}

def selected(n, shortcut=False):
    k = K[n % 9]
    return ((16 if shortcut else 32)*n, (16*n-2**(4-k))//3)


for n in range(1, 10001):
    check(all(C(x) == n for x in predC(n)), ('C inverse', n))
    check(all(T(x) == n for x in predT(n)), ('T inverse', n))
    if n % 3:
        k = K[n % 9]
        check(pow(2, k)*n % 9 in (4,7), ('mod9 port', n))
        for shortcut, f, d in [(False,C,5), (True,T,4)]:
            a,b = selected(n,shortcut)
            check(a != b and a > n and b > n, ('strict selection',n))
            check(a % 3 and b % 3, ('unit selection',n))
            check(iterate(f,a,d) == n == iterate(f,b,d), ('forward check',n))

# Independent route: explicitly enumerate full inverse fibres rather than
# trusting the affine formulas. Both clocks fail uniformly one step sooner.
for n in range(1, 501):
    if n % 3:
        for shortcut, pred, d in [(False,predC,5),(True,predT,4)]:
            check(set(selected(n,shortcut)) <= fibre(pred,n,d), ('fibre',n))
check({n for n in fibre(predC,7,4) if n % 3} == {112}, 'C sharp clock')
check({n for n in fibre(predT,7,3) if n % 3} == {56}, 'T sharp clock')
check(fibre(predC,7,4) == {112,18}, 'C full hostile fibre')
check(fibre(predT,7,3) == {56,18}, 'T full hostile fibre')

rows = []
for root in (1,2,5,7,8,31):
    for shortcut in (False,True):
        layer, seen = {root}, {root}
        for r in range(1,15):
            layer = {p for n in layer for p in selected(n,shortcut)}
            check(len(layer) == 2**r, ('level cardinality',root,r,shortcut))
            check(not (layer & seen), ('cross-level collision',root,r,shortcut))
            seen.update(layer)
        rows.append({'root':root,'map':'T' if shortcut else 'C','depth':14,'leaves':len(layer)})

# Independent odd-map branching, supplied concurrently by the HN lane.
for m in range(1,2000,2):
    if m % 3:
        children = [(2**a*m-1)//3 for a in range(1,7)
                    if (2**a*m) % 9 in (4,7)]
        check(len(children) == 2, ('odd two ports',m))
        check(all(n % 2 and n % 3 and U(n) == m for n in children), ('odd parent',m))
        exponents = (3,5,7) if m % 3 == 2 else (4,6,8)
        increasing = [(2**a*m-1)//3 for a in exponents
                      if (2**a*m) % 9 in (4,7)]
        check(len(increasing)==2 and all(n>m and U(n)==m for n in increasing),
              ('increasing odd subtree',m))

# A concrete directed hexagon routes the ordinary 1,4,2 cycle with length two.
def honeycomb_adj(p,q):
    x,y = p
    xx,yy = q
    return (y == yy and abs(x-xx) == 1) or (x == xx and abs(y-yy) == 1 and (x+min(y,yy)) % 2 == 0)
hex_cycle = [(0,0),(1,0),(2,0),(2,1),(1,1),(0,1)]
check(all(honeycomb_adj(hex_cycle[i],hex_cycle[(i+1)%6]) for i in range(6)), 'hex cycle')
check(honeycomb_adj((3,0),(2,0)), 'incoming port from 8')

# Brick-wall H: horizontal edges, and vertical (x,y)--(x,y+1) for x+y even.
# Contract vertical matching. A matching edge with lower endpoint (x,y)
# corresponds to square-grid point ((x+y)//2,(x-y)//2).
for x in range(-10,11):
    for y in range(-10,11):
        if (x+y) % 2 == 0:
            p = ((x+y)//2,(x-y)//2)
            neigh = {((xx+yy)//2,(xx-yy)//2)
                     for xx,yy in [(x-1,y-1),(x-1,y+1),(x+1,y-1),(x+1,y+1)]}
            expected = {(p[0]-1,p[1]),(p[0]+1,p[1]),(p[0],p[1]-1),(p[0],p[1]+1)}
            check(neigh == expected, ('grid matching contraction',x,y))

out = {'status':'PASS','arithmetic_universe':'1 <= n <= 10000',
       'inverse_fibre_universe':'1 <= n <= 500',
       'sharp_uniform_unit_clocks':{'C':5,'T':4,'hostile_root':7},
       'full_hostile_fibres':{'C_inverse_depth4_at7':sorted(fibre(predC,7,4)),
                             'T_inverse_depth3_at7':sorted(fibre(predT,7,3))},
       'selected_macro_tree_checks':rows,
       'independent_odd_branch_universe':'odd m < 2000, 3 not dividing m',
       'brick_wall_matching_universe':'-10 <= x,y <= 10; x+y even',
       'scope':'finite exact controls, not a convergence test or proof'}
print(json.dumps(out,sort_keys=True,indent=2))
