#!/usr/bin/env python3
"""Independent exact audit; no assertions, floats, or production imports.

Finite universes are printed. General statements rest on REPORT.md proofs.
"""
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt
from pathlib import Path
from hashlib import sha256
import json
import sympy as s

checks = 0
def require(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)

ROOT = Path(__file__).resolve().parent
for pin in json.loads((ROOT/'SOURCE_PINS.json').read_text()):
    if 'snapshot' in pin:
        require(sha256((ROOT/pin['snapshot']).read_bytes()).hexdigest()==pin['sha256'],
                'source hash '+pin['snapshot'])

def rank(rows):
    if not rows:
        return 0
    a = [[F(x) for x in row] for row in rows]
    r = 0
    for j in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][j]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        d = a[r][j]
        a[r] = [x / d for x in a[r]]
        for i in range(len(a)):
            if i != r:
                d = a[i][j]
                a[i] = [x-d*y for x,y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r

def matmul(a, b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]

def join(a,b):
    return [x+y for x,y in zip(a,b)]

def gaussian(a,b,c,d):
    return a*c-b*d, a*d+b*c

def norm_test(n):
    return all(p%4 != 3 or e%2 == 0 for p,e in s.factorint(n).items())

rep = {}
for n in range(1,1001):
    pairs = [(a,isqrt(n-a*a)) for a in range(isqrt(n)+1)
             if isqrt(n-a*a)**2 == n-a*a]
    require(bool(pairs) == norm_test(n), f'norm criterion {n}')
    if n in (2,5,44,65):
        rep[n] = pairs
require(rep[44] == [], '44 not a norm')
require((1,2) in rep[5], '5 is a norm')
for a,b,c,d in product(range(-5,6), repeat=4):
    x,y=gaussian(a,b,c,d)
    require(x*x+y*y == (a*a+b*b)*(c*c+d*d), 'Brahmagupta norm')

# Exact inherited raw Phi identities, with primitive/opposite-parity filters
# applied only to the primitive-triple consequence.
phi=lambda a,b:(a*a-b*b,2*a*b,a*a+b*b)
for a,b,c,d in product(range(-5,6), repeat=4):
    A,B,C=phi(a,b); D,E,G=phi(c,d)
    require(C*G-A*D-B*E == 2*(a*d-b*c)**2, 'Lorentz determinant')
for a,b in product(range(-20,21), repeat=2):
    if gcd(a,b)==1:
        A,B,C=phi(a,b)
        require(gcd(gcd(A,B),C)==(2 if a%2 and b%2 else 1), 'parity content')

# Rank gate: A is 2x1, B is 2xm, entries in {-1,0,1}.
# L=C beta uses either beta=(1,3^(1/4)) or beta=(1,phi).
# Cross-check every case by direct 2-equation left-kernel compatibility.
candidates = [([[2,0],[0,1]],'dudeney'),
              ([[1,0],[0,1],[1,1]],'golden')]
rank_counts={}
for C,label in candidates:
    m=len(C); passed=failed=0
    for av in product(range(-1,2),repeat=2):
        A=[[av[0]],[av[1]]]
        for bs in product(range(-1,2),repeat=2*m):
            B=[list(bs[:m]),list(bs[m:])]
            D=matmul(B,C)
            gate=rank(A)==rank(join(A,D))
            if av==(0,0):
                direct=not any(x for row in D for x in row)
            else:
                direct=all(av[1]*D[0][j]-av[0]*D[1][j]==0 for j in range(2))
            require(gate==direct, f'{label} labelled rank equivalence')
            passed+=gate; failed+=not gate
            if sum(abs(x) for x in bs)<=1:
                require(rank(join(A,D))==s.Matrix(join(A,D)).rank(), 'independent SymPy rank')
    rank_counts[label]={'pass':passed,'fail':failed,'total':passed+failed}

# Equal relation-space dimension does not preserve the labelled annihilator.
C=[[1,0],[0,1],[1,1]]
good=[[-1,-1,1]]; bad=[[1,-1,0]]
require(rank(good)==rank(bad)==1, 'equal rank hostile')
require(matmul(good,C)==[[0,0]], 'correct labelled relation')
require(matmul(bad,C)==[[1,-1]], 'wrong labelled relation')
# Equality feasibility is not positive-length feasibility.
require(rank([[1]])==rank([[1,-1]]), 'positive-length hostile rank passes')
require(F(-1)<0, 'positive-length hostile unique solution negative')

# Exact target geometry and splitting field checks (not cut coordinates).
x=s.Symbol('x'); alpha=s.real_root(3,4)
require(s.Poly(x**4-3,x).is_irreducible, 'Eisenstein degree four')
G,_=s.polys.numberfields.galois_group(x**4-3,x)
require(G.order()==8 and not G.is_abelian, 'quartic splitting group D4')
require(s.simplify(alpha**2-s.sqrt(3))==0, 'square area')
require(s.simplify((s.I*alpha/2)*s.conjugate(s.I*alpha/2)-s.sqrt(3)/4)==0,
        'complex squared absolute value remains positive')
require(s.simplify((s.I*alpha/2)**2+s.sqrt(3)/4)==0,
        'r sends w c(w) to negative conjugate')
require(s.sqrt(3).is_irrational is True, 'sqrt3 irrational')
require(s.Poly(x**4-3,x).degree()==4, 'alpha rank-two with 1')

print(json.dumps({'status':'PASS','checks':checks,
    'sympy_version':s.__version__,
    'norm_integer_universe':[1,1000],
    'gaussian_identity_box':[-5,5],
    'primitive_parameter_box':[-20,20],
    'representations_nonnegative_ordered':rep,
    'rank_universe':rank_counts,
    'rank_scope':'unconstrained linear equalities only; no embedding or piece count',
    'galois_order':G.order(),
    'field_scope':'target vertices only, not dissection cuts'},sort_keys=True,indent=2))
