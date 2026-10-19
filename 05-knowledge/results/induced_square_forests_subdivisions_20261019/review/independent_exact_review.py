#!/usr/bin/env python3
"""Independent finite controls; explicit checks survive Python -O."""
from fractions import Fraction
from itertools import product
from math import isqrt
import json

CHECKS = 0

def require(test, msg):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(msg)

def square(t):
    return t >= 0 and isqrt(t) ** 2 == t

def tau_sieve(n):
    a = [0] * (n + 1)
    for d in range(1, n + 1):
        for k in range(d, n + 1, d):
            a[k] += 1
    return a

def interval_graph(M):
    return {v: {r*r-v for r in range(isqrt(M+v), isqrt(2*M+v)+1)
                if M <= r*r-v <= 2*M and r*r-v != v}
            for v in range(M, 2*M+1)}

def factor_neighbors(a, b, M):
    if a > b: a, b = b, a
    d = b-a
    ans = set()
    for u in range(1, isqrt(d)+1):
        if d % u: continue
        v = d//u
        if v <= u or (v-u) % 2: continue
        r = (v-u)//2
        x = r*r-a
        if M <= x <= 2*M and x not in (a,b): ans.add(x)
    return ans

def append_child(labels, parent, M, H):
    b = labels[parent]
    r = isqrt(M+b)
    while r*r < M+b: r += 1
    while r*r-b <= H:
        x = r*r-b
        if x not in labels and all(not square(x+v) for i,v in enumerate(labels) if i != parent):
            return x
        r += 1
    raise RuntimeError(('greedy exhausted', labels, parent, M, H))

def verify_embedding(parents, labels):
    require(len(labels) == len(set(labels)), 'collision')
    for i in range(len(labels)):
        for j in range(i):
            require(square(labels[i]+labels[j]) == (parents[i] == j), ('induced mismatch',i,j))

def forest_embedding(parents, M):
    # Independent auxiliary-root implementation; 0 is an artificial vertex.
    aug = [-1] + [p+1 if p >= 0 else 0 for p in parents]
    labels = [M]
    for p in aug[1:]: labels.append(append_child(labels, p, M, 2*M))
    verify_embedding(aug, labels)
    out = labels[1:]
    verify_embedding(parents, out)
    return out

report = {'status': 'PASS', 'universe_notes': [
    'Finite arithmetic controls are not an all-n proof.',
    'All ordered parent arrays through n=8 are exhausted; these include every unlabeled forest shape through n=8.',
    'Prime-power maxima are checked through exponent 100; the ratio monotonicity proof covers the remaining tail.'
]}

# Exact fourth powers avoid numerical fourth-root comparisons.
tau = tau_sieve(1000000)
for d in range(1, len(tau)):
    require(tau[d]**4 <= 9**4*d, ('tau quartic failure', d))
report['tau_bound_checked_through'] = len(tau)-1
expected = {2:5, 3:3, 5:2, 7:1, 11:1, 13:1}
prod4 = Fraction(1)
for p,a in expected.items():
    vals = [Fraction((k+1)**4, p**k) for k in range(101)]
    require(vals[a] == max(vals), ('prime maximum', p))
    require(Fraction((a+2)**4,(a+1)**4) < p, ('decreasing tail', p))
    prod4 *= vals[a]
require(prod4 == Fraction(127401984,25025), 'product mismatch')
require(prod4 < 9**4, 'product exceeds constant')
report['prime_power_product_fourth_power'] = str(prod4)

# Independent direct adjacency versus divisor-factor reconstruction.
intervals = list(range(1,81)) + [100,127,128,255,256,511,512,1000]
rows = []
for M in intervals:
    g = interval_graph(M)
    delta = min(map(len,g.values()))
    # Exact sufficient weaker degree lower bound: delta+2 >= sqrt(M)/4.
    require(16*(delta+2)**2 >= M, ('degree lower bound',M))
    lam = 0
    for a in g:
        require(a not in g[a], ('loop',M,a))
        for b in range(a+1,2*M+1):
            direct = g[a] & g[b]
            require(direct == factor_neighbors(a,b,M), ('factor mismatch', M,a,b))
            require(len(direct) <= tau[b-a], ('codegree exceeds tau', M,a,b))
            lam = max(lam,len(direct))
    rows.append({'M':M,'minimum_degree':delta,'maximum_codegree':lam})
report['interval_rows'] = rows

# Entire ordered-forest universe through order eight, at certified quartic height.
forest_counts = {}
for n in range(1,9):
    M = (64*(n+1))**4
    require(isqrt(M)//4-2 > (n-1)*9*(64*(n+1))+n, ('quartic criterion',n))
    count = 0
    for parents in product(*(range(-1,i) for i in range(n))):
        labels = forest_embedding(parents,M)
        require(M <= min(labels) <= max(labels) <= 2*M, 'range')
        count += 1
    forest_counts[n] = count
report['forest_parent_arrays'] = forest_counts

# Longer exact controls, neither statistical nor a purported universal test.
long_rows = []
for n in [16,32,64,128,256]:
    families = {
        'star':[-1]+[0]*(n-1),
        'path':[-1]+list(range(n-1)),
        'binary':[-1]+[(i-1)//2 for i in range(1,n)],
        'isolated':[-1]*n,
        'comb':[(-1 if i==0 else (i-1 if i%2 else max(0,i-2))) for i in range(n)],
    }
    for name,parents in families.items():
        M=(64*(n+1))**4
        labels=forest_embedding(parents,M)
        long_rows.append({'n':n,'family':name,'minimum':min(labels),'maximum':max(labels)})
report['long_controls'] = long_rows

# Exact seed extension: arbitrary fixed seed can be retained, including a triangle.
seed=[482,3362,6242]
labels=seed[:]
parents={}
for step in range(60):
    p=(step*7+1)%len(labels)
    parents[len(labels)] = p
    labels.append(append_child(labels,p,1,10**12))
require(labels[:3] == seed, 'seed moved')
for i in range(3,len(labels)):
    for j in range(i):
        require(square(labels[i]+labels[j]) == (parents[i]==j), 'seed extension mismatch')
report['preserved_seed'] = seed
report['seed_extension_vertices'] = 60
# Difference 1 cannot be difference of two distinct positive squares.
require(factor_neighbors(1,2,1) == set(), 'two-parent hostile')
# Simpler non-induced greedy really can create extra edges.
labels=[1]
for i in range(1,30):
    b=labels[-1]; r=isqrt(max(labels)+b+1)
    if r*r < max(labels)+b+1: r+=1
    labels.append(r*r-b)
extra=[(i,j,labels[i],labels[j]) for i in range(len(labels)) for j in range(i-1) if square(labels[i]+labels[j])]
require(bool(extra), 'ordinary greedy hostile absent')
report['ordinary_greedy_extra_edge'] = extra[0]
report['always_active_checks'] = CHECKS
print(json.dumps(report,indent=2,sort_keys=True))
