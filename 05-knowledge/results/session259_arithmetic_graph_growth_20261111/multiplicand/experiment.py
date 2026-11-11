#!/usr/bin/env python3
"""Independent F2 boundary reduction, not the relative-complex proof algorithm."""
import json
from itertools import combinations
from math import prod, isqrt

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def factors(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out

B = 300
primes = [n for n in range(2, B+1) if factors(n) == {n: 1}]
faces = []
def visit(face, first, weight):
    for j in range(first, len(primes)):
        p = primes[j]
        if weight*p > B:
            break
        s = face+(p,)
        faces.append(s)
        visit(s, j+1, weight*p)
visit((), 0, 1)
faces.sort(key=lambda s: (prod(s), len(s), s))
index = {s:i for i,s in enumerate(faces)}
columns = {}
zeros = set()
pairs = {}
for j,s in enumerate(faces):
    col = 0
    if len(s) > 1:
        for k in range(len(s)):
            col ^= 1 << index[s[:k]+s[k+1:]]
    while col:
        pivot = col.bit_length()-1
        if pivot not in columns:
            break
        col ^= columns[pivot]
    if col:
        pivot = col.bit_length()-1
        columns[pivot] = col
        pairs[pivot] = j
    else:
        zeros.add(j)
actual = sorted((len(faces[i])-1,prod(faces[i]),prod(faces[pairs[i]]) if i in pairs else None)
                for i in zeros)
expected = [(0,2,None)]
for m in range(3,B+1,2):
    f = factors(m)
    if all(e == 1 for e in f.values()):
        expected.append((len(f)-1,m,2*m if 2*m <= B else None))
expected.sort()
require(actual == expected, ('barcode', actual, expected))

# Independent reachability using integer bitsets; no factoring in adjacency.
for N in range(2,121):
    weak = {a:[b for b in range(a+1,N+1) if b%a == 0] for a in range(2,N+1)}
    strict = {a:[b for b in weak[a] if b != a*a] for a in weak}
    reach = {}
    for a in range(N,1,-1):
        r = 0
        for b in strict[a]:
            r |= (1<<b) | reach[b]
        reach[a] = r
    for a in range(2,N+1):
        for b in range(a+1,N+1):
            predicted = b%a == 0 and not (b == a*a and factors(a) == {a:1})
            require(bool(reach[a] & (1<<b)) == predicted, ('reach',N,a,b))
    require(sum(map(len,weak.values()))-sum(map(len,strict.values())) == isqrt(N)-1,
            ('square_edge_count',N))
    neighbors = {a:set(weak[a]) for a in weak}
    for a,bs in weak.items():
        for b in bs:
            neighbors[b].add(a)
    unseen = set(weak)
    components = 0
    while unseen:
        components += 1
        stack = [unseen.pop()]
        while stack:
            for b in neighbors[stack.pop()] & unseen:
                unseen.remove(b)
                stack.append(b)
    if N >= 4:
        require(components == 1+sum(N//2 < p <= N for p in primes), ('components',N))
    for p in primes:
        if p*p > N:
            break
        incident = bool(strict[p*p]) or any(p*p in bs for bs in strict.values())
        require(incident == (N >= 2*p*p), ('prime_square_isolation',N,p))

# Face Euler check, independently against divisor-recursive Mobius.
mu = {1:1}
for n in range(2,B+1):
    mu[n] = -sum(mu[d] for d in range(1,n) if n%d == 0)
for N in range(2,B+1):
    chi = -1+sum((-1)**(len(s)-1) for s in faces if prod(s) <= N)
    require(chi == -sum(mu[n] for n in range(1,N+1)), ('Euler',N))
receipt = {
    'status':'FINITE-EXACT PASS',
    'boundary_field':'F2',
    'prime_product_cutoff':B,
    'nonempty_simplices':len(faces),
    'persistent_intervals_including_essential_H0':len(actual),
    'completed_intervals':sum(d is not None for _,_,d in actual),
    'positive_degree_intervals':[x for x in actual if x[0] > 0],
    'weak_graph_and_strict_reachability_cutoff':120,
    'hostiles':{'4':'square diagonal; nonsquarefree topologically silent',
                '12':'nonsquare but topologically silent',
                '15':'H1 birth', '30':'H1 death'},
    'all_order_status':'proved in RESEARCH.md, not inferred from this finite computation'
}
print(json.dumps(receipt,indent=2))
