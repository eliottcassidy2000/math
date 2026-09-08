#!/usr/bin/env python3
"""Exact local/Euler passports; no asserted affine-complement presentation."""
from itertools import permutations, combinations
from collections import deque
import hashlib
import json

gates = 0

def need(value, name):
    global gates
    gates += 1
    if not value:
        raise RuntimeError(name)

def mul(p, q):
    return tuple(p[q[i]] for i in range(len(p)))

def inverse(p):
    out = [0]*len(p)
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)

def cyc(d, *cycles):
    p = list(range(d))
    for c in cycles:
        for i, j in zip(c, c[1:]+c[:1]):
            p[i-1] = j-1
    return tuple(p)

def fixed(p):
    return {i for i, j in enumerate(p) if i == j}

def support(p):
    return set(range(len(p))) - fixed(p)

def lengths(p):
    seen = set()
    out = []
    for i in range(len(p)):
        if i in seen:
            continue
        j, size = i, 0
        while j not in seen:
            seen.add(j)
            size += 1
            j = p[j]
        out.append(size)
    return tuple(sorted(out))

def braid(s, t):
    return mul(mul(s, t), s) == mul(mul(t, s), t)

def local(s, t, A):
    need(braid(s, t), 'actual marked ordinary braid')
    need(A <= fixed(s), 'retained subset is inertia fixed')
    g = mul(s, t)
    B = {g[i] for i in A}
    need(B <= fixed(t), 'simultaneous conjugation and reaccess')
    U, V = A-B, B-A
    O = set(range(len(s))) - (A | B)
    need({t[i] for i in U} <= O, 'ordinary-cusp injection')
    need(2*len(A & B) >= 3*len(A)-len(s), 'uniform actual-cusp count')
    if len(U) == len(O):
        need({t[i] for i in U} == O and {t[i] for i in O} == U,
             'tau exchanges equal blocks')
        need({s[i] for i in O} == V and {s[i] for i in V} == O,
             'sigma exchanges equal blocks')
        need(fixed(s) == A, 'no deleted fixed sheets at equality')
        need(all(k == 1 or k % 2 == 0 for k in lengths(s)),
             'equality gives even nontrivial cycles')
    return len(A & B)

def group(gens):
    identity = tuple(range(len(gens[0])))
    seen, todo = {identity}, deque([identity])
    while todo:
        p = todo.popleft()
        for g in gens:
            q = mul(g, p)
            if q not in seen:
                seen.add(q)
                todo.append(q)
    return seen

def node(s, t, a):
    need(mul(s, t) == mul(t, s), 'node meridians commute')
    need(len(fixed(s)) == len(fixed(t)) == a, 'actual node subsets are full fixed sets')
    omega = len(support(s) & support(t))
    need(len(fixed(s) & fixed(t)) == 2*a-len(s)+omega,
         'actual node count versus deleted overlap')
    return omega

local_rows = []
for d in range(2, 5):
    pp = list(permutations(range(d)))
    accepted = 0
    for s in pp:
        for t in pp:
            if not braid(s, t):
                continue
            fs = sorted(fixed(s))
            for a in range(1, min(d-1, len(fs))+1):
                for subset in combinations(fs, a):
                    local(s, t, set(subset))
                    accepted += 1
            if mul(s, t) == mul(t, s) and all(k == 1 or k % 2 == 0 for k in lengths(s)):
                need(len(support(s) & support(t)) % 2 == 0,
                     'commuting support intersection is even')
    local_rows.append([d, accepted])

# Equality does not force transpositions.
s4 = cyc(6, (1,5,2,6))
t4 = cyc(6, (1,3,2,4))
need(local(s4, t4, {2,3}) == 0 and lengths(s4) == (1,1,4),
     'four-cycle equality hostile')

# Relation/fixation without the actual reaccess can fail the inequality.
id2 = tuple(range(2)); A, B = {0}, {1}
need(braid(id2, id2) and A <= fixed(id2) and B <= fixed(id2),
     'identity false transport has all inertia predicates')
need(B != {mul(id2,id2)[i] for i in A} and 2*len(A & B) < 3*len(A)-2,
     'actual reaccess is load-bearing')

# Full bounded scalar ledger before and after the actual node-parity sidecar.
scalar_rows = []
n2_boundary = []
for N in (2,3,4,5):
    survived = set()
    for d in range(2,17):
        for a in range(1,d):
            q = d-2*a
            for n1 in range(a+1):
                for n2 in range(a+1):
                    W = a+1-n1-n2
                    if (2*n1 < 3*a-d or 2*n2 < 3*a-d or
                        W < N*max(q,0) or W > N*(d-a)):
                        continue
                    if N == 2 and q == 1:
                        n2_boundary.append((d,a,n1,n2,W))
                    equality = (2*n1 == 3*a-d or 2*n2 == 3*a-d)
                    if equality:
                        even_floor = 2*((max(q,0)+1)//2)
                        if W % 2 or W < N*even_floor:
                            continue
                    need(q in (-1,0) and W == 0, 'uniform final degree and zero-overlap bounds')
                    survived.add((d,a,n1,n2,W))
    expected = set()
    for d in range(2,17):
        if d % 4 == 0:
            k = d//4
            expected.update([(d,2*k,k,k+1,0),(d,2*k,k+1,k,0)])
        if d % 4 == 1:
            k = d//4
            expected.add((d,2*k+1,k+1,k+1,0))
        if d % 4 == 2:
            k = d//4
            expected.add((d,2*k+1,k+1,k+1,0))
    need(survived == expected, 'complete bounded table, without hidden degree filter')
    scalar_rows.append([N,len(survived)])
need((3,1,0,0,2) in n2_boundary, 'N2 scalar boundary is retained before node parity')
ps3 = [p for p in permutations(range(3)) if lengths(p) == (1,2)]
need({len(support(p)&support(q)) for p in ps3 for q in ps3 if mul(p,q)==mul(q,p)} == {2},
     'N2 scalar omega1 has no commuting-transposition node realization')

def passport(d, a, gens, cusps, node_pair, expected_group, expected_n, label):
    G = group(gens)
    need(len(G) == expected_group, label+' exact generated group order')
    need({p[0] for p in G} == set(range(d)), label+' actual transitivity')
    conjugates = {mul(mul(h,gens[0]),inverse(h)) for h in G}
    nn = []
    for s, t in cusps:
        need(s in conjugates and t in conjugates, label+' cusp meridians in common conjugacy class')
        nn.append(local(s,t,fixed(s)))
        need(len(fixed(s)) == a, label+' actual generic count')
    need(nn == expected_n, label+' complete cusp counts')
    s, t = node_pair
    need(s in conjugates and t in conjugates, label+' node meridians in common class')
    omega = node(s,t,a)
    need(omega == 0 and -a+sum(nn) == 1, label+' Euler for every N>=2')
    return dict(d=d,a=a,cusps=nn,omega=omega,group_order=len(G),cycle_type=lengths(gens[0]))

g4 = [cyc(4,(1,2)),cyc(4,(2,3)),cyc(4,(3,4))]
control4 = passport(4,2,g4,[(g4[0],g4[1]),(g4[0],g4[0])],
                    (g4[0],g4[2]),24,[1,2],'degree4 transposition')
g5 = [cyc(5,(1,2)),cyc(5,(2,3)),cyc(5,(3,4))]
need(local(g5[0],g5[1],fixed(g5[0])) == 2, 'degree5 exact cusp count')
need(node(g5[0],g5[2],3) == 0, 'degree5 exact node overlap')
need(-3+2+2 == 1 and len({p[0] for p in group(g5)}) == 4,
     'degree5 local Euler survives while global action is intransitive')
need(3 < 5-1, 'three-transposition edge bound excludes every degree5 action')

g6 = [cyc(6,(1,2,3)),cyc(6,(1,3,4)),cyc(6,(4,5,6))]
control6 = passport(6,3,g6,[(g6[0],g6[1]),(g6[2],cyc(6,(4,6,3)))],
                    (g6[0],g6[2]),360,[2,2],'degree6 A6')
A = fixed(g6[0])
need({mul(g6[1],g6[0])[i] for i in A} != fixed(g6[1]),
     'reversing product without transporting gauge fails')
need(all(len(support(p)) == 3 for p in g6), 'degree6 three-cycle support')

payload = dict(local_universe=local_rows,scalar_universe=scalar_rows,
               N2_scalar_boundary=n2_boundary,degree4=control4,degree6=control6,gates=gates)
digest = hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('Ordinary cusp injection/equality: exhaustive local rows',local_rows)
print('Exact scalar+node-parity tables (N,count):',scalar_rows)
print('Scalar N=2 q=1 boundary retained, then excluded by actual node parity')
print('Abstract transposition passport:',json.dumps(control4,sort_keys=True))
print('Abstract A6 passport:',json.dumps(control6,sort_keys=True))
print('No actual affine-complement representation or Keller realization asserted')
print('Always-active exact gates:',gates)
print('Semantic SHA256:',digest)
print('RESULT: PASS')
