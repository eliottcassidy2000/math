#!/usr/bin/env python3
"""Independent exhaustive small-state diameter controls; no producer imports."""
from collections import deque, defaultdict
from itertools import combinations, product
from functools import lru_cache
from random import Random
import json

checks = systems = pairs = moves = 0
profiles = defaultdict(int)

def ck(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)

@lru_cache(None)
def span(vs):
    S={0}
    for v in vs:
        S |= {x^v for x in tuple(S)}
    return frozenset(S)

def flat(vs): return span(tuple(sorted(set(vs))))
def rank(vs): return len(flat(vs)).bit_length()-1

def g(n): return n*n-n+1 if n else 0
def c(n): return 1 if n == 1 else (n-1)*(n-2)+1

def examine(classes, base=(), label='linear'):
    global systems, pairs, moves
    s=len(classes); br=rank(base)
    def union(ids): return tuple(v for i in ids for v in classes[i])
    def qr(vs): return rank(tuple(base)+tuple(vs))-br
    ck(all(qr(A)>=s-i for i,A in enumerate(classes)), 'staircase')
    tight=[]
    for mask in range(1,1<<s):
        J=[i for i in range(s) if mask>>i&1]
        rr=qr(union(J))
        ck(rr>=len(J), 'weak Hall')
        if rr==len(J):
            ck(J==list(range(min(J),s)), 'tight set is suffix')
            tight.append(min(J))
    tight.sort(); starts=tight+[s]; a=tight[0] if tight else s
    blocks=list(zip(starts,starts[1:])); bound=g(a)+sum(c(j-i) for i,j in blocks)
    ck(bound<=g(s), 'sum bound')
    for i,j in blocks:
        inner=tuple(base)+union(range(j,s)); ir=rank(inner)
        rq=lambda vs: rank(inner+tuple(vs))-ir
        ck(rq(union(range(i,j)))==j-i, 'block rank')
        ck(rq(classes[i])==j-i, 'buffer full rank')
        for mask in range(1,1<<(j-i-1)):
            J=[k for k in range(i+1,j) if mask>>(k-i-1)&1]
            ck(rq(union(J))>len(J), 'block strict Hall')
    prefixbase=tuple(base)+union(range(a,s)); pr=rank(prefixbase)
    for mask in range(1,1<<a):
        J=[i for i in range(a) if mask>>i&1]
        ck(rank(prefixbase+union(J))-pr>len(J),'prefix strict Hall')
    states=[]; flats=[]
    for st in product(*(range(len(A)) for A in classes)):
        vs=tuple(base)+tuple(classes[i][j] for i,j in enumerate(st))
        if rank(vs)==br+s:
            states.append(st); flats.append(flat(vs))
    ck(bool(states), 'nonempty states')
    adj=[[] for st in states]
    for u in range(len(states)):
        for v in range(u):
            diff=sum(x!=y for x,y in zip(states[u],states[v]))
            if diff==1 or (diff==2 and flats[u]==flats[v]):
                adj[u].append(v); adj[v].append(u); moves+=1
    for u in range(len(states)):
        distances={u:0}; todo=deque([u])
        while todo:
            v=todo.popleft()
            for w in adj[v]:
                if w not in distances:
                    distances[w]=distances[v]+1;todo.append(w)
        ck(len(distances)==len(states),'connected allowed-move graph')
        for v in range(u):
            ck(distances[v]<=bound,'all-pairs decomposition diameter')
            if a==0:
                uniform=1 if s<=2 else s*s-3*s+3
                ck(distances[v]<=uniform,'all-pairs tight-total diameter')
            pairs+=1
    systems+=1; profiles[label]+=1

# Exhaustive distinct vector subsets in the last two colors of F_2^3;
# actual elements remain color-labelled, so repeated vector values across
# classes are distinguishable parallel elements. First color is a fixed basis.
small=[tuple(x) for k in (2,3) for x in combinations(range(1,8),k)]
last=[tuple(x) for k in (1,2) for x in combinations(range(1,8),k)]
for B in small:
    for C in last:
        examine(((1,2,4),B,C),label='exhaustive_F2_3')

# Deterministic controls with loops, parallel copies and non-tight prefixes.
rng=Random(19462026)
for dim in (3,4,5):
    for trial in range(40):
        colors=[]
        for need in (3,2,1):
            A=[]
            while rank(A)<need:
                A.append(rng.randrange(1,1<<dim))
            if trial%4==0:A.append(0)
            if trial%5==0:A.append(A[0])
            colors.append(tuple(A))
        examine(tuple(colors),label='loops_parallel_prefix')

# Multiple tight blocks and an ambient unselected direction. Repeated bases
# within each color are distinct realizations.
for block_sizes in ((1,1,1),(2,1),(1,2),(2,2),(3,1),(1,3),(2,1,2)):
    s=sum(block_sizes); vectors=[1<<i for i in range(s)]; starts=[];q=0
    for b in block_sizes:starts.append((q,q+b));q+=b
    colors=[None]*s
    for i,j in starts:
        # A color has the entire span of its own block plus all inner blocks.
        # Each suffix rank matches exactly at the specified boundaries.
        options=tuple(vectors[i:])
        for k in range(i,j):colors[k]=options+(options[0],0)
    examine(tuple(colors),label='chain_blocks')

# Every graceful independent transversal and every upper-tail fiber for
# labels 0..m, m<=6; complete allowed-move graphs and all state pairs.
for m in range(1,7):
    colors=tuple(tuple((1<<x)^(1<<(x+d)) for x in range(m+1-d)) for d in range(1,m+1))
    states=[st for st in product(*(range(len(A)) for A in colors))
            if rank(tuple(colors[i][j] for i,j in enumerate(st)))==m]
    for r in range(1,m+1):
        tails=sorted({st[r:] for st in states})
        for tail in tails:
            base=tuple(colors[r+i][j] for i,j in enumerate(tail))
            examine(colors[:r],base,label='all_graceful_through_6')

# Exact maximum of block sums over all integer compositions; the recurrence
# enumerates all compositions implicitly, without assuming superadditivity.
maximum=[0]
for n in range(1,501):
    maximum.append(max(c(k)+maximum[n-k] for k in range(1,n+1)))
    if n>=3:ck(maximum[n]<=c(n),'all composition block sums')
    ck(max(g(a)+maximum[n-a] for a in range(n+1))<=g(n),
       'all prefix-and-block decompositions')

print(json.dumps({'status':'PASS','checks':checks,'systems':systems,
                  'profiles':dict(sorted(profiles.items())),
                  'all_pairs_tested':pairs,'allowed_undirected_edges':moves,
                  'composition_bounds_through':500,
                  'scope':'Finite independent controls only; no producer imports. Universal conclusions rest on the reviewed proof and cited strict-rank theorem.'},indent=2,sort_keys=True))
