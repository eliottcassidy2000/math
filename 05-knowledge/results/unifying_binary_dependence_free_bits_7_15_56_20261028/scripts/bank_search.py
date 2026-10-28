#!/usr/bin/env python3
"""bank_search.py -- least prime banks realizing the excluded minors of graphic matroids and R10 (klein, 2026-10-28).
A prime p has square-class vector v(p) = squarefree part of p + 1 (as an F2 vector: exponent parities); for a set of
primes the square-class matroid is the binary matroid of these vectors; XOR of vectors = squarefree part of the product,
a (+) b = a b / gcd(a,b)^2.  Among all primes p < LIMIT we keep, for each squarefree class d, the least prime L(d).
A bank realizing a target matroid T can always be replaced by least primes, so the least possible maximum prime of a
T-bank is  min over T-configurations {d_1..d_n} of max_i L(d_i);  the search is exhaustive inside D(N) = {d : L(d) < N}.
Targets: F7 (Fano), F7* = AG(3,2) minus a point, M*(K3,3) = the 3x3 grid (hyperbolic quadric of PG(3,2)),
M*(K5) (10 elements, rank 6, affine), R10 (10 elements, rank 5, affine).  Standard library only."""
import sys, itertools
from math import gcd
from collections import defaultdict
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
out = []
LOG = open('bank_search.log', 'a')
def say(s): out.append(s); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
# ---- sieve and squarefree parts
sieve = bytearray([1]) * (LIMIT + 2); sieve[0] = sieve[1] = 0
for i in range(2, int((LIMIT + 1) ** 0.5) + 1):
    if sieve[i]: sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
primes = [i for i in range(2, LIMIT) if sieve[i]]
from array import array
spf = array('i', [0]) * (LIMIT + 2)
for i in range(2, LIMIT + 2):
    if spf[i] == 0:
        for j in range(i, LIMIT + 2, i):
            if spf[j] == 0: spf[j] = i
def sqf(n):
    r = 1
    while n > 1:
        q = spf[n]; e = 0
        while n % q == 0: n //= q; e += 1
        if e & 1: r *= q
    return r
L = {}   # squarefree class -> least prime
for p in primes:
    d = sqf(p + 1)
    if d not in L: L[d] = p
say('primes below %d: %d; distinct square classes: %d; loop class d=1 at p=%d' % (LIMIT, len(primes), len(L), L[1]))
def X(a, b): g = gcd(a, b); return a * b // (g * g)
def D_upto(N): return {d for d, p in L.items() if p < N and d != 1}
def maxp(S): return max(L[d] for d in S)
# ---- matroid utilities on a list of class labels
def rank_of(S):
    basis = []
    for v in S:
        w = v
        for b in basis:
            w2 = X(w, b)
            if w2 < w: w = w2     # not a proper reduction; use linear algebra over supports instead
        basis.append(v)
    return None
def fvec(d):            # F2 vector as frozenset of prime factors
    s = set(); n = d; q = 2
    while q * q <= n:
        if n % q == 0:
            s.add(q)
            while n % q == 0: n //= q
        q += 1
    if n > 1: s.add(n)
    return frozenset(s)
def rank(S):
    rows = [set(fvec(d)) for d in S]; r = 0
    while rows:
        piv = rows.pop()
        if not piv: continue
        r += 1; x = min(piv)
        rows = [row ^ piv if x in row else row for row in rows]
    return r
def circuits(S):
    """all circuits of the binary matroid on the labels S (|S| <= 10): subsets with XOR 1 that are minimal"""
    S = list(S); n = len(S); deps = []
    for mask in range(1, 1 << n):
        v = 1
        for i in range(n):
            if mask >> i & 1: v = X(v, S[i])
        if v == 1: deps.append(mask)
    depset = set(deps)
    circ = [m for m in deps if not any((m2 != m and m2 & m == m2) for m2 in deps)]
    return [frozenset(S[i] for i in range(n) if m >> i & 1) for m in circ]
def circuit_profile(S):
    from collections import Counter
    return dict(sorted(Counter(len(c) for c in circuits(S)).items()))
def isomorphic(S, T):
    """backtracking isomorphism of binary matroids given by labels; compares circuit systems"""
    S, T = list(S), list(T)
    if len(S) != len(T) or rank(S) != rank(T): return False
    CS, CT = circuits(S), circuits(T)
    if sorted(map(len, CS)) != sorted(map(len, CT)): return False
    CTset = set(CT)
    # degree signature: number of circuits of each size through each element
    def sig(e, C): 
        from collections import Counter
        return tuple(sorted(Counter(len(c) for c in C if e in c).items()))
    sS = {e: sig(e, CS) for e in S}; sT = {e: sig(e, CT) for e in T}
    if sorted(sS.values()) != sorted(sT.values()): return False
    order = sorted(S, key=lambda e: -sum(1 for c in CS if e in c))
    def rec(i, phi):
        if i == len(order): return True
        e = order[i]; used = set(phi.values())
        for f in T:
            if f in used or sS[e] != sT[f]: continue
            phi[e] = f
            ok = True
            for c in CS:
                if all(x in phi for x in c):
                    if frozenset(phi[x] for x in c) not in CTset: ok = False; break
            if ok and rec(i + 1, phi): return True
            del phi[e]
        return False
    return rec(0, {})
# ---- reference matroids (as sets of F2 vectors encoded by integers over "abstract primes" 2,3,5,7,11,13 as coordinates)
coord = [2, 3, 5, 7, 11, 13]
def vec(bits): 
    d = 1
    for i, b in enumerate(bits):
        if b: d *= coord[i]
    return d
REF = {}
REF['F7'] = [vec(b) for b in itertools.product((0, 1), repeat=3) if any(b)]
REF['F7*'] = [vec((1,) + b) for b in itertools.product((0, 1), repeat=3)][:-1]   # AG(3,2) = (1,x,y,z), minus the point (1,1,1,1)
REF['R10'] = [vec(b) for b in itertools.product((0, 1), repeat=5) if sum(b) == 3]  # weight-3 vectors of F2^5
# M*(K5), M*(K3,3): columns = edges; rows = a cycle basis (fundamental cycles)
def cographic(nv, edges):
    import itertools as it
    # spanning tree by BFS
    adj = defaultdict(list)
    for i, (a, b) in enumerate(edges): adj[a].append((b, i)); adj[b].append((a, i))
    parent = {0: (None, None)}; order = [0]
    for u in order:
        for v, i in adj[u]:
            if v not in parent: parent[v] = (u, i); order.append(v)
    tree = {parent[v][1] for v in parent if parent[v][1] is not None}
    def path(u, v):   # tree path edge set between u and v
        pu, pv = [], []
        x = u
        while x != 0: pu.append(parent[x][1]); x = parent[x][0]
        x = v
        while x != 0: pv.append(parent[x][1]); x = parent[x][0]
        return set(pu) ^ set(pv)
    chords = [i for i in range(len(edges)) if i not in tree]
    cycles = [path(*edges[c]) | {c} for c in chords]      # fundamental cycles, one per chord: a basis of the cycle space
    cols = []
    for e in range(len(edges)):
        d = 1
        for j, cyc in enumerate(cycles):
            if e in cyc: d *= coord[j]
        cols.append(d)
    return cols
K5 = [(a, b) for a in range(5) for b in range(a + 1, 5)]
K33 = [(a, b) for a in range(3) for b in range(3, 6)]
REF['M*(K5)'] = cographic(5, K5)
REF['M*(K3,3)'] = cographic(6, K33)
for name, S in REF.items():
    say('reference %-9s elements %2d rank %d circuits by size %s' % (name, len(S), rank(S), circuit_profile(S)))
# ---- searches (exhaustive inside D(N))
def search_F7(N):
    D = sorted(D_upto(N)); Ds = set(D); best = None
    for i, a in enumerate(D):
        for b in D[i+1:]:
            c = X(a, b)
            if c not in Ds or c <= b: continue
            for e in D:
                if e <= b or e in (a, b, c): continue      # a < b are the two least points; any fourth point e > b spans the plane
                pts = [a, b, c, e, X(a, e), X(b, e), X(c, e)]
                if all(x in Ds for x in pts) and len(set(pts)) == 7:
                    m = maxp(pts)
                    if best is None or m < best[0]: best = (m, sorted(pts))
    return best
def search_F7star(N):
    """AG(3,2) minus a point: 8 points t (+) W, W a 3-dim subspace, with 7 of them in D"""
    D = sorted(D_upto(N)); Ds = set(D); best = None
    for i, p0 in enumerate(D):
        for p1 in D[i+1:]:
            w1 = X(p0, p1)
            for p2 in D:
                if p2 <= p1: continue
                w2 = X(p0, p2)
                if w2 == w1: continue
                for p3 in D:
                    if p3 <= p2: continue
                    w3 = X(p0, p3)
                    if w3 in (w1, w2, X(w1, w2)): continue
                    W = [1, w1, w2, w3, X(w1, w2), X(w1, w3), X(w2, w3), X(X(w1, w2), w3)]
                    pts = [X(p0, w) for w in W]
                    if 1 in pts: continue          # p0 in W: the coset is W itself (a Fano plane), not AG(3,2)
                    inD = [x for x in pts if x in Ds]
                    if len(inD) >= 7:
                        for drop in ([None] if len(inD) == 7 else pts):
                            cfg = [x for x in inD if x != drop]
                            if len(cfg) != 7: continue
                            m = maxp(cfg)
                            if (best is None or m < best[0]) and rank(cfg) == 4 and isomorphic(cfg, REF['F7*']): best = (m, sorted(cfg))
    return best
def search_K33(N):
    """3x3 grid: a1,a2,b1,b2 independent, points {a1,a2,a1a2, b1,b2,b1b2, a1b1, a2b2, a1a2b1b2}"""
    D = sorted(D_upto(N)); Ds = set(D); best = None
    lines = []   # 3-circuits as (a,b,c) with a<b<c
    for i, a in enumerate(D):
        for b in D[i+1:]:
            c = X(a, b)
            if c in Ds and c > b: lines.append((a, b, c))
    say('  N=%d: |D|=%d, 3-circuits=%d' % (N, len(D), len(lines)))
    for (a1, a2, a3) in lines:
        for bl in lines:
            if len(set((a1, a2, a3) + bl)) < 6: continue
            for (b1, b2, b3) in itertools.permutations(bl):      # every row-to-row matching
                c1, c2, c3 = X(a1, b1), X(a2, b2), X(a3, b3)
                if c1 in Ds and c2 in Ds and c3 in Ds:
                    pts = [a1, a2, a3, b1, b2, b3, c1, c2, c3]
                    if len(set(pts)) == 9 and rank(pts) == 4:
                        m = maxp(pts)
                        if (best is None or m < best[0]) and isomorphic(pts, REF['M*(K3,3)']): best = (m, sorted(pts))
    return best
def affine_search(N, target, dims, thresholds):
    """exhaustive search for an affine target (all circuits even) of rank r = dims inside D(N):
    grow subspaces from 4-circuits; at full dimension test the affine hyperplane sections."""
    D = sorted(D_upto(N)); Ds = set(D); best = None
    # 4-circuits via XOR-pair collisions
    pairs = defaultdict(list)
    for i, a in enumerate(D):
        for b in D[i+1:]:
            pairs[X(a, b)].append((a, b))
    quads = set()
    for v, lst in pairs.items():
        if len(lst) >= 2:
            for (a, b), (c, d) in itertools.combinations(lst, 2):
                if len({a, b, c, d}) == 4: quads.add(frozenset((a, b, c, d)))
    say('  N=%d: |D|=%d, 4-circuits=%d' % (N, len(D), len(quads)))
    def span_points(basis):
        pts = set()
        for mask in range(1, 1 << len(basis)):
            v = 1
            for i in range(len(basis)):
                if mask >> i & 1: v = X(v, basis[i])
            pts.add(v)
        return pts
    seen_spaces = set(); tested = 0
    def extend(basis, pts_in_D, depth):
        nonlocal best, tested
        if depth == dims:
            key = frozenset(span_points(basis))
            if key in seen_spaces: return
            seen_spaces.add(key)
            P = sorted(key & Ds)
            if len(P) < len(target): return
            # affine hyperplanes: a functional is determined by its values on the basis; elements with value 1
            for vals in itertools.product((0, 1), repeat=dims):
                if not any(vals): continue
                # value of a point = parity of basis coordinates; recover coordinates by solving: brute force over masks
                H = []
                for mask in range(1, 1 << dims):
                    v = 1
                    for i in range(dims):
                        if mask >> i & 1: v = X(v, basis[i])
                    if v in Ds and sum(vals[i] for i in range(dims) if mask >> i & 1) % 2 == 1: H.append(v)
                if len(H) < len(target): continue
                for cfg in itertools.combinations(sorted(H), len(target)):
                    m = maxp(cfg)
                    if best is not None and m >= best[0]: continue
                    tested += 1
                    if rank(cfg) == dims and circuit_profile(cfg) == circuit_profile(target) and isomorphic(cfg, target):
                        best = (m, sorted(cfg))
            return
        cur = span_points(basis)
        for x in D:
            if x in cur or (best is not None and L[x] >= best[0]): continue
            nb = basis + [x]
            pts = span_points(nb) & Ds
            if len(pts) >= thresholds[depth]:
                extend(nb, pts, depth + 1)
    for q in quads:
        a, b, c, d = sorted(q)
        if best is not None and maxp(q) >= best[0]: continue
        extend([a, b, c], span_points([a, b, c]) & Ds, 3)
    say('  configurations tested: %d' % tested)
    return best
say('--- F7 (Fano) ---'); r = search_F7(200); say('  least max prime %d: classes %s primes %s' % (r[0], r[1], [L[d] for d in r[1]]))
say('--- F7* = AG(3,2) minus a point ---'); r = search_F7star(400); say('  least max prime %d: classes %s primes %s' % (r[0], r[1], [L[d] for d in r[1]]))
say('--- M*(K3,3) = 3x3 grid ---'); r = search_K33(6000); say('  least max prime %d: classes %s primes %s' % (r[0], r[1], [L[d] for d in r[1]]))
say('--- R10 = weight-3 vectors of F2^5 (affine, rank 5) ---'); r = affine_search(15000, REF['R10'], 5, {3: 4, 4: 6, 5: 10}); say('  least max prime %d: classes %s primes %s' % (r[0], r[1], [L[d] for d in r[1]]))
def star_search_K5(N):
    """M*(K5): its four-circuits are exactly the five vertex stars S_0..S_4 of K5 with |S_i ∩ S_j| = 1 (the edge ij)."""
    D = sorted(D_upto(N)); Ds = set(D); best = None
    pairs = defaultdict(list)
    for i, a in enumerate(D):
        for b in D[i+1:]:
            pairs[X(a, b)].append((a, b))
    quads = set()
    for v, lst in pairs.items():
        if len(lst) >= 2:
            for (a, b), (c, d) in itertools.combinations(lst, 2):
                if len({a, b, c, d}) == 4: quads.add(frozenset((a, b, c, d)))
    byel = defaultdict(list)
    for q in quads:
        for e in q: byel[e].append(q)
    say('  N=%d: |D|=%d, 4-circuits=%d, max per element=%d' % (N, len(D), len(quads), max((len(v) for v in byel.values()), default=0)))
    tested = 0
    for S0 in quads:
        if best is not None and maxp(S0) >= best[0]: continue
        e01, e02, e03, e04 = sorted(S0)
        for S1 in byel[e01]:
            if S1 == S0 or len(S1 & S0) != 1: continue
            rest1 = sorted(S1 - {e01})           # e12, e13, e14 in some order
            for (e12, e13, e14) in itertools.permutations(rest1):
                for S2 in byel[e02]:
                    if e12 not in S2 or len(S2 & S0) != 1 or len(S2 & S1) != 1: continue
                    rest2 = sorted(S2 - {e02, e12})   # e23, e24
                    for (e23, e24) in itertools.permutations(rest2):
                        S3 = frozenset((e03, e13, e23))
                        for S3q in byel[e03]:
                            if not S3 <= S3q or len(S3q & S0) != 1 or len(S3q & S1) != 1 or len(S3q & S2) != 1: continue
                            e34 = next(iter(S3q - S3))
                            S4 = frozenset((e04, e14, e24, e34))
                            if S4 not in quads: continue
                            cfg = [e01, e02, e03, e04, e12, e13, e14, e23, e24, e34]
                            if len(set(cfg)) != 10: continue
                            m = maxp(cfg)
                            if best is not None and m >= best[0]: continue
                            tested += 1
                            if rank(cfg) == 6 and circuit_profile(cfg) == circuit_profile(REF['M*(K5)']) and isomorphic(cfg, REF['M*(K5)']):
                                best = (m, sorted(cfg)); say('  candidate max prime %d' % m)
    say('  configurations tested: %d' % tested)
    return best
say('--- M*(K5): star search ---')
for N in [int(x) for x in (sys.argv[2:] or ['5400'])]:   # a construction gives a bank with max prime 5393, so the least bank lies in D(5394); exhaustive there
    r = star_search_K5(N)
    if r: say('  least max prime %d: classes %s primes %s' % (r[0], r[1], [L[d] for d in r[1]])); break
    say('  none below %d' % N)
open('bank_search.out', 'w').write('\n'.join(out) + '\n')
