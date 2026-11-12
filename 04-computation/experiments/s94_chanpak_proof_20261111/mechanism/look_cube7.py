#!/usr/bin/env python3
"""Inspect additive 3-cubes found at n = 7 (arc_cube_7.out): adjacency, scores, odd cycles through the digit arcs,
and the sufficient condition (1) no odd cycle through two digit pairs, (2) cycles through different digit pairs meet."""
import itertools, sys
n = 7
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
def adj_of(m):
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if (m >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    return A
def ham(A):
    N = len(A); dp = [[0]*N for _ in range(1 << N)]
    for v in range(N): dp[1 << v][v] = 1
    for mask in range(1, 1 << N):
        for l in range(N):
            c = dp[mask][l]
            if not c: continue
            for w in range(N):
                if A[l][w] and not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
    return sum(dp[(1 << N) - 1])
def odd_cycles(A):
    N = len(A); out = []
    for L in (3, 5, 7):
        for S in itertools.combinations(range(N), L):
            f = S[0]
            for perm in itertools.permutations(S[1:]):
                cyc = (f,) + perm
                if all(A[cyc[i]][cyc[(i + 1) % L]] for i in range(L)): out.append(cyc)
    return out
def arcs_of(c): return {(c[i], c[(i + 1) % len(c)]) for i in range(len(c))}
def report(m, E):
    A = adj_of(m)
    print("mask", m, "digit arcs", [pairs[e] for e in E])
    for r in A: print("   ", "".join(map(str, r)), " score", sum(r))
    vals = {}
    for s in itertools.product((0, 1), repeat=len(E)):
        mm = m
        for b, e in zip(s, E):
            if b: mm ^= (1 << e)
        vals[s] = ham(adj_of(mm))
    print("   values", vals)
    # condition check over all 8 tournaments
    ok1 = ok2 = True
    for s in itertools.product((0, 1), repeat=len(E)):
        mm = m
        for b, e in zip(s, E):
            if b: mm ^= (1 << e)
        A2 = adj_of(mm); cyc = odd_cycles(A2)
        through = []
        for e in E:
            i, j = pairs[e]
            through.append([c for c in cyc if (i, j) in arcs_of(c) or (j, i) in arcs_of(c)])
        for a in range(len(E)):
            for b in range(a + 1, len(E)):
                sa = set(map(tuple, through[a])); sb = set(map(tuple, through[b]))
                if sa & sb: ok1 = False
                for c in through[a]:
                    for d in through[b]:
                        if not (set(c) & set(d)): ok2 = False
    # hub: vertex on all cycles through any digit arc
    print("   (1) no odd cycle through two digit pairs:", ok1, "  (2) cycles of different digit pairs always meet:", ok2)
    A2 = adj_of(m); cyc = odd_cycles(A2)
    print("   #odd cycles", len(cyc), " 3-cycles:", [c for c in cyc if len(c) == 3])
E = (0, 11, 18)
for m in (50546, 53616):
    report(m, E)
