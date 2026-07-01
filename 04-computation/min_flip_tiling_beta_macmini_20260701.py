#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S81 -- MINIMUM-FLIP TILING INVARIANT beta(T).

Tiling model: fix a base Hamiltonian path P; a tournament containing P is encoded by its m=C(n-1,2) TILE
bits (the non-path arcs), all-forward tiling = transitive. To express a tournament T we flip the tiles that
point backward relative to P. Minimizing over the choice of base Ham path:

    beta(T) := min over Hamiltonian paths P of T  of  #{backward tiles rel. P}
             =  C(n,2) - max over Ham-path orders sigma of  #forward-arcs(sigma)

= the MINIMUM number of tile-flips to express T's iso class. Covering radius R(n) = max_iso-class beta.
Also compute the unconstrained feedback version:  beta_order(T) = C(n,2) - max over ALL orders of #forward
= the MINIMUM FEEDBACK ARC SET (MFAS). beta_order <= beta (Ham-path constraint costs >= 0).

Per iso class we report beta, beta_order, #3-cycles c3, score sequence, |Aut|, H(#Ham paths).
"""
import itertools
from math import comb
from collections import defaultdict

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        A = [[0]*n for _ in range(n)]
        for (i, j), b in zip(pairs, bits):
            if b: A[i][j] = 1          # i -> j
            else: A[j][i] = 1          # j -> i
        yield A

def canon(A, n):
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or key < best: best = key
    return best

def forward_arcs(A, order):
    n = len(order); f = 0
    for a in range(n):
        for b in range(a+1, n):
            if A[order[a]][order[b]]: f += 1
    return f

def is_ham_path(A, order):
    return all(A[order[k]][order[k+1]] for k in range(len(order)-1))

def invariants(A, n):
    C2 = comb(n, 2)
    best_path = -1; best_any = -1; H = 0
    for order in itertools.permutations(range(n)):
        f = forward_arcs(A, order)
        if f > best_any: best_any = f
        if is_ham_path(A, order):
            H += 1
            if f > best_path: best_path = f
    beta = C2 - best_path
    beta_order = C2 - best_any
    # 3-cycles
    c3 = 0
    for i, j, k in itertools.combinations(range(n), 3):
        # cyclic iff not transitive on the triple
        s = A[i][j] + A[j][k] + A[k][i]
        if s == 0 or s == 3: c3 += 1
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    return beta, beta_order, c3, scores, H

def aut_size(A, n):
    cnt = 0
    for p in itertools.permutations(range(n)):
        if all(A[i][j] == A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j): cnt += 1
    return cnt

print(f"{'n':>2} {'m':>3} {'#iso':>5} {'R=maxβ':>7} {'R_ord':>6} {'β-distribution (count by value)':>32}")
data = {}
for n in range(3, 7):
    m = comb(n-1, 2)
    seen = {}
    for A in all_tournaments(n):
        c = canon(A, n)
        if c in seen: continue
        seen[c] = A
    # compute invariants per iso rep
    betas = []; recs = []
    for c, A in seen.items():
        beta, beta_order, c3, scores, H = invariants(A, n)
        aut = aut_size(A, n)
        betas.append(beta)
        recs.append((beta, beta_order, c3, scores, H, aut))
    R = max(betas); R_ord = max(b for _, b, *_ in recs)
    dist = defaultdict(int)
    for b in betas: dist[b] += 1
    diststr = " ".join(f"{k}:{dist[k]}" for k in sorted(dist))
    data[n] = (m, len(seen), R, R_ord, dict(dist), recs)
    print(f"{n:>2} {m:>3} {len(seen):>5} {R:>7} {R_ord:>6}   {diststr}")


print("\nR(n) = covering radius (max min-flips over iso classes):", [data[n][2] for n in sorted(data)])
print("R_order(n) = max MFAS:                                   ", [data[n][3] for n in sorted(data)])
print("m = C(n-1,2):                                            ", [data[n][0] for n in sorted(data)])

# detail for n<=6: the extremal (max-beta) tournaments and their signature
for n in [4, 5, 6]:
    if n not in data: continue
    m, niso, R, R_ord, dist, recs = data[n]
    print(f"\n--- n={n}: max-beta (least-compressible) iso classes (beta={R}) ---")
    for beta, beta_order, c3, scores, H, aut in sorted(recs, reverse=True)[:6]:
        print(f"   beta={beta} beta_order(MFAS)={beta_order} c3={c3} scores={scores} H(#hampaths)={H} |Aut|={aut}")
