#!/usr/bin/env python3
"""
excess8_prediction_klein.py  --  klein-2026-07-01-S84

Compute the HYP-3817 PREDICTED flip-rank excess(8) = #{iso classes C on 8 vertices : C self-complementary
AND |Aut(C)| > 8}.  (excess = 0,0,0,1,3 for n=3..7; conjecture excess(n) = #{SC & |Aut|>n}.)

On 8 points, |Aut|>8 with |Aut| ODD (Moon) forces |Aut| in {9,15,21} (11,13,17,19 need a >8-cycle;
27,45,... need >8 points for a faithful action).  Each such class has:
  - |Aut| in {9,21}: an order-3 automorphism of type 3+3+1+1  => appears in fix(sigma1), sigma1=(012)(345)(6)(7)
  - |Aut| = 15 (C3xC5, orbits 3+5): an order-5 automorphism (5-cycle) => appears in fix(tau5), tau5=(01234)(5)(6)(7)
fix(sigma1) and fix(tau5) together capture EVERY |Aut|>8 class.  Enumerate both (2^10 + 2^8), dedup by a
score-refined canonical form, compute exact |Aut| and SC, count distinct SC classes with |Aut|>8.

Self-check: the same machinery on n=7 fix((012)(345)(6)) must reproduce S83 (|Aut| {9:4,21:1}, SC {9:2,21:1}
=> #{SC & |Aut|>7}=3=excess(7)).
"""
from itertools import permutations, product
from collections import Counter

def scores(adj, n): return [sum(adj[i]) for i in range(n)]

def canon(adj, n):
    sc = scores(adj, n)
    order = sorted(range(n), key=lambda v: sc[v])
    so = [sc[v] for v in order]
    blocks = []; i = 0
    while i < n:
        j = i
        while j < n and so[j] == so[i]: j += 1
        blocks.append(order[i:j]); i = j
    best = None
    for combo in product(*[permutations(b) for b in blocks]):
        relabel = [v for b in combo for v in b]     # relabel[newpos] = oldvertex
        ser = tuple(adj[relabel[a]][relabel[b]] for a in range(n) for b in range(n))
        if best is None or ser < best: best = ser
    return best

def aut_order(adj, n):
    sc = scores(adj, n)
    groups = {}
    for v, s in enumerate(sc): groups.setdefault(s, []).append(v)
    gl = [groups[k] for k in sorted(groups)]
    cnt = 0
    for combo in product(*[permutations(g) for g in gl]):
        p = [0]*n
        for g, pm in zip(gl, combo):
            for o, nw in zip(g, pm): p[o] = nw
        if all(adj[i][j] == adj[p[i]][p[j]] for i in range(n) for j in range(n)): cnt += 1
    return cnt

def is_SC(adj, n):
    sc = scores(adj, n)
    bysc = {}
    for v, s in enumerate(sc): bysc.setdefault(s, []).append(v)
    for s in list(bysc):
        t = n-1-s
        if t not in bysc or len(bysc[t]) != len(bysc[s]): return False
    scs = sorted(bysc)
    for combo in product(*[permutations(bysc[n-1-s]) for s in scs]):
        p = [0]*n
        for s, perm in zip(scs, combo):
            for o, nw in zip(bysc[s], perm): p[o] = nw
        if all(adj[i][j] == adj[p[j]][p[i]] for i in range(n) for j in range(n)): return True
    return False

def orbit_reps_of_pairs(n, sigma):
    """Representative unordered pairs, one per sigma-orbit on unordered pairs."""
    seen = set(); reps = []
    for i in range(n):
        for j in range(i+1, n):
            if (i, j) in seen: continue
            orb = set(); a, b = i, j
            while True:
                u, w = (a, b) if a < b else (b, a)
                if (u, w) in orb: break
                orb.add((u, w)); a, b = sigma[a], sigma[b]
            seen |= orb
            reps.append((i, j))
    return reps

def fix_sigma_classes(n, sigma):
    """All sigma-invariant tournaments, deduped to distinct classes (score-refined canon)."""
    reps = orbit_reps_of_pairs(n, sigma)
    classes = {}
    for bits in product((0, 1), repeat=len(reps)):
        adj = [[0]*n for _ in range(n)]
        for k, (i0, j0) in enumerate(reps):
            a, b = (i0, j0) if bits[k] else (j0, i0)          # orientation of the representative
            aa, bb = a, b                                      # propagate around the sigma-orbit
            while True:
                adj[aa][bb] = 1
                aa, bb = sigma[aa], sigma[bb]
                if (aa, bb) == (a, b): break
        k = canon(adj, n)
        if k not in classes: classes[k] = adj
    return list(classes.values()), len(reps)

def analyze(n, sigma, label):
    reps_classes, norb = fix_sigma_classes(n, sigma)
    dist = Counter(); scdist = Counter(); overN = 0; over_sc = 0; over_details = []
    for adj in reps_classes:
        a = aut_order(adj, n)
        dist[a] += 1
        sc = is_SC(adj, n)
        if sc: scdist[a] += 1
        if a > n:
            overN += 1
            if sc: over_sc += 1
            over_details.append((a, sc, tuple(sorted(scores(adj, n)))))
    print(f"  {label}: fix has 2^{norb} tournaments -> {len(reps_classes)} distinct classes")
    print(f"    |Aut| dist: {dict(sorted(dist.items()))};  SC among them: {dict(sorted(scdist.items()))}")
    print(f"    #{{|Aut|>{n}}} = {overN};  #{{|Aut|>{n} AND SC}} = {over_sc}")
    return reps_classes, over_details

if __name__ == "__main__":
    print("="*78)
    print("SELF-CHECK n=7 (must reproduce S83: |Aut| {9:4,21:1}, SC {9:2,21:1}, excess=3)")
    print("="*78)
    analyze(7, [1,2,0,4,5,3,6], "n=7 fix((012)(345)(6))")

    print("\n" + "="*78)
    print("PREDICT excess(8) = #{SC classes on 8 vtx with |Aut|>8}")
    print("="*78)
    c1, d1 = analyze(8, [1,2,0,4,5,3,6,7], "n=8 fix((012)(345)(6)(7))  [captures |Aut| in {9,21}]")
    c2, d2 = analyze(8, [1,2,3,4,0,5,6,7], "n=8 fix((01234)(5)(6)(7))  [captures |Aut| = 15]")

    # union the distinct |Aut|>8 classes across both enumerations (dedup by canon)
    over = {}
    for adj in c1 + c2:
        a = aut_order(adj, 8)
        if a > 8:
            k = canon(adj, 8)
            if k not in over: over[k] = (a, is_SC(adj, 8), tuple(sorted(scores(adj, 8))))
    total_over = len(over)
    total_over_sc = sum(1 for v in over.values() if v[1])
    print("\n" + "-"*78)
    print(f"UNION (dedup): distinct classes with |Aut|>8 = {total_over}")
    adist = Counter(v[0] for v in over.values()); ascd = Counter(v[0] for v in over.values() if v[1])
    print(f"  |Aut| dist: {dict(sorted(adist.items()))};  SC among them: {dict(sorted(ascd.items()))}")
    for k, (a, sc, ss) in sorted(over.items(), key=lambda kv: (-kv[1][0], kv[1][2])):
        print(f"    |Aut|={a:2d}  SC={sc}  scores={ss}")
    print("\n" + "="*78)
    print(f"PREDICTED excess(8) = #{{SC & |Aut|>8}} = {total_over_sc}")
    print(f"  (excess sequence would be 0,0,0,1,3,{total_over_sc} for n=3..8; verifiable once rho(8) is computed)")
    print("="*78)
