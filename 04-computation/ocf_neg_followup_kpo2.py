#!/usr/bin/env python3
"""
ocf_neg_followup_kpo2.py -- THREAD C follow-up (HYP-2380), kind-pasteur-2026-06-10-S2
======================================================================================
Boundary probes after the main census (ocf_negative_eval_kpo2.py):

 Q1. n=7, larger sample: what alpha-vectors give I(-2) > 0?  Is I(-2)=+3
     realized (predicted by the block construction strong4 (+) C3)?
     What alpha-vectors give I(-1) = +1?  Is I(-1) <= 1 on the sample?
 Q2. Explicit constructions:
       (a) strong-4-tournament (+) C3, block-transitive, n=7:
           predict alpha = [1,3,2], I(-2) = 3, I(-1) = 0, H = 15.
       (b) C3 (+) C3 (+) singleton, n=7: predict alpha = [1,2,1], I(-2)=1, H=9.
       (c) k disjoint triangles, block-transitive, n = 3k (k=2,3):
           predict I(x) = (1+x)^k, H = 3^k, I(-2) = (-1)^k, I(-1) = 0.
           k=3 (n=9) is the EXPLICIT counterexample to CAND-H:
           I(-2) = -1 != H - 4*alpha_1 = 27 - 12 = 15  (alpha_3 = 1 > 0).
 Q3. n=8 random sample: I(2)=H spot-verification; I(-2) = H - 4*alpha_1 exact
     (alpha_3 = 0 at n=8: three disjoint odd cycles need >= 9 vertices);
     sign distribution of I(-2); max I(-1); positive I(-2) values.

Same fresh cycle machinery as the main script (MISTAKE-001/023/054 guards).
Pure integer arithmetic.
"""

import sys, itertools, random

# ---- core machinery (duplicated from ocf_negative_eval_kpo2.py, fresh code) ----

def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def adj_from_mask(n, mask, plist):
    adj = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(plist):
        if (mask >> b) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def directed_odd_cycles(adj):
    n = len(adj)
    out = []
    for k in range(3, n + 1, 2):
        for S in itertools.combinations(range(n), k):
            v0 = S[0]
            rest = S[1:]
            for perm in itertools.permutations(rest):
                if not adj[v0][perm[0]]:
                    continue
                ok = True
                for a in range(k - 2):
                    if not adj[perm[a]][perm[a + 1]]:
                        ok = False
                        break
                if ok and adj[perm[-1]][v0]:
                    m = (1 << v0)
                    for v in perm:
                        m |= (1 << v)
                    out.append((m, k))
    return out

def alpha_vector(cycles, n):
    masks = [m for m, _ in cycles]
    alpha = [1, len(masks)]
    kmax = n // 3
    if kmax >= 2:
        cnt = {}
        L = len(masks)
        def dfs(start, used, depth):
            if depth >= 2:
                cnt[depth] = cnt.get(depth, 0) + 1
            if depth == kmax:
                return
            for t in range(start, L):
                if not (masks[t] & used):
                    dfs(t + 1, used | masks[t], depth + 1)
        for s in range(L):
            dfs(s + 1, masks[s], 1)
        for k in range(2, kmax + 1):
            alpha.append(cnt.get(k, 0))
    while alpha and alpha[-1] == 0:
        alpha.pop()
    return alpha

def I_eval(alpha, x):
    v = 0
    for k in reversed(range(len(alpha))):
        v = v * x + alpha[k]
    return v

def held_karp_H(adj):
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if not c:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[full])

def block_transitive(blocks):
    """blocks = list of adjacency matrices; arcs between blocks: earlier -> later."""
    sizes = [len(b) for b in blocks]
    n = sum(sizes)
    adj = [[0] * n for _ in range(n)]
    off = []
    o = 0
    for s in sizes:
        off.append(o)
        o += s
    for bi, B in enumerate(blocks):
        for i in range(len(B)):
            for j in range(len(B)):
                adj[off[bi] + i][off[bi] + j] = B[i][j]
    for bi in range(len(blocks)):
        for bj in range(bi + 1, len(blocks)):
            for i in range(sizes[bi]):
                for j in range(sizes[bj]):
                    adj[off[bi] + i][off[bj] + j] = 1
    return adj

C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
# the strong 4-tournament: 0->1->2->3->0 with 0->2, 1->3 (scores 2,2,1,1)
S4 = [[0, 1, 1, 0],
      [0, 0, 1, 1],
      [0, 0, 0, 1],
      [1, 0, 0, 0]]
SINGLE = [[0]]

def analyze_print(name, adj, expect=None):
    n = len(adj)
    cyc = directed_odd_cycles(adj)
    alpha = alpha_vector(cyc, n)
    H = held_karp_H(adj)
    bylen = {}
    for _, k in cyc:
        bylen[k] = bylen.get(k, 0) + 1
    a1 = alpha[1] if len(alpha) > 1 else 0
    ok = "OK" if I_eval(alpha, 2) == H else "*** OCF FAIL ***"
    print(f"{name:<28} n={n} bylen={bylen} alpha={alpha} H={H} "
          f"I(2)={I_eval(alpha,2)}[{ok}] I(-2)={I_eval(alpha,-2)} "
          f"I(-1)={I_eval(alpha,-1)} I(1)={I_eval(alpha,1)} "
          f"H-4a1={H-4*a1}")
    if expect:
        print(f"  expected: {expect}  ->  "
              f"{'MATCH' if expect == (alpha, I_eval(alpha,-2)) else 'CHECK BY HAND'}")
    return alpha, H

def main():
    rng = random.Random(46064202)

    print("=" * 78)
    print("Q2: explicit constructions")
    print("=" * 78)
    analyze_print("strong4 (+) C3   [n=7]", block_transitive([S4, C3]),
                  ([1, 3, 2], 3))
    analyze_print("C3 (+) strong4   [n=7]", block_transitive([C3, S4]),
                  ([1, 3, 2], 3))
    analyze_print("C3 (+) C3 (+) pt [n=7]", block_transitive([C3, C3, SINGLE]),
                  ([1, 2, 1], 1))
    analyze_print("C3 (+) C3        [n=6]", block_transitive([C3, C3]),
                  ([1, 2, 1], 1))
    analyze_print("C3 x3 blocks     [n=9]", block_transitive([C3, C3, C3]),
                  ([1, 3, 3, 1], -1))
    print()
    print("CAND-H explicit n=9 counterexample: for C3(+)C3(+)C3,")
    print("  I(-2) = -1 but H - 4*alpha_1 = 27 - 12 = 15.  (alpha_3 = 1 > 0)")

    print()
    print("=" * 78)
    print("Q1: n=7 larger sample (5000) -- positive I(-2) and I(-1)=1 anatomy")
    print("=" * 78)
    n = 7
    plist = pairs_of(n)
    m = len(plist)
    pos_alpha = {}
    im1_eq1_alpha = {}
    max_im1 = None
    max_im2 = None
    sgn = [0, 0, 0]
    NS = 5000
    for _ in range(NS):
        mask = rng.randrange(1 << m)
        adj = adj_from_mask(n, mask, plist)
        cyc = directed_odd_cycles(adj)
        alpha = alpha_vector(cyc, n)
        Im2 = I_eval(alpha, -2)
        Im1 = I_eval(alpha, -1)
        sgn[0 if Im2 < 0 else (1 if Im2 == 0 else 2)] += 1
        if Im2 > 0:
            key = (tuple(alpha), Im2)
            pos_alpha[key] = pos_alpha.get(key, 0) + 1
        if Im1 == 1:
            key = tuple(alpha)
            im1_eq1_alpha[key] = im1_eq1_alpha.get(key, 0) + 1
        if max_im1 is None or Im1 > max_im1[0]:
            max_im1 = (Im1, tuple(alpha), mask)
        if max_im2 is None or Im2 > max_im2[0]:
            max_im2 = (Im2, tuple(alpha), mask)
    print(f"n=7 sample {NS}: sign I(-2): neg={sgn[0]} zero={sgn[1]} pos={sgn[2]}")
    print(f"positive-I(-2) cases by (alpha, I(-2)):")
    for k, v in sorted(pos_alpha.items()):
        print(f"   alpha={list(k[0])} I(-2)={k[1]} : {v} samples")
    print(f"I(-1)=+1 cases by alpha: "
          f"{ {tuple(k): v for k, v in sorted(im1_eq1_alpha.items())} }")
    print(f"max I(-1) seen: {max_im1}")
    print(f"max I(-2) seen: {max_im2}")

    print()
    print("=" * 78)
    print("Q3: n=8 random sample (400)")
    print("=" * 78)
    n = 8
    plist = pairs_of(n)
    m = len(plist)
    sgn = [0, 0, 0]
    ok = 0
    mirror_ok = True
    disc_ok = True
    max_im1 = None
    max_im2 = None
    pos_alpha = {}
    im1_dist = {}
    NS = 400
    for _ in range(NS):
        mask = rng.randrange(1 << m)
        adj = adj_from_mask(n, mask, plist)
        cyc = directed_odd_cycles(adj)
        alpha = alpha_vector(cyc, n)
        H = held_karp_H(adj)
        a1 = alpha[1] if len(alpha) > 1 else 0
        a2 = alpha[2] if len(alpha) > 2 else 0
        a3 = alpha[3] if len(alpha) > 3 else 0
        if a3 != 0:
            print(f"  *** alpha_3 != 0 at n=8?! mask={mask} alpha={alpha}")
        if I_eval(alpha, 2) == H:
            ok += 1
        else:
            print(f"  *** OCF FAIL n=8 mask={mask}")
        Im2 = I_eval(alpha, -2)
        Im1 = I_eval(alpha, -1)
        if Im2 != H - 4 * a1:
            mirror_ok = False
            print(f"  *** mirror FAIL n=8 mask={mask}")
        if a1 * a1 < 4 * a2:
            disc_ok = False
            print(f"  *** DISC FAIL n=8 mask={mask} alpha={alpha}")
        sgn[0 if Im2 < 0 else (1 if Im2 == 0 else 2)] += 1
        if Im2 > 0:
            key = (tuple(alpha), Im2)
            pos_alpha[key] = pos_alpha.get(key, 0) + 1
        im1_dist[Im1] = im1_dist.get(Im1, 0) + 1
        if max_im1 is None or Im1 > max_im1[0]:
            max_im1 = (Im1, tuple(alpha))
        if max_im2 is None or Im2 > max_im2[0]:
            max_im2 = (Im2, tuple(alpha))
    print(f"n=8 sample {NS}: I(2)=H: {ok}/{NS} | exact mirror I(-2)=H-4a1: "
          f"{'ALL OK' if mirror_ok else 'FAIL'} | disc>=0: "
          f"{'ALL OK' if disc_ok else 'FAIL'}")
    print(f"  sign I(-2): neg={sgn[0]} zero={sgn[1]} pos={sgn[2]}")
    print(f"  positive-I(-2) cases: "
          f"{ {(tuple(k[0]), k[1]): v for k, v in sorted(pos_alpha.items())} }")
    print(f"  max I(-1): {max_im1} | max I(-2): {max_im2}")
    print(f"  I(-1) > 1 occurrences: "
          f"{sum(v for k, v in im1_dist.items() if k > 1)}")
    print()
    print("done.")

if __name__ == "__main__":
    main()
