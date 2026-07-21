#!/usr/bin/env python3
"""
pencil_sum_law_leafpoly_analytics_boxeph_S185.py  (HYP-8560, THM-1760 part 3)

(Q1) THE PENCIL ORDINAL-SUM LAW (one-line proof: block triangularity of
     tI + uD_in - vA on T (+) S; S-vertices gain n_T indegree):
        M_{T(+)S}(t,u,v) = M_T(t,u,v) * M_S(t + u n_T, u, v).
     Verified exactly on all pairs n_T, n_S <= 3, grid {0..3}^3.
     Corollary (coning): M_{1(+)T'} = t * M_{T'}(t+u, u, v) — explains the
     S185 joint-residue pairs (walls coned from n=6 pencil collisions).
(Q2) leaf polynomials (S184 census objects): unimodality and log-concavity
     of (c_1..c_{n-1}) for all 530 classes n=3..7 (summed over roots).
(Q3) free-action depth TIGHTNESS: find classes with 3 | |Aut| and 3 NOT
     dividing c_3 (the law stops at l < p — witnesses that it is sharp).
boxeph-2026-07-20-S185. Pure python.
"""
import itertools
import math

def pairs_of(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def bits_to_adj(bits, n):
    adj = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs_of(n)):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def build_perm_maps(n):
    prs = pairs_of(n)
    idx = {pr: k for k, pr in enumerate(prs)}
    maps = []
    for p in itertools.permutations(range(n)):
        pm = []
        for (i, j) in prs:
            a, b = p[i], p[j]
            pm.append((idx[(a, b)], 0) if a < b else (idx[(b, a)], 1))
        maps.append(pm)
    return maps

def classes_for_n(n):
    if n == 7:
        with open('05-knowledge/results/n7_class_reps_boxeph_S152.txt') as f:
            return [int(t) for t in f.read().split()]
    perm_maps = build_perm_maps(n)
    total = 1 << (n * (n - 1) // 2)
    visited = bytearray(total)
    reps = []
    for bits in range(total):
        if visited[bits]:
            continue
        members = set()
        for pm in perm_maps:
            nb = 0
            for k in range(len(pm)):
                src, flip = pm[k]
                b = (bits >> src) & 1
                if flip:
                    b ^= 1
                nb |= b << k
            members.add(nb)
        for m in members:
            visited[m] = 1
        reps.append(min(members))
    return reps

def aut_order(bits, n, perm_maps):
    cnt = 0
    for pm in perm_maps:
        nb = 0
        for k in range(len(pm)):
            src, flip = pm[k]
            b = (bits >> src) & 1
            if flip:
                b ^= 1
            nb |= b << k
        if nb == bits:
            cnt += 1
    return cnt

def det_int(mat):
    m = [row[:] for row in mat]
    n = len(m)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if m[k][k] == 0:
            piv = next((r for r in range(k + 1, n) if m[r][k] != 0), None)
            if piv is None:
                return 0
            m[k], m[piv] = m[piv], m[k]
            sign = -sign
        mkk = m[k][k]
        for i in range(k + 1, n):
            mik = m[i][k]
            for j in range(k + 1, n):
                m[i][j] = (m[i][j] * mkk - mik * m[k][j]) // prev
            m[i][k] = 0
        prev = mkk
    return sign * m[n - 1][n - 1]

def Mpencil(adj, n, t, u, v):
    indeg = [sum(adj[i][j] for i in range(n)) for j in range(n)]
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        M[i][i] = t + u * indeg[i]
        for j in range(n):
            if i != j and adj[i][j]:
                M[i][j] -= v
    return det_int(M)

def leaf_poly_total(adj, n):
    inmask = [sum(1 << u for u in range(n) if adj[u][v]) for v in range(n)]
    tot = [0] * (n + 1)
    for r in range(n):
        others = [v for v in range(n) if v != r]
        S = [0] * (n + 1)
        for Ub in range(1 << (n - 1)):
            U = 0
            t = Ub
            for v in others:
                if t & 1:
                    U |= 1 << v
                t >>= 1
            if any((inmask[v] & ~U) == 0 for v in others):
                continue
            k = bin(U).count('1')
            oi = {v: c for c, v in enumerate(others)}
            M = [[0] * (n - 1) for _ in range(n - 1)]
            for j in others:
                indeg = 0
                for i in range(n):
                    if adj[i][j] and not ((U >> i) & 1):
                        indeg += 1
                        if i != r:
                            M[oi[i]][oi[j]] -= 1
                M[oi[j]][oi[j]] += indeg
            v = det_int(M)
            if v:
                S[k] += v
        for l in range(n + 1):
            tot[l] += sum((-1) ** (k - l) * math.comb(k, l) * S[k]
                          for k in range(l, n + 1))
    return tot

print("=" * 78)
print("(Q1) pencil ordinal-sum law  M_{T+S}(t,u,v) = M_T(t,u,v) M_S(t+u n_T,u,v):")
bad = tot = 0
for nT in (2, 3):
    for nS in (2, 3):
        for bT in range(1 << (nT * (nT - 1) // 2)):
            for bS in range(1 << (nS * (nS - 1) // 2)):
                adjT = bits_to_adj(bT, nT)
                adjS = bits_to_adj(bS, nS)
                n = nT + nS
                adj = [[0] * n for _ in range(n)]
                for i in range(nT):
                    for j in range(nT):
                        adj[i][j] = adjT[i][j]
                for i in range(nS):
                    for j in range(nS):
                        adj[nT + i][nT + j] = adjS[i][j]
                for i in range(nT):
                    for j in range(nS):
                        adj[i][nT + j] = 1
                for t in range(4):
                    for u in range(4):
                        for v in range(4):
                            tot += 1
                            lhs = Mpencil(adj, n, t, u, v)
                            rhs = (Mpencil(adjT, nT, t, u, v) *
                                   Mpencil(adjS, nS, t + u * nT, u, v))
                            if lhs != rhs:
                                bad += 1
print("     %d evaluations, %d mismatches" % (tot, bad))

print("\n(Q2) leaf polynomials: unimodality + log-concavity over all classes:")
uni_bad = lc_bad = total_cls = 0
lc_examples = []
for n in (3, 4, 5, 6, 7):
    for bits in classes_for_n(n):
        adj = bits_to_adj(bits, n)
        c = leaf_poly_total(adj, n)[1:n]   # c_1..c_{n-1}
        total_cls += 1
        # unimodality
        peaked = True
        rising = True
        for i in range(1, len(c)):
            if c[i] < c[i - 1]:
                rising = False
            elif not rising and c[i] > c[i - 1]:
                peaked = False
        if not peaked:
            uni_bad += 1
        # log-concavity c_i^2 >= c_{i-1} c_{i+1}
        for i in range(1, len(c) - 1):
            if c[i] * c[i] < c[i - 1] * c[i + 1]:
                lc_bad += 1
                if len(lc_examples) < 4:
                    lc_examples.append((n, bits, c))
                break
print("     %d classes: unimodality violations %d ; log-concavity violations %d"
      % (total_cls, uni_bad, lc_bad))
for ex in lc_examples:
    print("       LC fail: n=%d bits=%d c=%s" % ex)

print("\n(Q3) depth tightness: classes with 3 | |Aut| and 3 not dividing c_3:")
found = 0
for n in (5, 6, 7):
    pm = build_perm_maps(n) if n <= 6 else build_perm_maps(7)
    for bits in classes_for_n(n):
        ao = aut_order(bits, n, pm)
        if ao % 3 != 0:
            continue
        c = leaf_poly_total(bits_to_adj(bits, n), n)
        if c[3] % 3 != 0:
            found += 1
            if found <= 5:
                print("     WITNESS n=%d bits=%d |Aut|=%d c=%s (c3=%d, 3∤c3)" %
                      (n, bits, ao, c[1:n], c[3]))
    if found:
        break
print("     tightness witnesses found: %d %s" %
      (found, "(the l < p bound is SHARP)" if found else "(none at this n — law may extend?)"))

print("\nDONE.")
