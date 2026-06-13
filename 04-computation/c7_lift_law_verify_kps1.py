#!/usr/bin/env python3
"""
c7_lift_law_verify_kps1.py — kind-pasteur-2026-06-09-S1, BRANCH A final checks

V1  pair-feature identity p(3,4,3) == 2 * p(3,3,2) on all 74 classes n=3..6
    (the rank deficiency found in k2_spectrum_functional_kps1 S2)
V2  clean integer lift law, all 74 classes:
      c7'(K2) == 128 c7 + 192 c6 + 80 c5 + 8 c4
                 + 64 p331 + 48 p332 + 64 p341 + 32 p342
V3  H(D) mod laws by n: census (n, H(D) mod 4) and (n, H(D) mod 8);
    candidate: n odd => H(D) == 1 (mod 4)
V4  H(SC) == 1 (mod 4) on all classes
Output: 05-knowledge/results/c7_lift_law_verify_kps1.out
"""
import sys, itertools, time
sys.path.insert(0, '04-computation')
import numpy as np
from skew_doubling_core_kps1 import all_tournaments, H_count, D_blowup, D_skew, D_scblowup

OUT = open('05-knowledge/results/c7_lift_law_verify_kps1.out', 'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

def iso_classes_fast(n):
    perms = list(itertools.permutations(range(n)))
    seen = set(); reps = []
    for A in all_tournaments(n):
        if A.tobytes() in seen:
            continue
        reps.append(A.copy())
        for p in perms:
            seen.add(np.ascontiguousarray(A[np.ix_(p, p)]).tobytes())
    return reps

def ham_cycles_subset(A, S):
    k = len(S)
    adj = [0] * k
    for i in range(k):
        m = 0
        for j in range(k):
            if A[S[i], S[j]]:
                m |= 1 << j
        adj[i] = m
    dp = [[0] * k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k, 2):
        row = dp[mask]
        for last in range(k):
            c = row[last]
            if not c:
                continue
            avail = adj[last] & ~mask
            while avail:
                b = avail & -avail
                nxt = b.bit_length() - 1
                dp[mask | b][nxt] += c
                avail ^= b
    full = (1 << k) - 1
    return sum(dp[full][l] for l in range(1, k) if adj[l] & 1)

def cycle_sets(A, lengths):
    n = A.shape[0]
    out = {}
    for k in lengths:
        if 3 <= k <= n:
            for S in itertools.combinations(range(n), k):
                c = ham_cycles_subset(A, S)
                if c:
                    out[S] = c
    return out

def pair_features(allc):
    items = list(allc.items())
    pf = {}
    for i in range(len(items)):
        S1, c1 = items[i]
        l1 = len(S1); set1 = set(S1)
        if c1 >= 2:
            pf[(l1, l1, l1)] = pf.get((l1, l1, l1), 0) + c1 * (c1 - 1) // 2
        for j in range(i + 1, len(items)):
            S2, c2 = items[j]
            s = len(set1 & set(S2))
            a, b = sorted((l1, len(S2)))
            pf[(a, b, s)] = pf.get((a, b, s), 0) + c1 * c2
    return pf

t0 = time.time()
w("=== c7_lift_law_verify_kps1 ===")
v1 = v2 = v4 = 0
tot = 0
tabD4, tabD8 = {}, {}
for n in (3, 4, 5, 6):
    for idx, A in enumerate(iso_classes_fast(n)):
        tot += 1
        allc = cycle_sets(A, range(3, n + 1))
        spec = tuple(sum(c for S, c in allc.items() if len(S) == k) for k in (3, 4, 5, 6, 7))
        pf = pair_features(allc)
        Ad = D_blowup(A)[0]
        c7p = sum(ham_cycles_subset(Ad, S)
                  for S in itertools.combinations(range(Ad.shape[0]), 7)) if Ad.shape[0] >= 7 else 0
        p331, p332 = pf.get((3, 3, 1), 0), pf.get((3, 3, 2), 0)
        p341, p342, p343 = pf.get((3, 4, 1), 0), pf.get((3, 4, 2), 0), pf.get((3, 4, 3), 0)
        ok1 = (p343 == 2 * p332)
        pred = (128 * spec[4] + 192 * spec[3] + 80 * spec[2] + 8 * spec[1]
                + 64 * p331 + 48 * p332 + 64 * p341 + 32 * p342)
        ok2 = (c7p == pred)
        v1 += ok1; v2 += ok2
        if not ok1:
            w(f"  V1 FAIL n={n} idx={idx}: p343={p343} 2*p332={2*p332}")
        if not ok2:
            w(f"  V2 FAIL n={n} idx={idx}: c7'={c7p} pred={pred} spec={spec} "
              f"p=({p331},{p332},{p341},{p342},{p343})")
        hD = H_count(D_skew(A)[0])
        hS = H_count(D_scblowup(A)[0])
        v4 += (hS % 4 == 1)
        tabD4.setdefault((n, hD % 4), 0); tabD4[(n, hD % 4)] += 1
        tabD8.setdefault((n, hD % 8), 0); tabD8[(n, hD % 8)] += 1
w(f"V1: p(3,4,3) == 2*p(3,3,2): {v1}/{tot}")
w(f"V2: c7'(K2) == 128c7+192c6+80c5+8c4+64p331+48p332+64p341+32p342: {v2}/{tot}")
w(f"V3: (n, H(D) mod 4) census: {dict(sorted(tabD4.items()))}")
w(f"    (n, H(D) mod 8) census: {dict(sorted(tabD8.items()))}")
w(f"V4: H(SC) == 1 (mod 4): {v4}/{tot}")
w(f"=== done in {time.time()-t0:.1f}s ===")
OUT.close()
