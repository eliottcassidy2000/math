#!/usr/bin/env python3
"""Strong 'transitive-plus-one' tournaments: u_1 -> ... -> u_m transitive, plus z with z -> u_i iff w_i = 1 (w_1 = 1, w_m = 0).
Hamiltonian paths <-> sets A of u's preceding z (each run increasing), valid iff (A empty or w_maxA = 0) and
(B empty or w_minB = 1).  Exact fractions.  Reports per m: min delta over the family, whether a balanced pair always
involves z, monotonicity of P(u_i before z) in i, and the maximal jump."""
from fractions import Fraction as F
from itertools import product
import sys

def stats(w):
    m = len(w)
    # enumerate valid A as (initial segments) + (crossed beta < alpha with free middle); weight per element membership
    H = 0
    inA = [0] * m          # number of valid A containing u_i
    pairB = {}             # for i<j: number of valid A with j in A and i in B  (u_j before u_i)
    cnt_ji = [[0] * m for _ in range(m)]
    # initial segments A = [0..k-1]
    for k in range(m + 1):
        ok = (k == 0 or w[k - 1] == 0) and (k == m or w[k] == 1)
        if ok:
            H += 1
            for i in range(k): inA[i] += 1
    # crossed: beta < alpha (0-based), w[beta] = 1, w[alpha] = 0; elements < beta in A, > alpha in B, between free
    for beta in range(m):
        if w[beta] != 1: continue
        for alpha in range(beta + 1, m):
            if w[alpha] != 0: continue
            free = alpha - beta - 1
            tot = 1 << free
            H += tot
            for i in range(beta): inA[i] += tot
            inA[alpha] += tot
            for i in range(beta + 1, alpha): inA[i] += tot // 2
            # pairs (i<j) with j in A, i in B
            # membership: <beta: A; beta: B; (beta,alpha): free; alpha: A; >alpha: B
            def memb(t):
                if t < beta: return 'A'
                if t == beta: return 'B'
                if t < alpha: return 'F'
                if t == alpha: return 'A'
                return 'B'
            for i in range(m):
                for j in range(i + 1, m):
                    mi, mj = memb(i), memb(j)
                    if mi == 'A' or mj == 'B': continue
                    c = tot
                    if mi == 'F': c //= 2
                    if mj == 'F': c //= 2
                    cnt_ji[i][j] += c
    return H, inA, cnt_ji

def delta(w):
    m = len(w)
    H, inA, cnt = stats(w)
    best = F(0); bestz = F(0)
    pz = [F(inA[i], H) for i in range(m)]         # P(u_i before z)
    for i in range(m):
        b = min(pz[i], 1 - pz[i]); bestz = max(bestz, b); best = max(best, b)
    for i in range(m):
        for j in range(i + 1, m):
            p = F(cnt[i][j], H)                        # P(u_j before u_i)
            best = max(best, min(p, 1 - p))
    return best, bestz, pz, H

M = int(sys.argv[1]) if len(sys.argv) > 1 else 14
for m in range(2, M + 1):
    mn = None; z_always = True; mono = True; maxjump = F(0); argmin = None
    for mid in product((0, 1), repeat=m - 2):
        w = (1,) + mid + (0,)
        d, dz, pz, H = delta(w)
        if mn is None or d < mn: mn, argmin = d, (w, H)
        if 3 * dz < 1: z_always = False
        if any(pz[i] > pz[i + 1] for i in range(m - 1)): mono = False
        maxjump = max(maxjump, max(abs(pz[i + 1] - pz[i]) for i in range(m - 1)), pz[0], 1 - pz[-1])
    print(f"m={m} (n={m+1}): {2**(m-2)} words; min delta {mn} = {float(mn):.5f} at w={''.join(map(str,argmin[0]))} (H={argmin[1]}); "
          f"balanced pair with z always: {z_always}; P(u_i before z) nondecreasing always: {mono}; max jump {maxjump} = {float(maxjump):.4f}")
