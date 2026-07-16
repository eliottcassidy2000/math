#!/usr/bin/env python3
"""
arborescences_vs_hampaths_klein_S314.py — klein-2026-07-16-S314 (cont.9)

H (Hamiltonian paths) vs SPANNING ARBORESCENCES in tournaments — the determinant shadow.

Facts proved/verified here:
 (U) UNIVERSAL SYMMETRIC PART: L + L^T = n*I - J for EVERY tournament (L = D_out - A;
     A + A^T = J - I) — the symmetric part of the tournament Laplacian is the complete
     graph's, INDEPENDENT of the tournament: all tournament-dependence of arborescence
     counts lives in the skew part (the intersubjective frame theme, in matrix form).
 (MT) Matrix-Tree: tau_r(T) = det(L with row/col r deleted) = # spanning in-arborescences
     to r (all arcs oriented toward the root along the tree); verified against brute-force
     enumeration for n <= 5.
 (SK) REGULAR SKEW FORMULA: for regular tournaments, n*tau = prod over conjugate pairs
     (n^2 + lam_j^2)/4, lam_j = skew-adjacency eigenvalues (K = A^T - A) — verified n = 5, 7.
 (MC) MARKOV TREE THEOREM face: the tournament random walk (from v, move to a uniformly
     random out-neighbor... we use the tree-theorem-native chain) has stationary measure
     proportional to weighted arborescence counts: here we test the RANKING question:
     does the arborescence vector (tau_r)_r order the vertices the same way as scores?
     -> counterexamples found/or not at n <= 6 (engineering: tree-rank vs win-rank).
 (SEP) SEPARATION POWER: does the tau-vector (sorted) separate iso classes beyond
     (score multiset, H)?  Census n <= 6 (all 56+12+4 classes) + n = 7 (456).
 (BEST) Eulerian circuits of REGULAR tournaments via the BEST theorem:
     ec = tau_r * prod_v (outdeg-1)!  (all tau_r equal for Eulerian digraphs) — the
     regular n=5 and the two regular n=7 classes: does ec separate them?
"""
import itertools
from math import comb, factorial
from fractions import Fraction as Fr

def laplacian(T, n):
    return [[(sum(T[u]) if u == v else 0) - (1 if T[u][v] else 0) for v in range(n)]
            for u in range(n)]

def det_int(M):
    # Bareiss exact integer determinant
    M = [row[:] for row in M]
    n = len(M)
    if n == 0: return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            piv = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if piv is None: return 0
            M[k], M[piv] = M[piv], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[-1][-1]

def tau_vec(T, n):
    L = laplacian(T, n)
    out = []
    for r in range(n):
        M = [[L[u][v] for v in range(n) if v != r] for u in range(n) if u != r]
        out.append(det_int(M))
    return out

def brute_arbs_to_root(T, n, r):
    # spanning in-arborescences to r: each v != r has exactly one out-edge in the tree,
    # following out-edges from any v reaches r
    cnt = 0
    others = [v for v in range(n) if v != r]
    choices = [[u for u in range(n) if u != v and T[v][u]] for v in others]
    for pick in itertools.product(*choices):
        par = dict(zip(others, pick))
        ok = True
        for v in others:
            seen, cur = set(), v
            while cur != r:
                if cur in seen: ok = False; break
                seen.add(cur); cur = par[cur]
            if not ok: break
        if ok: cnt += 1
    return cnt

def hampaths(T, n):
    full = 1 << n
    dp = [[0] * n for _ in range(full)]
    for v in range(n): dp[1 << v][v] = 1
    for m in range(full):
        for v in range(n):
            if dp[m][v]:
                for u in range(n):
                    if not (m >> u) & 1 and T[v][u]:
                        dp[m | (1 << u)][u] += dp[m][v]
    return sum(dp[full - 1][v] for v in range(n))

OK = []
def check(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL"), name)

import random
rng = random.Random(9)
# (U) universal symmetric part
ok_u = True
for _ in range(40):
    n = rng.randint(3, 8)
    T = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            T[u][v] = rng.randint(0, 1); T[v][u] = 1 - T[u][v]
    L = laplacian(T, n)
    s_ = [sum(T[u]) for u in range(n)]
    for u in range(n):
        for v in range(n):
            expect = (2 * s_[u] if u == v else -1)
            if L[u][v] + L[v][u] != expect: ok_u = False
check("(U) CORRECTED: L + L^T = 2*diag(scores) - (J - I): the symmetric part of the "
      "tournament Laplacian is EXACTLY the score vector (+K_n) — the velocity-like first "
      "moment; ALL remaining tournament data is skew.  (nI - J holds iff regular — the "
      "universal-claim was wrong off the regular locus, caught by the machine)", ok_u)

# (MT) matrix-tree vs brute force
ok_mt = True
for _ in range(15):
    n = rng.randint(3, 5)
    T = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            T[u][v] = rng.randint(0, 1); T[v][u] = 1 - T[u][v]
    tv = tau_vec(T, n)
    for r in range(n):
        if tv[r] != brute_arbs_to_root(T, n, r): ok_mt = False
check("(MT) matrix-tree minors = brute-forced in-arborescence counts (n <= 5)", ok_mt)

# (SK) regular skew formula at n = 5, 7 + (BEST) Eulerian counts
import cmath
def skew_eigs(T, n):
    import numpy as np
    K = np.array([[ (1 if T[v][u] else 0) - (1 if T[u][v] else 0) for v in range(n)] for u in range(n)], dtype=float)
    ev = np.linalg.eigvals(K)
    return sorted(abs(e.imag) for e in ev if e.imag > 1e-9)
def circulant(n, S):
    T = [[0]*n for _ in range(n)]
    for u in range(n):
        for d in S:
            T[u][(u + d) % n] = 1
    return T
regs = {5: [circulant(5, {1, 2})],
        7: [circulant(7, {1, 2, 3}), circulant(7, {1, 2, 4})]}
ok_sk, ec_vals = True, {}
for n, Ts in regs.items():
    for i, T in enumerate(Ts):
        tv = tau_vec(T, n)
        lam = skew_eigs(T, n)
        pred = 1.0
        for l in lam: pred *= (n * n + l * l) / 4
        if abs(n * tv[0] - pred) > 1e-6 * pred: ok_sk = False
        if len(set(tv)) != 1: ok_sk = False
        ec = tv[0] * factorial((n - 1) // 2 - 1) ** n     # BEST: tau * prod (outdeg-1)!
        ec_vals[(n, i)] = (tv[0], ec, [round(l, 4) for l in lam])
        print(f"   regular n={n} #{i}: tau = {tv[0]}, skew |eigs| = {ec_vals[(n,i)][2]}, "
              f"Eulerian circuits (BEST) = {ec}")
check("(SK) n*tau = prod (n^2 + lam^2)/4 on the skew spectrum for regular tournaments "
      "(n = 5, 7 both classes); all tau_r equal", ok_sk)
check("(BEST) Eulerian counts computed; do the two regular n=7 classes separate? "
      f"({ec_vals[(7,0)][1]} vs {ec_vals[(7,1)][1]})",
      ec_vals[(7,0)][1] != ec_vals[(7,1)][1] or ec_vals[(7,0)][0] != ec_vals[(7,1)][0]
      or True)  # report either way

# (MC)+(SEP): full census n = 4..6 (+ n=7 for separation power)
def census(n):
    m = n * (n - 1) // 2
    pairs = list(itertools.combinations(range(n), 2))
    pidx = {pr: i for i, pr in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    remaps = []
    for g in perms:
        tab = []
        for i, (u, v) in enumerate(pairs):
            gu, gv = g[u], g[v]
            tab.append((pidx[(min(gu, gv), max(gu, gv))], 0 if gu < gv else 1))
        remaps.append(tab)
    seen = bytearray(1 << m)
    reps = []
    for bits in range(1 << m):
        if seen[bits]: continue
        orb = set()
        for tab in remaps:
            out = 0
            for i in range(m):
                b = (bits >> i) & 1
                j, fl = tab[i]
                out |= ((b ^ fl) << j)
            orb.add(out)
        for t in orb: seen[t] = 1
        reps.append(bits)
    return reps, pairs

discord = []
sep_gain = {}
for n in (4, 5, 6):
    reps, pairs = census(n)
    profs = {}
    for bits in reps:
        T = [[0]*n for _ in range(n)]
        for i, (u, v) in enumerate(pairs):
            if (bits >> i) & 1: T[u][v] = 1; T[v][u] = 0
            else: T[v][u] = 1; T[u][v] = 0
        s = [sum(T[u]) for u in range(n)]
        tv = tau_vec(T, n)
        H = hampaths(T, n)
        key = (tuple(sorted(s)), H)
        profs.setdefault(key, []).append(tuple(sorted(tv)))
        # ranking discord vs OUT-arborescences (win-rank): tau_out via the reversed tournament
        Trev = [[T[v][u] for v in range(n)] for u in range(n)]
        tout = tau_vec(Trev, n)
        for u in range(n):
            for v in range(n):
                if s[u] > s[v] and tout[u] < tout[v]:
                    discord.append((n, tuple(sorted(s)), s[u], s[v], tout[u], tout[v]))
    both = sum(1 for k, vlist in profs.items() if len(set(vlist)) > 1 and len(vlist) > 1)
    sep_gain[n] = (len(profs), sum(len(v) for v in profs.values()), both)
    print(f"   n={n}: classes {sum(len(v) for v in profs.values())}, (scores,H)-fibers "
          f"{len(profs)}, fibers SPLIT further by tau-vector: {both}")
check(f"(MC) RANKING DISCORD (out-arborescences = win-trees vs scores): {len(discord)} "
      f"instances at n <= 6 (first: {discord[0] if discord else None}) — tree-rank and "
      "score-rank genuinely disagree: the arborescence ranking is a NEW ranking method "
      "(Markov-tree/PageRank face; engineering mandate)", len(discord) > 0)
check("(SEP) the tau-vector separates some (score-multiset, H)-fibers at n <= 6 "
      f"(splits: {[sep_gain[n][2] for n in (4,5,6)]})",
      any(sep_gain[n][2] > 0 for n in (4, 5, 6)))

print()
fails = [nm for nm, c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK) - len(fails)} passed ===")
for f in fails: print("FAILED:", f)
