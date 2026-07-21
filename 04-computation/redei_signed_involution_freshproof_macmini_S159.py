#!/usr/bin/env python3
"""
Toward a FRESH involution/signed-sum proof of Redei's theorem (# Ham paths is ODD).  (mac-mini-S159)
====================================================================================================
Engine: "transitive survive, intransitive cancel."  Redei is about ONE tournament's Hamiltonian
paths.  We look for the SIGNED count and its determinant/involution structure.

Objects computed per tournament T (adjacency A, A[u][v]=1 iff u->v):
  h(T)  = # Hamiltonian paths (Redei: always odd)
  R(T)  = sum over Ham paths pi of sgn(pi as a permutation)   [the SIGNED count -- engine-aligned]
  candidates mod 2:  det(A), det(A+I), Pf/det(A-A^T), ... vs h(T) mod 2
  induction:  h(T) = h(T - v) mod 2 ?
Distributions: h mod 4, h mod 8; |R(T)| distribution; is R(T) always a unit?
"""
import sys, itertools
import numpy as np
from itertools import combinations, permutations, product

def out(*a):
    print(*a); sys.stdout.flush()

def all_tournaments(n):
    pairs = list(combinations(range(n), 2))
    for bits in product((0, 1), repeat=len(pairs)):
        A = [[0]*n for _ in range(n)]
        for (i, j), b in zip(pairs, bits):
            if b: A[i][j] = 1          # i->j
            else: A[j][i] = 1          # j->i
        yield A

def ham_count(n, A):
    """# Hamiltonian paths via Held-Karp DP."""
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            c = dp[mask][last]
            if not c: continue
            for nxt in range(n):
                if mask & (1 << nxt): continue
                if A[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def sgn_perm(p):
    """sign of permutation p (one-line, sequence) via inversions."""
    inv = 0; n = len(p)
    for i in range(n):
        for j in range(i+1, n):
            if p[i] > p[j]: inv += 1
    return -1 if inv & 1 else 1

def ham_count_and_signed(n, A):
    """(h, R) with R = sum_{Ham path pi} sgn(pi).  Enumerate perms (small n)."""
    h = 0; R = 0
    for pi in permutations(range(n)):
        ok = all(A[pi[i]][pi[i+1]] for i in range(n-1))
        if ok:
            h += 1
            R += sgn_perm(pi)
    return h, R

def det_mod2(M, n):
    """det over F2 via Gaussian elimination."""
    M = [row[:] for row in M]
    det = 1
    for c in range(n):
        piv = -1
        for r in range(c, n):
            if M[r][c] & 1: piv = r; break
        if piv < 0: return 0
        if piv != c: M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r != c and (M[r][c] & 1):
                M[r] = [(M[r][k] ^ M[c][k]) for k in range(n)]
    return det  # 1 (nonsingular) since we never returned 0

# ===================================================================== part 1
out("=" * 84)
out("PART 1 -- h(T) odd (Redei), and the SIGNED count R(T) = sum_{Ham paths} sgn(pi)")
out("=" * 84)
for n in (3, 4, 5):
    allh, allR = [], []
    for A in all_tournaments(n):
        h, R = ham_count_and_signed(n, A)
        allh.append(h); allR.append(R)
    allh = np.array(allh); allR = np.array(allR)
    out(f"  n={n}: #tournaments={len(allh)} | h always odd: {bool(np.all(allh % 2 == 1))} | "
        f"h range [{allh.min()},{allh.max()}]")
    out(f"        h mod 4 counts: {dict(zip(*np.unique(allh % 4, return_counts=True)))}")
    out(f"        R(T) range [{allR.min()},{allR.max()}] | |R| values: "
        f"{sorted(set(abs(int(r)) for r in allR))} | R always odd: {bool(np.all(allR % 2 == 1))}")
    out(f"        R(T) == +/-1 always? {bool(np.all(np.abs(allR) == 1))}")

# ===================================================================== part 2
out()
out("=" * 84)
out("PART 2 -- mod-2 matrix candidates vs h(T) mod 2 (=1). Which are ALWAYS odd/nonsingular?")
out("=" * 84)
for n in (3, 4, 5):
    cnt = {"det(A)": 0, "det(A+I)": 0, "det(A-A^T)": 0, "det(A+A^T+I)": 0, "tot": 0}
    for A in all_tournaments(n):
        M = np.array(A)
        AI = (M + np.eye(n, dtype=int)) % 2
        S  = (M - M.T) % 2
        SAT = (M + M.T + np.eye(n, dtype=int)) % 2
        cnt["det(A)"]      += det_mod2(M.tolist(), n)
        cnt["det(A+I)"]    += det_mod2(AI.tolist(), n)
        cnt["det(A-A^T)"]  += det_mod2(S.tolist(), n)
        cnt["det(A+A^T+I)"]+= det_mod2(SAT.tolist(), n)
        cnt["tot"] += 1
    frac = {k: f"{cnt[k]}/{cnt['tot']}" for k in cnt if k != "tot"}
    out(f"  n={n}: #nonsingular mod2 over all T: {frac}")

# ===================================================================== part 3
out()
out("=" * 84)
out("PART 3 -- the induction h(T) = h(T - v) mod 2 (Redei descent), all T, all v")
out("=" * 84)
for n in (4, 5):
    ok = True; checked = 0
    for A in all_tournaments(n):
        hT = ham_count(n, A)
        for v in range(n):
            keep = [u for u in range(n) if u != v]
            B = [[A[i][j] for j in keep] for i in keep]
            hTv = ham_count(n-1, B)
            if (hT % 2) != (hTv % 2):
                ok = False
            checked += 1
    out(f"  n={n}: h(T) = h(T-v) mod 2 for all T,v: {ok}  (checked {checked})")

# ===================================================================== part 4
out()
out("=" * 84)
out("PART 4 -- the SIGNED count R(T): closed form? R(T) = sum_{Ham paths} sgn(pi)")
out("=" * 84)
out("  Testing R(T) against: sign only depends on endpoints? relation to # of 3-cycles?")
for n in (3, 4, 5):
    vals = {}
    ncyc3 = {}
    for A in all_tournaments(n):
        h, R = ham_count_and_signed(n, A)
        # count 3-cycles
        c3 = 0
        for a, b, c in combinations(range(n), 3):
            arcs = [(a, b), (b, c), (c, a)]
            def has(u, v): return A[u][v] == 1
            if (has(a,b) and has(b,c) and has(c,a)) or (has(a,c) and has(c,b) and has(b,a)):
                c3 += 1
        vals.setdefault(R, 0); vals[R] += 1
        ncyc3.setdefault((R, c3), 0); ncyc3[(R, c3)] += 1
    out(f"  n={n}: R distribution {dict(sorted(vals.items()))}")
