#!/usr/bin/env python3
r"""
lrc_chi_vs_margin_s582b.py    oracle-2026-06-03-S582o (overnight cycle 2)

Cycle 1 pinned: R_m is the unique round/chi=2 regular tournament => the tight regular orbit
is the AP. Question now: does chi of the OPTIMAL-TIME round tournament separate TIGHT
(M=1/n) from LOOSE (M>1/n)? I.e. is chi=2 <=> tight and chi>=3 => strictly lonely?

For n=8 (runner size m=7), sample configs, compute exact M (pinch) and the half-turn
runner tournament at the optimal time (ties resolved both ways), its dichromatic chi, and
whether it is round. Correlate chi with the margin M*n.
"""
from fractions import Fraction as Fr
from functools import reduce
from math import gcd
from itertools import combinations
import random

def circ(r, C): r %= C; return min(r, C - r)
def pinchM(S):
    best = Fr(0); arg = None
    for C in set(a + b for i, a in enumerate(S) for b in S[i+1:]):
        for mm in range(1, C):
            v = Fr(min(circ(x*mm, C) for x in S), C)
            if v > best: best, arg = v, (mm, C)
    return best, arg

def runner_tournament(S, mm, C, sign):
    k = len(S); pos = [(x*mm) % C for x in S]
    adj = [[0]*k for _ in range(k)]; half = Fr(C, 2)
    for a in range(k):
        for b in range(a+1, k):
            diff = (pos[b]-pos[a]) % C
            if diff < half: ab = True
            elif diff > half: ab = False
            else: ab = (sign*S[a] < sign*S[b])
            if ab: adj[a][b] = 1
            else: adj[b][a] = 1
    return adj

def is_acyclic(adj, Sset):
    for i, j, k in combinations(Sset, 3):
        if adj[i][j]+adj[j][k]+adj[k][i] in (0, 3): return False
    return True
def dichromatic(adj, m):
    for k in range(1, m+1):
        color = [-1]*m
        def bt(v):
            if v == m: return True
            for c in range(k):
                color[v] = c
                if is_acyclic(adj, [u for u in range(v+1) if color[u]==c]) and bt(v+1): return True
            color[v] = -1; return False
        if bt(0): return k
    return m
def locally_transitive(adj, m):
    def trans(Sx):
        Sx = list(Sx)
        return sorted(sum(1 for b in Sx if a!=b and adj[a][b]) for a in Sx) == list(range(len(Sx)))
    return all(trans([w for w in range(m) if w!=v and adj[v][w]]) and
               trans([w for w in range(m) if w!=v and adj[w][v]]) for v in range(m))

def chi_at_opt(S):
    M, (mm, C) = pinchM(S); k = len(S)
    chis = set(); rounds = set()
    for sign in (+1, -1):
        adj = runner_tournament(S, mm, C, sign)
        chis.add(dichromatic(adj, k)); rounds.add(locally_transitive(adj, k))
    return M, chis, rounds

def main():
    n = 8; k = n - 1; rnd = random.Random(582)
    print("=" * 74)
    print(f"n={n}: chi of optimal-time round tournament vs loneliness margin M*n")
    print("=" * 74)
    pool = [tuple(range(1, n)),                          # AP (tight)
            (1, 2, 3, 4, 5, 7, 12), (1, 4, 5, 6, 7, 11, 13)]  # tight sporadics
    while len(pool) < 60:
        S = tuple(sorted(rnd.sample(range(1, 40), k)))
        if reduce(gcd, S) == 1: pool.append(S)
    from collections import defaultdict
    byband = defaultdict(list)
    rows = []
    for S in pool:
        M, chis, rounds = chi_at_opt(S)
        mn = float(M) * n
        rows.append((mn, sorted(chis), sorted(rounds), S))
    rows.sort()
    print(f"  {'M*n':>6}  chi(opt)   round?   example")
    for mn, chis, rounds, S in rows[:8] + rows[-6:]:
        tag = "TIGHT" if abs(mn-1.0) < 1e-9 else ""
        print(f"  {mn:6.3f}   {chis}    {rounds}   {S} {tag}")
    # correlation: chi at tight vs loose
    tight = [r for r in rows if abs(r[0]-1.0) < 1e-9]
    loose = [r for r in rows if r[0] > 1.0 + 1e-9]
    def chiset(rs):
        s = set()
        for _, chis, _, _ in rs: s.update(chis)
        return sorted(s)
    print(f"\n  TIGHT (M*n=1): {len(tight)} configs, chi values = {chiset(tight)}")
    print(f"  LOOSE (M*n>1): {len(loose)} configs, chi values = {chiset(loose)}")
    print(f"  loose configs with chi=2 (would refute 'chi=2 <=> tight'): "
          f"{sum(1 for _,chis,_,_ in loose if 2 in chis)}/{len(loose)}")

    print("\n" + "=" * 74)
    print("READING: if TIGHT => chi has 2 and LOOSE => chi has 3 (cleanly), chi separates")
    print("  the boundary. If many loose configs also show chi=2, chi does NOT separate")
    print("  tight from loose (chi=2 is necessary-not-sufficient for tightness).")

if __name__ == "__main__":
    main()
