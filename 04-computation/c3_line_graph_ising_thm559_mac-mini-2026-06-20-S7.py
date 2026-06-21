#!/usr/bin/env python3
"""
c3_line_graph_ising_thm559 — mac-mini-2026-06-20-S7  (verifies THM-559)

c3(T) = n(n^2-1)/24 - (1/2) E_score, where E_score = sum_v (s_v - sbar)^2 is an EXACT
frustrated 2-body Ising energy on the arc spins (line graph L(K_n)):
  E_score = C(n,2)/2 + (1/2) sum_{cherries {e,f}@v} eps(v,e)eps(v,f) sigma_e sigma_f,
  coupling = +1 if shared vertex is an EXTREME (min/max) of the cherry, -1 if MIDDLE; zero field.
Ground state (E_score=0) = regular tournament = max c3 = n(n^2-1)/24.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb

def eps(v, e): return 1 if v == e[0] else -1   # e=(low,high); +1 if v is low endpoint

def check(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    cherries = []
    for a in range(len(edges)):
        for b in range(a+1, len(edges)):
            e, f = edges[a], edges[b]; sh = set(e) & set(f)
            if len(sh) == 1:
                v = sh.pop(); cherries.append((a, b, eps(v, e)*eps(v, f)))
    sbar = F(n-1, 2); mism_ising = mism_c3 = 0
    for orient in product([1, -1], repeat=len(edges)):
        sig = list(orient); s = [0]*n
        for a, (i, j) in enumerate(edges):
            if sig[a] == 1: s[i] += 1
            else: s[j] += 1
        E = sum((F(s[v]) - sbar)**2 for v in range(n))
        ising = F(n*(n-1), 2)/2 + sum(F(J, 2)*sig[a]*sig[b] for a, b, J in cherries)
        if E != ising: mism_ising += 1
        c3 = comb(n, 3) - sum(comb(s[v], 2) for v in range(n))
        if F(c3) != F(n*(n*n-1), 24) - F(1, 2)*E: mism_c3 += 1
    pos = sum(1 for *_, J in cherries if J > 0); neg = len(cherries) - pos
    return mism_ising, mism_c3, len(cherries), pos, neg

print("THM-559 verification (exhaustive over all 2^C(n,2) tournaments):")
for n in [4, 5, 6]:
    mi, mc, nc, p, ng = check(n)
    print(f"  n={n}: E_score==line-graph-Ising mismatches={mi}; "
          f"c3==n(n^2-1)/24-E/2 mismatches={mc}; #couplings={nc}=|E(L(K_n))| (+:{p} -:{ng}, ratio 2:1)")
print("max c3 (regular, E_score=0) = n(n^2-1)/24:", {n: n*(n*n-1)//24 for n in [5,7,9]})
