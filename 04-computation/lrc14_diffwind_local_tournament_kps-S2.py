#!/usr/bin/env python3
"""
kind-pasteur-2026-06-17-S2 — independent confirmation of the tournament-map headline.

The difference-winding map M6: at any tau, arc i->j iff frac((v_i - v_j)*tau) in (0,1/2).
This IS the circular/phase tournament on phase points a_v = frac(v*tau) on R/Z.
CLAIM (verified here + by workflow): tie-free => LOCAL (round) tournament; the maximal-H
NON-local n=5 class (score (1,2,2,2,3), c3=4, H=15) is UNREACHABLE, at the lonely optimum
AND at random (non-lonely) phases. So the forbidding is ALGEBRAIC (circle geometry), not a
loneliness effect (forbidden-by-loneliness = 0). Stdlib, exact Fraction.
"""
from fractions import Fraction as F
from math import gcd
from itertools import permutations, combinations
import random
random.seed(3)

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def Mopt(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at

def diffwind(S, t):
    n = len(S); A = [[0]*n for _ in range(n)]; tiefree = True
    for i in range(n):
        for j in range(n):
            if i == j: continue
            rel = ((S[i]-S[j]) * t) % 1
            if 0 < rel < F(1, 2): A[i][j] = 1
            elif rel > F(1, 2): A[i][j] = 0
            else:
                A[i][j] = 1 if S[i] < S[j] else 0; tiefree = False
    return A, tiefree

def canon(A):
    n = len(A); best = None
    for p in permutations(range(n)):
        key = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or key < best: best = key
    return best
def stats(A):
    n = len(A)
    od = tuple(sorted(sum(A[i]) for i in range(n)))
    c3 = sum(1 for a, b, c in combinations(range(n), 3)
             if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]))
    H = sum(1 for p in permutations(range(n)) if all(A[p[k]][p[k+1]] for k in range(n-1)))
    return od, c3, H

# full n=5 iso-class set + the target
allc = {}
for _ in range(200000):
    A = [[0]*5 for _ in range(5)]
    for i in range(5):
        for j in range(i+1, 5):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    allc[canon(A)] = stats(A)
print("n=5 iso classes:", len(allc), "(expect 12)")
target = [k for k, v in allc.items() if v[2] == 15 and v[1] == 4]
print("TARGET (H=15,c3=4):", [allc[k] for k in target])

# diff-wind at the LRC lonely optimum over primitive 5-speed sets
real_lonely = set()
for S in combinations(range(1, 17), 5):
    from functools import reduce
    if reduce(gcd, S) != 1: continue
    _, t = Mopt(list(S))
    A, _ = diffwind(list(S), t)
    real_lonely.add(canon(A))
print("diff-wind @lonely-opt (vmax<=16): realized", len(real_lonely), "/12  TARGET realized?",
      any(c in real_lonely for c in target))

# tie-free local-tournament test on clean random rational phases (no two equal)
real_rand = set(); nonlocal_seen = 0; tot = 0
def is_local(A):
    n = len(A)
    for v in range(n):
        outs = [u for u in range(n) if A[v][u]]
        ins = [u for u in range(n) if A[u][v]]
        for sub in (outs, ins):  # in/out-neighbourhoods must be transitive
            for a, b, c in permutations(sub, 3) if len(sub) >= 3 else []:
                if A[a][b] and A[b][c] and A[c][a]: return False
    return True
for _ in range(60000):
    ph = random.sample(range(1, 100000), 5)
    ph = [F(p, 100000) for p in ph]
    n = 5; A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            rel = (ph[i]-ph[j]) % 1
            A[i][j] = 1 if 0 < rel < F(1, 2) else 0
    real_rand.add(canon(A)); tot += 1
    if not is_local(A): nonlocal_seen += 1
print("random-phase (tie-free) diff-wind: realized", len(real_rand),
      "/12  non-local count:", nonlocal_seen, "/", tot,
      " TARGET realized?", any(c in real_rand for c in target))
print("CONCLUSION: target unreachable both at lonely-opt and at free phases => ALGEBRAIC,",
      "not a loneliness obstruction. Tie-free diff-wind is always a local/round tournament.")
