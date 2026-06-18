#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_tourmap_probe_kps-S2-wf.py

Standalone control probe for the switch-parity tournament maps. For each map,
distinguish:
  (A) classes forbidden by the MAP'S ALGEBRA (free - reachable-at-ANY-tau), and
  (B) classes forbidden by LONELINESS (reachable-at-any-tau - reachable-at-LONELY-tau).

Only (B) is a genuine "loneliness constrains the tournament" find. (A) just means
the construction is structurally incapable of that class regardless of LRC.

Unbuffered prints; smaller dense grid for speed.
"""
import sys
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
from collections import defaultdict

def P(*a): print(*a); sys.stdout.flush()

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2): C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def all_optima(S):
    b = F(0)
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v
    return b, sorted(t for t in cand(S) if gmin(S, t) == b)

def ham(adj, V):
    n = len(V); dp = [[0]*n for _ in range(1 << n)]
    for k in range(n): dp[1 << k][k] = 1
    full = (1 << n) - 1
    for mask in range(1 << n):
        for last in range(n):
            c = dp[mask][last]
            if not c: continue
            for nx in range(n):
                if mask & (1 << nx): continue
                if adj[(V[last], V[nx])]: dp[mask | (1 << nx)][nx] += c
    return sum(dp[full][k] for k in range(n))
def scores(adj, V):
    return tuple(sorted(sum(1 for u in V if u != v and adj[(v, u)]) for v in V))
def c3(adj, V):
    c = 0
    for a, b, cc in combinations(V, 3):
        s = sum(1 for (x, y) in [(a, b), (b, cc), (cc, a)] if adj[(x, y)])
        if s in (0, 3): c += 1
    return c
def canon(adj, V):
    n = len(V)
    base = [[1 if (i != j and adj[(V[i], V[j])]) else 0 for j in range(n)] for i in range(n)]
    best = None
    for p in permutations(range(n)):
        m = tuple(tuple(base[p[i]][p[j]] for j in range(n)) for i in range(n))
        if best is None or m < best: best = m
    return best
def is_t(adj, V):
    for i in range(len(V)):
        for j in range(i+1, len(V)):
            if adj[(V[i], V[j])] == adj[(V[j], V[i])]: return False
    return True
def lab(k, V):
    n = len(V); a = {}
    for i in range(n):
        for j in range(n):
            if i != j: a[(i, j)] = bool(k[i][j])
    return (scores(a, V), c3(a, V), ham(a, V))
def free_set(n):
    V = list(range(n)); pr = list(combinations(range(n), 2)); seen = set()
    for bits in range(1 << len(pr)):
        adj = {}
        for idx, (i, j) in enumerate(pr):
            b = (bits >> idx) & 1
            adj[(i, j)] = bool(b); adj[(j, i)] = not b
        seen.add(canon(adj, V))
    return seen
def prim(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1
def ssets(n, ms):
    return [S for S in combinations(range(1, ms+1), n) if prim(S)]
def frac(x):
    r = x - int(x); return r + 1 if r < 0 else r
def sdev(v, t):
    r = frac(v*t); return r-1 if r > F(1,2) else r
def tslope(v, t):
    r = frac(v*t)
    if r == 0 or r == F(1,2): return 0
    return v if r < F(1,2) else -v

# ---- maps under test (the non-trivial / interesting ones) ----
def m3(S, t):  # meeting parity x speed
    V = list(range(len(S))); sp = list(S); adj = {}
    for i, j in combinations(V, 2):
        vi, vj = sp[i], sp[j]
        kstar = ((vi+vj)*t).__floor__()
        cond = ((kstar & 1) == 0) == (vi < vj)
        adj[(i, j)] = cond; adj[(j, i)] = not cond
    return V, adj

STATE_BEATS = {}
for a in range(4):
    for b in range(4):
        if a != b: STATE_BEATS[(a, b)] = (a < b)
STATE_BEATS[(3, 1)] = True; STATE_BEATS[(1, 3)] = False
def _state(v, t):
    w = (v*t).__floor__(); p = sdev(v, t)
    return 2*(w & 1) + (1 if p > 0 else 0), p
def m4(S, t):  # section-winding rps4
    V = list(range(len(S))); sp = list(S)
    st = {i: _state(sp[i], t) for i in V}; adj = {}
    for i, j in combinations(V, 2):
        si, pi = st[i]; sj, pj = st[j]
        if si != sj: iw = STATE_BEATS[(si, sj)]
        else: iw = (abs(pi), sp[i]) > (abs(pj), sp[j])
        adj[(i, j)] = iw; adj[(j, i)] = not iw
    return V, adj

def m5(S, t):  # separation parity
    V = list(range(len(S))); sp = list(S)
    pos = {i: frac(sp[i]*t) for i in V}; adj = {}
    order = sorted(V, key=lambda i: (pos[i], sp[i]))
    rank = {v: r for r, v in enumerate(order)}; n = len(V)
    for i, j in combinations(V, 2):
        ri, rj = rank[i], rank[j]
        between = (rj-ri-1) if ri < rj else ((n-ri-1)+rj)
        par = between & 1
        if (n-2) % 2 == 0: iw = (par == 0) ^ (rank[i] > rank[j])
        else: iw = (par == 0)
        adj[(i, j)] = iw; adj[(j, i)] = not iw
    return V, adj

MAPS = {"M3_meeting_parity_x_speed": m3,
        "M4_section_winding_rps4": m4,
        "M5_separation_parity": m5}

def probe(n, ms, dense_den=36):
    V = list(range(n))
    free = {lab(k, V) for k in free_set(n)}
    Ss = ssets(n, ms)
    grid = sorted({F(k, d) for d in range(2, dense_den+1) for k in range(1, d) if gcd(k, d) == 1})
    P("#"*72)
    P(f"PROBE n={n}, speeds 1..{ms}  (#sets={len(Ss)}, free={len(free)}, grid_den<= {dense_den})")
    for name, fn in MAPS.items():
        opt = set(); anyt = set()
        for S in Ss:
            Mv, opts = all_optima(S)
            if Mv > 0:
                for t in opts:
                    Vv, adj = fn(S, t)
                    if is_t(adj, Vv): opt.add(canon(adj, Vv))
            Cset = set(cand(S)) | set(grid)
            for t in Cset:
                if t <= 0 or t >= 1: continue
                Vv, adj = fn(S, t)
                if is_t(adj, Vv): anyt.add(canon(adj, Vv))
        optL = {lab(k, V) for k in opt}
        anyL = {lab(k, V) for k in anyt}
        forb_map = free - anyL
        forb_lonely = anyL - optL
        P("-"*60)
        P(f"[{name}]")
        P(f"  free={len(free)}  any-tau={len(anyL)}  lonely-tau={len(optL)}")
        P(f"  forbidden by MAP algebra (free - anytau): {len(forb_map)}")
        if forb_map: P(f"     {sorted(forb_map)}")
        P(f"  forbidden by LONELINESS (anytau - lonely): {len(forb_lonely)}")
        if forb_lonely: P(f"     {sorted(forb_lonely)}")

if __name__ == "__main__":
    for n, ms in [(3, 14), (4, 12), (5, 9)]:
        probe(n, ms)
        P("")
