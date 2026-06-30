#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""R-odd first Betti number b_1^- of the arc-flip tournament metagraph  (klein-2026-06-29-S6).

Tasks (owner): (1) does 7 | b_1^-(n) for n>=7 (apex-7 divisibility persist)?  (2) is b_1^-(5)=7
and are the 7 R-odd cycles the 7 Fano lines under a Singer Z_7 (octonion structure)?

DEFINITION. Arc-flip metagraph G_n: vertices = iso classes (A000568(n)), simple edge i~j iff a
single dominance reversal sends a rep of i to j (i!=j). Connected. Complement involution R (T->T^op)
is a graph automorphism (THM-584); it acts on H_1(G_n;Q) = H_1^+ (+) H_1^-. b_1^- = dim H_1^-.

LEFSCHETZ (corrected). trace(R_*|H_1) = 1 - #Rfix_vertices + signed #Rfix_edges. R-fixed vertices =
SC classes. An edge {u,v} is R-fixed setwise iff (a) u,v both SC (orientation kept, +1) OR (b) v=R(u)
(a complement PAIR that is single-flip adjacent; orientation flipped, -1). Case (b) is NONempty in
general (a class one flip from its own complement). Hence
    b_1^- = (beta1 - trace)/2 = (E - V + SC - E_SCSC + E_comppair)/2,
beta1=E-V+1, E_SCSC=#edges both SC, E_comppair=#edges {c, R(c)}. Cross-checked vs a direct eigenspace
computation of R_* on the cycle space (n<=6).  Metagraph built by BFS over arc-flips (reaches n=7).
"""
from __future__ import annotations
import itertools
from math import comb
from collections import deque
import numpy as np

def pairs(n): return [(i, j) for i in range(n) for j in range(i+1, n)]
def perm_tables(n):
    P = pairs(n); idx = {p: k for k, p in enumerate(P)}; T = []
    for perm in itertools.permutations(range(n)):
        row = []
        for (i, j) in P:
            a, b = perm[i], perm[j]
            row.append((idx[(a, b)], False) if a < b else (idx[(b, a)], True))
        T.append(row)
    return T
def make_canon(T):
    def canonical(bits):
        best = None
        for row in T:
            v = 0
            for q, (tgt, inv) in enumerate(row):
                bit = (bits >> tgt) & 1
                if inv: bit ^= 1
                v |= bit << q
            if best is None or v < best: best = v
        return best
    return canonical

def build_metagraph_bfs(n):
    d = comb(n, 2); full = (1 << d) - 1
    canonical = make_canon(perm_tables(n))
    start = canonical(0)
    seen = {start}; q = deque([start]); edges = set(); sigma = {}
    while q:
        c = q.popleft()
        sigma[c] = canonical(c ^ full)
        for a in range(d):
            nb = canonical(c ^ (1 << a))
            if nb != c:
                edges.add((min(c, nb), max(c, nb)))
            if nb not in seen:
                seen.add(nb); q.append(nb)
    classes = sorted(seen); cidx = {c: i for i, c in enumerate(classes)}
    # ensure sigma computed for all (BFS sets it when popped; all popped)
    for c in classes:
        if c not in sigma: sigma[c] = canonical(c ^ full)
    E = [(cidx[u], cidx[v]) for (u, v) in edges]
    sig_v = [cidx[sigma[c]] for c in classes]
    SCset = set(i for i, c in enumerate(classes) if sigma[c] == c)
    return classes, cidx, sig_v, SCset, E

def analyze(n, do_direct=True):
    classes, cidx, sig_v, SCset, E = build_metagraph_bfs(n)
    V = len(classes); Ec = len(E); beta1 = Ec - V + 1
    Eset = set(E)
    E_scsc = sum(1 for (u, v) in E if u in SCset and v in SCset)
    E_comp = sum(1 for (u, v) in E if sig_v[u] == v)         # edges {c, R(c)}
    trace = 1 - len(SCset) + E_scsc - E_comp
    bm = (beta1 - trace) // 2
    chk = None
    if do_direct and V <= 60:
        eidx = {e: k for k, e in enumerate(E)}
        Sg = np.zeros((Ec, Ec)); B = np.zeros((V, Ec))
        for k, (u, v) in enumerate(E):
            su, sv = sig_v[u], sig_v[v]
            tu, tv = (su, sv) if su < sv else (sv, su)
            Sg[eidx[(tu, tv)], k] = 1.0 if su < sv else -1.0
            B[u, k] -= 1.0; B[v, k] += 1.0
        s_ = np.linalg.svd(B, compute_uv=False)
        rank = int(np.sum(s_ > 1e-9));
        _, _, vt = np.linalg.svd(B); ker = vt[rank:].T
        M = ker.T @ Sg @ ker
        chk = int(round(np.sum(np.linalg.eigvals(M).real < 0)))
    return dict(V=V, E=Ec, beta1=beta1, SC=len(SCset), E_scsc=E_scsc, E_comp=E_comp,
                trace=trace, b1m=bm, b1p=beta1-bm, direct=chk)

if __name__ == "__main__":
    print("="*92)
    print(" R-odd first Betti b_1^- of the arc-flip metagraph  (klein-S6; BFS + corrected Lefschetz)")
    print("="*92)
    print(f" {'n':>2} {'V':>5} {'E':>6} {'beta1':>6} {'SC':>4} {'E_SCSC':>7} {'E_comp':>7} "
          f"{'b1-':>6} {'b1- mod 7':>9} {'direct':>7}")
    vals = {}
    for n in (3, 4, 5, 6, 7):
        r = analyze(n, do_direct=(n <= 6))
        vals[n] = r['b1m']
        dc = r['direct'] if r['direct'] is not None else '-'
        print(f" {n:>2} {r['V']:>5} {r['E']:>6} {r['beta1']:>6} {r['SC']:>4} {r['E_scsc']:>7} "
              f"{r['E_comp']:>7} {r['b1m']:>6} {r['b1m']%7:>9} {str(dc):>7}")
    print(f"\n b_1^-(n) = {[vals[n] for n in (3,4,5,6,7)]}  (n=3..7)")
    print(f" b_1^-(5) = {vals[5]} (owner's claim 7: {'CONFIRMED' if vals[5]==7 else 'NO'})")
    print(f" 7 | b_1^-(7)? {vals[7]} mod 7 = {vals[7]%7}  -> {'YES' if vals[7]%7==0 else 'NO'}")
