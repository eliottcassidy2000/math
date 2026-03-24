#!/usr/bin/env python3
"""
hoffman_chromatic_s314.py — Hoffman bound and Lie rank connection to χ = n-1
opus-2026-03-24-S314

HOFFMAN BOUND: χ ≥ 1 + λ_max / |λ_min|
If tight: n-1 = 1 + λ_max/|λ_min| → λ_max/|λ_min| = n-2

LIE ALGEBRA RANK: A_{n-1} has rank n-1 = χ.
Is χ = rank(A_{n-1}) a coincidence or deep?

STRIPS: δ_{n-2} has n-2 strips. χ = n-1 = #strips + 1.
"""

import sys, subprocess, numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def build_metagraph(n):
    m_total = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1, n)]
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    lines_raw = [l.strip() for l in (r.stdout or '').split('\n')
                 if len(l.strip()) == m_total and all(c in '01' for c in l.strip())]
    def b2a(bits):
        a = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if bits & (1 << (m_total-1-k)): a[i][j] = 1
            else: a[j][i] = 1
        return a
    def ss(a): return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))
    def c3(a):
        c = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if a[i][j] and a[j][k] and a[k][i]: c += 1
                    if a[i][k] and a[k][j] and a[j][i]: c += 1
        return c
    def H(a):
        dp = {}
        for v in range(n): dp[(1<<v,v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S&(1<<v)): continue
                val = dp.get((S,v), 0)
                if val == 0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if a[v][u]: dp[(S|(1<<u),u)] = dp.get((S|(1<<u),u),0) + val
        return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
    def cp(a):
        sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
        sg = defaultdict(list)
        for v in range(n): sg[sc[v]].append(v)
        gs = [sg[s] for s in sorted(sg.keys())]
        best = None; cnt = 0
        def gp(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for rest in gp(gs2[1:]): yield list(p)+rest
        for p in gp(gs):
            f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or f < best: best = f
            cnt += 1
            if cnt > 500000: break
        return best
    cls = []; hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines_raw):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'score':s, 'c3':cc, 'H':h, 'sc':sc, 'comp_canon':cc2, 'canon':canon})
        hm[(s,cc,h)].append(i); cm[canon] = i
    for cl in cls: cl['comp_cid'] = cm.get(cl['comp_canon'], -1)
    def lk(adj):
        s = ss(adj); cc = c3(adj); h = H(adj)
        cids = hm.get((s,cc,h))
        if not cids: return -1
        if len(cids) == 1: return cids[0]
        return cm.get(cp(adj), -1)
    mn = {}; mid = 0
    for cl in cls:
        ci = cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci] = mid; mid += 1
        else:
            mn[ci] = mid; comp = cl['comp_cid']
            if comp >= 0 and comp != ci: mn[comp] = mid
            mid += 1
    mi = {}
    for mm in range(mid):
        cids = [cl['cid'] for cl in cls if mn[cl['cid']] == mm]
        cl0 = cls[cids[0]]; mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}
    edges = set()
    for cl in cls:
        b = int(lines_raw[cl['cid']], 2)
        for arc in range(m_total):
            fb = b ^ (1 << (m_total-1-arc))
            fa = b2a(fb)
            nb = lk(fa)
            if nb >= 0 and nb != cl['cid']:
                a2, b2 = mn[cl['cid']], mn[nb]
                if a2 != b2: edges.add((min(a2,b2), max(a2,b2)))
    V = mid
    A = np.zeros((V, V))
    for (a, b) in edges: A[a][b] = 1; A[b][a] = 1
    return V, mi, edges, A

print("=" * 72)
print("  HOFFMAN BOUND AND LIE RANK")
print("  opus-2026-03-24-S314")
print("=" * 72)

for n in [4, 5, 6, 7]:
    print(f"\n{'='*72}")
    print(f"  n = {n}")
    print(f"{'='*72}")

    V, mi, edges, A = build_metagraph(n)
    E = len(edges)
    evals = sorted(np.linalg.eigvalsh(A))

    lam_max = evals[-1]
    lam_min = evals[0]
    hoffman = 1 + lam_max / abs(lam_min)
    hoffman_ceil = int(np.ceil(hoffman))

    print(f"  V={V}, E={E}")
    print(f"  λ_max = {lam_max:.4f}")
    print(f"  λ_min = {lam_min:.4f}")
    print(f"  Hoffman bound: χ ≥ 1 + {lam_max:.4f}/{abs(lam_min):.4f} = {hoffman:.4f}")
    print(f"  ⌈Hoffman⌉ = {hoffman_ceil}")
    print(f"  Actual χ = n-1 = {n-1}")

    # Is it Hoffman-colorable? (χ = Hoffman bound?)
    is_hoffman = (hoffman_ceil == n-1)
    print(f"  Hoffman colorable? {is_hoffman}")
    if is_hoffman:
        print(f"  *** χ = Hoffman bound! The spectrum DETERMINES χ. ***")

    # Check: λ_max/|λ_min| = n-2?
    ratio = lam_max / abs(lam_min)
    print(f"  λ_max/|λ_min| = {ratio:.4f}")
    print(f"  n-2 = {n-2}")
    print(f"  Ratio / (n-2) = {ratio/(n-2):.4f}")

    # The Lie algebra A_{n-1} has:
    # rank = n-1, Coxeter number h = n, dim = n²-1
    # Weyl group = S_n
    print(f"\n  LIE ALGEBRA A_{n-1}:")
    print(f"    rank = {n-1} = χ")
    print(f"    Coxeter number h = {n}")
    print(f"    |W| = {factorial(n)} = n!")
    print(f"    dim = {n**2 - 1}")

    # The n-2 strips of δ_{n-2}
    m = comb(n-1, 2)
    print(f"\n  STAIRCASE δ_{n-2}:")
    print(f"    {n-2} strips")
    print(f"    {m} tiles")
    print(f"    χ = n-1 = #strips + 1")

    # Spectrum details
    print(f"\n  FULL SPECTRUM:")
    for i in range(min(5, V)):
        print(f"    λ_{i} = {evals[-(i+1)]:.4f}")
    print(f"    ...")
    for i in range(min(3, V)):
        print(f"    λ_{V-1-i} = {evals[i]:.4f}")

print(f"\n{'='*72}")
print("  SYNTHESIS: WHY χ = n-1")
print(f"{'='*72}")
print("""
  FIVE CONVERGING EXPLANATIONS:

  1. HOFFMAN BOUND: χ ≥ 1 + λ_max/|λ_min|.
     If this is tight, χ is determined by the spectrum.
     Testing: is the Hoffman bound exactly n-1?

  2. LIE ALGEBRA RANK: A_{n-1} has rank n-1.
     The metagraph = Q_m / W(A_{n-1}). The chromatic number
     equals the rank of the Lie algebra whose Weyl group acts.
     Rank = number of simple roots = number of "independent directions."

  3. STRIP COUNT + 1: The staircase δ_{n-2} has n-2 strips.
     χ = n-1 = #strips + 1. The +1 comes from the base path
     (which is NOT a strip but adds one dimension).

  4. KNESER ANALOGY: K(n,2) has χ = n-2 (Lovász 1978).
     Our graph has χ = n-1 = K(n,2) + 1. The "+1" might come
     from the orientation structure that tournaments add to graphs.

  5. FLIP GRAPH UNIVERSALITY: Flip graphs of combinatorial objects
     (matchings, triangulations, tournaments) tend to have χ linear in n.
     The coefficient ≈ 1 and offset ≈ -1 or +1 seems universal.
""")

print("DONE.")
