#!/usr/bin/env python3
"""
chromatic_n7_deep_s314.py — Deep clique search + smarter coloring for n=7
opus-2026-03-24-S314

At n=4,5,6: ω = χ = n-1. Is the graph PERFECT?
Need: ω(G_7/Z_2) = 6 AND a 6-coloring.
"""

import sys, subprocess
from math import comb
from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

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
        cl0 = cls[cids[0]]; mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'c3': cl0['c3'], 'score': cl0['score']}
    edges = set()
    def lk(adj):
        s = ss(adj); cc = c3(adj); h = H(adj)
        cids = hm.get((s,cc,h))
        if not cids: return -1
        if len(cids) == 1: return cids[0]
        return cm.get(cp(adj), -1)
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
    adj_list = defaultdict(set)
    for (a, b) in edges: adj_list[a].add(b); adj_list[b].add(a)
    return V, mi, edges, adj_list

def bron_kerbosch_max(V, adj, edges, target_size=None):
    """Bron-Kerbosch with pivoting for maximum clique"""
    best = []
    count = [0]
    
    def bk(R, P, X):
        nonlocal best
        count[0] += 1
        if count[0] > 10000000: return  # safety limit
        if not P and not X:
            if len(R) > len(best):
                best = list(R)
                print(f"    New clique size {len(best)}: {best[:10]}")
                if target_size and len(best) >= target_size:
                    return
            return
        
        # Pruning: can we possibly beat best?
        if len(R) + len(P) <= len(best): return
        
        # Pivot selection: vertex in P∪X with most neighbors in P
        pivot_candidates = P | X
        if not pivot_candidates: return
        pivot = max(pivot_candidates, key=lambda u: len(adj[u] & P))
        
        for v in list(P - adj[pivot]):
            new_P = P & adj[v]
            new_X = X & adj[v]
            bk(R | {v}, new_P, new_X)
            if target_size and len(best) >= target_size: return
            P.remove(v)
            X.add(v)
    
    # Order vertices by degree (descending) for better pruning
    vertices = sorted(range(V), key=lambda v: -len(adj[v]))
    P = set(vertices)
    bk(set(), P, set())
    print(f"    BK iterations: {count[0]}")
    return best

def backtrack_color_ordered(V, adj, k, order=None, timeout=20000000):
    """Backtracking k-coloring with custom vertex ordering and forward checking"""
    colors = [-1] * V
    if order is None:
        order = list(range(V))
    count = [0]
    
    # Precompute neighbor lists in ordering
    adj_in_order = [[] for _ in range(V)]
    pos = [0] * V
    for i, v in enumerate(order):
        pos[v] = i
    for i, v in enumerate(order):
        for u in adj[v]:
            if pos[u] < i:  # already colored
                adj_in_order[i].append(pos[u])
    
    def solve(idx):
        if count[0] > timeout: return False
        if idx == V: return True
        count[0] += 1
        v = order[idx]
        used = set()
        for u in adj[v]:
            if colors[u] >= 0:
                used.add(colors[u])
        for c in range(k):
            if c not in used:
                colors[v] = c
                if solve(idx + 1): return True
                colors[v] = -1
        return False
    
    if solve(0):
        return colors
    return None

print("=" * 72)
print("  DEEP ANALYSIS: χ(G_7/Z_2)")
print("  opus-2026-03-24-S314")
print("=" * 72)

n = 7
print(f"\nBuilding metagraph for n={n}...")
V, mi, edges, adj = build_metagraph(n)
E = len(edges)
print(f"V={V}, E={E}")

# Degree statistics
degrees = [len(adj[v]) for v in range(V)]
print(f"\nDegree stats: min={min(degrees)}, max={max(degrees)}, mean={np.mean(degrees):.1f}")
print(f"Vertices with max degree: {sum(1 for d in degrees if d == max(degrees))}")

# 1. EXACT MAX CLIQUE via Bron-Kerbosch
print(f"\n{'='*72}")
print(f"  1. EXACT MAX CLIQUE (Bron-Kerbosch)")
print(f"{'='*72}")
clique = bron_kerbosch_max(V, adj, edges, target_size=6)
omega = len(clique)
print(f"\n  ω(G_7/Z_2) = {omega}")
print(f"  Clique: {clique}")
H_vals = [mi[v]['H'] for v in clique]
print(f"  H values: {sorted(H_vals)}")
sc_vals = [mi[v]['sc'] for v in clique]
print(f"  SC types: {sc_vals}")

# 2. Try multiple vertex orderings for coloring
print(f"\n{'='*72}")
print(f"  2. COLORING ATTEMPTS")
print(f"{'='*72}")

target_colors = n - 1  # = 6

# Ordering 1: by degree (descending)
order_deg = sorted(range(V), key=lambda v: -len(adj[v]))
print(f"\n  Attempt 1: degree ordering, k={target_colors}")
result = backtrack_color_ordered(V, adj, target_colors, order_deg, timeout=20000000)
if result:
    print(f"  SUCCESS: {target_colors}-coloring found!")
else:
    print(f"  Failed within timeout")

if not result:
    # Ordering 2: by H value
    order_h = sorted(range(V), key=lambda v: mi[v]['H'])
    print(f"\n  Attempt 2: H-value ordering, k={target_colors}")
    result = backtrack_color_ordered(V, adj, target_colors, order_h, timeout=20000000)
    if result:
        print(f"  SUCCESS: {target_colors}-coloring found!")
    else:
        print(f"  Failed within timeout")

if not result:
    # Ordering 3: by c3 value
    order_c3 = sorted(range(V), key=lambda v: -mi[v]['c3'])
    print(f"\n  Attempt 3: c3-value ordering, k={target_colors}")
    result = backtrack_color_ordered(V, adj, target_colors, order_c3, timeout=20000000)
    if result:
        print(f"  SUCCESS: {target_colors}-coloring found!")
    else:
        print(f"  Failed within timeout")

if not result:
    # Try k=7 to at least bound
    print(f"\n  Attempt 4: trying k={target_colors+1}")
    result7 = backtrack_color_ordered(V, adj, target_colors+1, order_deg, timeout=5000000)
    if result7:
        print(f"  {target_colors+1}-coloring found!")
    else:
        print(f"  Failed")

# 3. Summary
print(f"\n{'='*72}")
print(f"  SUMMARY FOR n=7")
print(f"{'='*72}")
print(f"  ω = {omega}")
print(f"  χ ≤ {target_colors if result else '?'}")
if omega == target_colors and result:
    print(f"  *** ω = χ = {target_colors} = n-1. PROVED! ***")
elif omega == target_colors:
    print(f"  ω = {omega} = n-1 but {target_colors}-coloring not found yet")
else:
    print(f"  Need larger clique (currently {omega}, target {target_colors})")

# 4. Check: is ω < n-1 possible?
# Fractional chromatic number from Hoffman
A = np.zeros((V, V))
for (a, b) in edges: A[a][b] = 1; A[b][a] = 1
evals = sorted(np.linalg.eigvalsh(A))
lam_max = evals[-1]; lam_min = evals[0]
hoffman = 1 + lam_max / abs(lam_min)
print(f"\n  Hoffman bound: χ ≥ {hoffman:.4f}")
print(f"  Lovász theta: need SDP, skipping")

print("\nDONE.")
