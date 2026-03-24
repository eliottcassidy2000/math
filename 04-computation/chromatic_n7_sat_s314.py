#!/usr/bin/env python3
"""
chromatic_n7_sat_s314.py — Prove χ(G_7/Z_2) = 6 using clique lower bound + greedy upper bound
opus-2026-03-24-S314

Strategy:
1. Find clique number ω — gives ω ≤ χ
2. Try (n-1)-coloring via DSATUR (best-first) — if found, χ ≤ n-1
3. If ω = n-1 and (n-1)-coloring exists, then χ = n-1 exactly
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

def find_max_clique(V, adj, edges):
    """Greedy + BFS clique search"""
    # Start from highest-degree vertex
    degrees = {v: len(adj[v]) for v in range(V)}
    
    # Try expanding cliques from every vertex
    best_clique = []
    for start in range(V):
        clique = [start]
        candidates = sorted(adj[start], key=lambda v: -degrees[v])
        for c in candidates:
            if all((min(c, x), max(c, x)) in edges for x in clique):
                clique.append(c)
        if len(clique) > len(best_clique):
            best_clique = clique
    
    # Also try random greedy starts
    rng = np.random.RandomState(42)
    for _ in range(1000):
        order = rng.permutation(V)
        clique = [order[0]]
        for v in order[1:]:
            if all((min(v, x), max(v, x)) in edges for x in clique):
                clique.append(v)
        if len(clique) > len(best_clique):
            best_clique = clique
    
    return best_clique

def dsatur_color(V, adj):
    """DSATUR greedy coloring — uses saturation degree heuristic"""
    colors = [-1] * V
    sat = [0] * V  # saturation: number of distinct colors in neighborhood
    adj_colors = [set() for _ in range(V)]
    degrees = [len(adj[v]) for v in range(V)]
    
    # Start with highest-degree vertex
    v0 = max(range(V), key=lambda v: degrees[v])
    colors[v0] = 0
    for u in adj[v0]:
        adj_colors[u].add(0)
        sat[u] = len(adj_colors[u])
    
    for step in range(1, V):
        # Pick uncolored vertex with highest saturation (break ties by degree)
        best = -1; best_sat = -1; best_deg = -1
        for v in range(V):
            if colors[v] >= 0: continue
            if sat[v] > best_sat or (sat[v] == best_sat and degrees[v] > best_deg):
                best = v; best_sat = sat[v]; best_deg = degrees[v]
        
        v = best
        # Assign smallest available color
        used = adj_colors[v]
        c = 0
        while c in used: c += 1
        colors[v] = c
        
        # Update neighbors
        for u in adj[v]:
            if colors[u] < 0:
                adj_colors[u].add(c)
                sat[u] = len(adj_colors[u])
    
    return max(colors) + 1, colors

def backtrack_color_k(V, adj, k, timeout=5000000):
    """Try to find k-coloring by backtracking with pruning"""
    colors = [-1] * V
    # Order vertices by degree (descending)
    order = sorted(range(V), key=lambda v: -len(adj[v]))
    count = [0]
    
    def solve(idx):
        if count[0] > timeout: return False
        if idx == V: return True
        count[0] += 1
        v = order[idx]
        used = {colors[u] for u in adj[v] if colors[u] >= 0}
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
print("  χ(G_n/Z_2) VIA CLIQUE + DSATUR")
print("  opus-2026-03-24-S314")
print("=" * 72)

for n in [4, 5, 6, 7]:
    print(f"\n{'='*72}")
    print(f"  n = {n}")
    print(f"{'='*72}")
    
    V, mi, edges, adj = build_metagraph(n)
    E = len(edges)
    print(f"  V={V}, E={E}")
    
    # 1. Find max clique
    clique = find_max_clique(V, adj, edges)
    omega = len(clique)
    print(f"\n  CLIQUE NUMBER: ω ≥ {omega}")
    print(f"    Clique: {clique[:10]}{'...' if len(clique) > 10 else ''}")
    H_vals = [mi[v]['H'] for v in clique]
    print(f"    H values: {sorted(H_vals)}")
    
    # 2. DSATUR coloring
    chi_upper, dsatur_colors = dsatur_color(V, adj)
    print(f"\n  DSATUR COLORING: χ ≤ {chi_upper}")
    
    # 3. Try to improve with exact backtracking
    if chi_upper > n - 1:
        print(f"  Trying {n-1}-coloring by backtracking...")
        bt = backtrack_color_k(V, adj, n - 1)
        if bt:
            chi_upper = n - 1
            print(f"  Found {n-1}-coloring! χ ≤ {n-1}")
        else:
            print(f"  No {n-1}-coloring found within timeout")
    
    # 4. Summary
    print(f"\n  BOUNDS: {omega} ≤ χ ≤ {chi_upper}")
    print(f"  n-1 = {n-1}")
    if omega == n - 1 and chi_upper == n - 1:
        print(f"  *** χ = {n-1} PROVED! ω = χ = n-1 ***")
    elif chi_upper == n - 1:
        print(f"  χ ≤ n-1 confirmed. Lower bound ω = {omega} < n-1.")
    
    # 5. Analyze the coloring structure
    if chi_upper <= n:
        color_classes = defaultdict(list)
        cols = dsatur_colors if chi_upper == max(dsatur_colors)+1 else (bt if 'bt' in dir() and bt else dsatur_colors)
        for v in range(V):
            color_classes[cols[v]].append(v)
        
        max_colors_used = max(cols) + 1
        print(f"\n  COLORING STRUCTURE ({max_colors_used} colors):")
        for c in range(max_colors_used):
            members = color_classes[c]
            h_vals = sorted([mi[v]['H'] for v in members])
            sc_count = sum(1 for v in members if mi[v]['sc'])
            print(f"    Color {c}: {len(members)} classes, SC={sc_count}, "
                  f"H range [{min(h_vals)}..{max(h_vals)}]")

print("\nDONE.")
