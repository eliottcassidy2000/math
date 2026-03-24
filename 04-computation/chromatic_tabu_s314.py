#!/usr/bin/env python3
"""
chromatic_tabu_s314.py — Tabu search for k-coloring of G_7/Z_2
opus-2026-03-24-S314

Tabu search is far superior to backtracking for graph coloring.
Also tries: independent set decomposition, fractional chromatic number.
"""

import sys, subprocess
from math import comb
from itertools import permutations
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

def tabu_color(V, adj, k, max_iter=2000000, seed=42):
    """Tabu search k-coloring (TabuCol algorithm)."""
    rng = np.random.RandomState(seed)
    
    # Random initial coloring
    colors = rng.randint(0, k, size=V)
    
    # Compute conflicts
    conflicts = 0
    for v in range(V):
        for u in adj[v]:
            if u > v and colors[u] == colors[v]:
                conflicts += 1
    
    if conflicts == 0:
        return list(colors)
    
    # Tabu tenure
    tabu = np.zeros((V, k), dtype=int)
    
    best_conflicts = conflicts
    
    for it in range(max_iter):
        # Find the best non-tabu move
        best_delta = float('inf')
        best_moves = []
        
        # Only consider conflicting vertices
        conflict_vertices = []
        for v in range(V):
            for u in adj[v]:
                if colors[u] == colors[v]:
                    conflict_vertices.append(v)
                    break
        
        if not conflict_vertices:
            return list(colors)
        
        for v in conflict_vertices:
            old_c = colors[v]
            # Count conflicts with current color
            cur_conflicts = sum(1 for u in adj[v] if colors[u] == old_c)
            
            for new_c in range(k):
                if new_c == old_c: continue
                # Count conflicts with new color
                new_conflicts = sum(1 for u in adj[v] if colors[u] == new_c)
                delta = new_conflicts - cur_conflicts
                
                # Accept if not tabu, or if aspiration (improves best)
                if tabu[v][new_c] <= it or (conflicts + delta < best_conflicts):
                    if delta < best_delta:
                        best_delta = delta
                        best_moves = [(v, new_c)]
                    elif delta == best_delta:
                        best_moves.append((v, new_c))
        
        if not best_moves:
            # All moves are tabu — pick random conflicting vertex
            v = conflict_vertices[rng.randint(len(conflict_vertices))]
            new_c = rng.randint(k)
            while new_c == colors[v]:
                new_c = rng.randint(k)
            best_moves = [(v, new_c)]
            best_delta = sum(1 for u in adj[v] if colors[u] == new_c) - sum(1 for u in adj[v] if colors[u] == colors[v])
        
        # Make the move
        v, new_c = best_moves[rng.randint(len(best_moves))]
        old_c = colors[v]
        colors[v] = new_c
        conflicts += best_delta
        
        # Update tabu
        tenure = 7 + rng.randint(10)
        tabu[v][old_c] = it + tenure
        
        if conflicts < best_conflicts:
            best_conflicts = conflicts
            if it % 100000 == 0 or conflicts <= 3:
                print(f"    iter {it}: conflicts = {conflicts}")
        
        if conflicts == 0:
            return list(colors)
    
    print(f"    Best conflicts after {max_iter} iterations: {best_conflicts}")
    return None

print("=" * 72)
print("  TABU SEARCH COLORING + STRUCTURAL ANALYSIS")
print("  opus-2026-03-24-S314")
print("=" * 72)

for n in [5, 6, 7]:
    print(f"\n{'='*72}")
    print(f"  n = {n}")
    print(f"{'='*72}")
    
    V, mi, edges, adj = build_metagraph(n)
    E = len(edges)
    print(f"  V={V}, E={E}")
    
    # Try k = n-2, n-1, n with tabu search
    for k in range(max(3, n-2), n+1):
        print(f"\n  Trying k={k}-coloring (tabu search):")
        result = None
        for seed in range(5):  # 5 random restarts
            result = tabu_color(V, adj, k, max_iter=500000, seed=seed*17+42)
            if result:
                # Verify
                valid = True
                for (a, b) in edges:
                    if result[a] == result[b]:
                        valid = False; break
                print(f"    k={k}: {'VALID' if valid else 'INVALID'} coloring found (seed {seed})")
                if valid:
                    break
                else:
                    result = None
        
        if result:
            # Analyze the coloring
            color_classes = defaultdict(list)
            for v in range(V):
                color_classes[result[v]].append(v)
            n_used = len(color_classes)
            print(f"    Colors used: {n_used}")
            for c in sorted(color_classes.keys()):
                members = color_classes[c]
                h_vals = [mi[v]['H'] for v in members]
                sc_count = sum(1 for v in members if mi[v]['sc'])
                print(f"      Color {c}: {len(members)} vertices, SC={sc_count}, "
                      f"H in [{min(h_vals)}..{max(h_vals)}]")
            
            if k < n - 1:
                print(f"    *** χ ≤ {k} < n-1 = {n-1}! Conjecture FAILS! ***")
            elif k == n - 1:
                print(f"    χ ≤ {n-1} = n-1 confirmed!")
            break
        else:
            print(f"    k={k}: no coloring found")
    
    # Independence number
    print(f"\n  INDEPENDENCE NUMBER:")
    # Greedy independent sets
    best_indep = []
    rng = np.random.RandomState(42)
    for trial in range(100):
        order = rng.permutation(V)
        indep = []
        for v in order:
            if all(v not in adj[u] for u in indep):
                indep.append(v)
        if len(indep) > len(best_indep):
            best_indep = indep
    alpha = len(best_indep)
    print(f"    α ≥ {alpha}")
    print(f"    V/α ≤ {V/alpha:.2f} (fractional χ lower bound)")
    print(f"    V/(n-1) = {V/(n-1):.2f} (expected class size if χ=n-1)")

print("\nDONE.")
