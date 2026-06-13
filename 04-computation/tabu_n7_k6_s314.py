#!/usr/bin/env python3
"""Quick tabu search: just k=6 for n=7, with 10 restarts and 2M iterations each."""
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

n = 7
print(f"Building metagraph n={n}...")
V, mi, edges, adj = build_metagraph(n)
print(f"V={V}, E={len(edges)}")

k = 6
print(f"\nTabu search for {k}-coloring with 10 restarts, 2M iters each:")

for seed in range(10):
    rng = np.random.RandomState(seed * 31 + 7)
    colors = rng.randint(0, k, size=V)
    
    conflicts = 0
    for v in range(V):
        for u in adj[v]:
            if u > v and colors[u] == colors[v]:
                conflicts += 1
    
    tabu = np.zeros((V, k), dtype=int)
    best_conflicts = conflicts
    best_colors = colors.copy()
    
    for it in range(2000000):
        conflict_verts = []
        for v in range(V):
            for u in adj[v]:
                if colors[u] == colors[v]:
                    conflict_verts.append(v)
                    break
        
        if not conflict_verts:
            print(f"  Seed {seed}: ZERO conflicts at iter {it}!")
            # Verify
            valid = True
            for (a, b) in edges:
                if colors[a] == colors[b]:
                    valid = False; break
            if valid:
                print(f"  *** VALID {k}-COLORING FOUND! ***")
                color_classes = defaultdict(list)
                for v in range(V):
                    color_classes[colors[v]].append(v)
                for c in range(k):
                    members = color_classes[c]
                    h_vals = [mi[v]['H'] for v in members]
                    sc_count = sum(1 for v in members if mi[v]['sc'])
                    print(f"    Color {c}: {len(members)} vertices, SC={sc_count}, H in [{min(h_vals)}..{max(h_vals)}]")
                sys.exit(0)
            break
        
        best_delta = float('inf')
        best_moves = []
        
        for v in conflict_verts:
            old_c = colors[v]
            cur = sum(1 for u in adj[v] if colors[u] == old_c)
            for new_c in range(k):
                if new_c == old_c: continue
                nc = sum(1 for u in adj[v] if colors[u] == new_c)
                delta = nc - cur
                if tabu[v][new_c] <= it or (conflicts + delta < best_conflicts):
                    if delta < best_delta:
                        best_delta = delta
                        best_moves = [(v, new_c)]
                    elif delta == best_delta:
                        best_moves.append((v, new_c))
        
        if not best_moves:
            v = conflict_verts[rng.randint(len(conflict_verts))]
            new_c = rng.randint(k)
            while new_c == colors[v]: new_c = rng.randint(k)
            best_moves = [(v, new_c)]
            best_delta = sum(1 for u in adj[v] if colors[u] == new_c) - sum(1 for u in adj[v] if colors[u] == colors[v])
        
        v, new_c = best_moves[rng.randint(len(best_moves))]
        old_c = colors[v]
        colors[v] = new_c
        conflicts += best_delta
        tabu[v][old_c] = it + 7 + rng.randint(10)
        
        if conflicts < best_conflicts:
            best_conflicts = conflicts
            best_colors = colors.copy()
    
    print(f"  Seed {seed}: best = {best_conflicts} conflicts")

# If none found with k=6, try k=7
print(f"\nTrying k=7:")
k = 7
for seed in range(3):
    rng = np.random.RandomState(seed * 53 + 11)
    colors = rng.randint(0, k, size=V)
    conflicts = sum(1 for v in range(V) for u in adj[v] if u > v and colors[u] == colors[v])
    tabu = np.zeros((V, k), dtype=int)
    best_conflicts = conflicts
    
    for it in range(1000000):
        conflict_verts = [v for v in range(V) if any(colors[u] == colors[v] for u in adj[v])]
        if not conflict_verts:
            valid = all(colors[a] != colors[b] for (a,b) in edges)
            if valid:
                print(f"  *** VALID {k}-COLORING at seed {seed}, iter {it} ***")
                break
            break
        
        best_delta = float('inf'); best_moves = []
        for v in conflict_verts:
            old_c = colors[v]
            cur = sum(1 for u in adj[v] if colors[u] == old_c)
            for new_c in range(k):
                if new_c == old_c: continue
                nc = sum(1 for u in adj[v] if colors[u] == new_c)
                delta = nc - cur
                if tabu[v][new_c] <= it or (conflicts + delta < best_conflicts):
                    if delta < best_delta: best_delta = delta; best_moves = [(v, new_c)]
                    elif delta == best_delta: best_moves.append((v, new_c))
        
        if not best_moves:
            v = conflict_verts[rng.randint(len(conflict_verts))]
            new_c = rng.randint(k)
            while new_c == colors[v]: new_c = rng.randint(k)
            best_moves = [(v, new_c)]
            best_delta = sum(1 for u in adj[v] if colors[u] == new_c) - sum(1 for u in adj[v] if colors[u] == colors[v])
        
        v, new_c = best_moves[rng.randint(len(best_moves))]
        colors[v] = new_c; conflicts += best_delta
        tabu[v][colors[v]] = it + 7 + rng.randint(10)
        if conflicts < best_conflicts: best_conflicts = conflicts
    
    print(f"  Seed {seed}: best = {best_conflicts}")

print("DONE.")
