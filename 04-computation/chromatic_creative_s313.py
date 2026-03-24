#!/usr/bin/env python3
"""
chromatic_creative_s313.py — Creative chromatic number attacks
opus-2026-03-24-S313

Find the EXACT chromatic number and an OPTIMAL coloring.
Try creative non-greedy strategies.
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

def exact_chromatic(V, edges, adj_list, max_colors):
    """Find exact chromatic number by backtracking."""
    for k in range(1, max_colors + 1):
        colors = [-1] * V
        if backtrack_color(0, V, adj_list, colors, k):
            return k, colors[:]
    return max_colors + 1, None

def backtrack_color(node, V, adj, colors, k):
    if node == V: return True
    used = {colors[nb] for nb in adj[node] if colors[nb] >= 0}
    for c in range(k):
        if c not in used:
            colors[node] = c
            if backtrack_color(node + 1, V, adj, colors, k):
                return True
            colors[node] = -1
    return False

print("=" * 72)
print("  CREATIVE CHROMATIC NUMBER ANALYSIS")
print("  opus-2026-03-24-S313")
print("=" * 72)

for n in [5, 6]:
    print(f"\n{'#'*72}")
    print(f"  n = {n}")
    print(f"{'#'*72}")

    V, mi, edges, adj = build_metagraph(n)
    E = len(edges)
    print(f"  V={V}, E={E}")

    # Exact chromatic number (brute force for small V)
    if V <= 40:
        print(f"\n  EXACT CHROMATIC NUMBER (backtracking):")
        chi, opt_colors = exact_chromatic(V, edges, adj, n+2)
        print(f"    χ = {chi}")
        print(f"    Conjectured n-1 = {n-1}")
        print(f"    Match? {'YES ✓' if chi == n-1 else 'NO ✗'}")

        if opt_colors:
            # Analyze the optimal coloring
            print(f"\n    OPTIMAL COLORING STRUCTURE:")
            color_classes = defaultdict(list)
            for v in range(V):
                color_classes[opt_colors[v]].append(v)

            for c in range(chi):
                members = color_classes[c]
                h_vals = [mi[v]['H'] for v in members]
                sc_count = sum(1 for v in members if mi[v]['sc'])
                print(f"      Color {c}: {len(members)} classes, "
                      f"H={sorted(h_vals)}, SC={sc_count}")

            # What FUNCTION of the class determines the color?
            # Check: is the color determined by H mod chi?
            h_mod = [mi[v]['H'] % chi for v in range(V)]
            color_matches_hmod = sum(1 for v in range(V) if opt_colors[v] == h_mod[v])
            print(f"\n    Color = H mod {chi}? {color_matches_hmod}/{V} match")

            # Check: is color determined by c3 mod chi?
            c3_mod = [mi[v]['c3'] % chi for v in range(V)]
            color_matches_c3mod = sum(1 for v in range(V) if opt_colors[v] == c3_mod[v])
            print(f"    Color = c3 mod {chi}? {color_matches_c3mod}/{V} match")

            # Check: is color determined by score type?
            score_types = defaultdict(list)
            for v in range(V):
                score_types[mi[v]['score']].append(v)
            n_score_types = len(score_types)
            print(f"    Distinct score types: {n_score_types}")

            # Check: is color determined by distance from transitive mod chi?
            from collections import deque
            dist = [-1] * V; dist[0] = 0; q = deque([0])
            while q:
                u = q.popleft()
                for v in adj[u]:
                    if dist[v] < 0: dist[v] = dist[u]+1; q.append(v)
            dist_mod = [d % chi for d in dist]
            color_matches_dmod = sum(1 for v in range(V) if opt_colors[v] == dist_mod[v])
            print(f"    Color = dist_from_trans mod {chi}? {color_matches_dmod}/{V} match")

            # Try: (H + c3) mod chi
            hc3_mod = [(mi[v]['H'] + mi[v]['c3']) % chi for v in range(V)]
            color_matches_hc3 = sum(1 for v in range(V) if opt_colors[v] == hc3_mod[v])
            print(f"    Color = (H+c3) mod {chi}? {color_matches_hc3}/{V} match")

            # Try: (H - c3) mod chi
            hc3d_mod = [(mi[v]['H'] - mi[v]['c3']) % chi for v in range(V)]
            color_matches_hc3d = sum(1 for v in range(V) if opt_colors[v] == hc3d_mod[v])
            print(f"    Color = (H-c3) mod {chi}? {color_matches_hc3d}/{V} match")

            # Actually check if H mod chi IS a valid coloring
            h_mod_valid = True
            for (a, b) in edges:
                if mi[a]['H'] % chi == mi[b]['H'] % chi:
                    h_mod_valid = False
                    break
            print(f"\n    H mod {chi} is valid coloring? {h_mod_valid}")

            # Check distance mod chi
            dist_mod_valid = True
            for (a, b) in edges:
                if dist[a] % chi == dist[b] % chi:
                    dist_mod_valid = False
                    break
            print(f"    dist mod {chi} is valid coloring? {dist_mod_valid}")

            # What about the METAGRAPH DISTANCE from transitive as coloring?
            # This gives at most diameter+1 colors.
            # But distance is well-defined: adjacent nodes have dist differ by ≤ 1
            # So dist mod k is valid iff all edges have |Δdist| not ≡ 0 mod k.
            # Adjacent classes: |Δdist| = 0 or 1 (in the BFS tree, adjacent nodes
            # can be at the same level or adjacent levels).
            # |Δdist| = 0: same-level edges. These break dist-based coloring.

            same_level = sum(1 for (a,b) in edges if dist[a] == dist[b])
            print(f"\n    Same-distance edges: {same_level}/{E} ({100*same_level/E:.1f}%)")

            # The clique number
            max_clique = 2
            for (a, b) in edges:
                for c in adj[a]:
                    if c > b and (min(b,c), max(b,c)) in edges:
                        max_clique = 3
                        break
            print(f"\n    Clique number ω ≥ {max_clique}")

print(f"\n{'='*72}")
print("  THE CHROMATIC STRUCTURE")
print(f"{'='*72}")
print("""
  THE OPTIMAL COLORING reveals:

  At n=5 (χ=4): The 10 nodes split into 4 color classes.
  At n=6: χ = ? (compute above).

  The coloring is NOT simply H mod χ or dist mod χ.
  Same-distance edges exist (nodes at the same BFS level are connected).
  These prevent distance-based coloring.

  THE DEEPER STRUCTURE:
  χ = n-1 would mean: tournament iso classes can be partitioned into
  n-1 independent sets (no edges within). This is a statement about
  the structure of SINGLE-ARC-FLIP connections.

  If two classes have the same color, they cannot be connected by ANY
  single arc flip. This means: no arc reversal can transform one into
  the other. They are "structurally incompatible" under single-arc changes.

  The n-1 colors correspond to n-1 "structural families" of tournaments
  that are mutually non-adjacent in the metagraph.
""")

print("DONE.")
