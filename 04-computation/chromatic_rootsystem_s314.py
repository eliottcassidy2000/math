#!/usr/bin/env python3
"""
chromatic_rootsystem_s314.py — Root system proof of χ = n-1
opus-2026-03-24-S314

THEORY:
Tournament arcs = positive roots of A_{n-1}.
An arc flip = reflecting through one root.
The metagraph G_n/Z_2 = the quotient of the root polytope flip graph.

In the root system A_{n-1}:
- n-1 simple roots α_1,...,α_{n-1}
- Each positive root α_{ij} = α_i + α_{i+1} + ... + α_{j-1}
- Flipping arc (i,j) = reflecting through α_{ij}

CLAIM: The n-1 simple roots give n-1 "independent directions" for coloring.
The color of a tournament class = which simple root was last flipped.

More precisely: define f(T) = Σ_arcs weight(arc) mod (n-1).
Check if f is a valid coloring.
"""

import sys, subprocess
from math import comb
from itertools import permutations
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_metagraph_with_weights(n):
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
    edge_arcs = {}  # (a,b) -> which arcs connect them
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
                if a2 != b2:
                    e = (min(a2,b2), max(a2,b2))
                    edges.add(e)
                    if e not in edge_arcs: edge_arcs[e] = set()
                    edge_arcs[e].add(arc)
    V = mid
    adj_list = defaultdict(set)
    for (a, b) in edges: adj_list[a].add(b); adj_list[b].add(a)
    
    # Arc weights: arc k corresponds to pair P[k] = (i,j)
    # Root weight = j - i (length of the root α_{ij})
    arc_weights = {}
    for k, (i,j) in enumerate(P):
        arc_weights[k] = j - i  # simple root has weight 1
    
    return V, mi, edges, adj_list, edge_arcs, arc_weights, P, cls, mn, lines_raw

print("=" * 72)
print("  ROOT SYSTEM COLORING ANALYSIS")
print("  opus-2026-03-24-S314")
print("=" * 72)

for n in [4, 5, 6]:
    print(f"\n{'='*72}")
    print(f"  n = {n}")
    print(f"{'='*72}")
    
    V, mi, edges, adj, edge_arcs, arc_weights, P, cls, mn, lines_raw = build_metagraph_with_weights(n)
    E = len(edges)
    m_total = comb(n, 2)
    print(f"  V={V}, E={E}")
    
    # For each edge, what is the root weight of the connecting arc?
    print(f"\n  EDGE ROOT WEIGHTS:")
    weight_dist = defaultdict(int)
    for e in edges:
        arcs = edge_arcs[e]
        for a in arcs:
            weight_dist[arc_weights[a]] += 1
    for w in sorted(weight_dist.keys()):
        i_j_pairs = [(i,j) for k,(i,j) in enumerate(P) if j-i == w]
        print(f"    weight {w} (roots of length {w}): {weight_dist[w]} arc-edge incidences, pairs: {i_j_pairs[:5]}")
    
    # TRY: coloring by score sum mod (n-1)
    print(f"\n  COLORING ATTEMPTS:")
    
    # f1: Σ (j-i) * direction(i→j) mod (n-1) 
    # where direction = +1 if i→j, -1 if j→i
    for cl in cls:
        cid = cl['cid']
        if cid >= len(lines_raw): continue
        bits = int(lines_raw[cid], 2)
        # Compute weight sum
        wsum = 0
        for k, (i,j) in enumerate(P):
            if bits & (1 << (m_total-1-k)):
                wsum += (j - i)
            else:
                wsum -= (j - i)
        cl['wsum'] = wsum
    
    # Map to merged classes
    for mm in range(V):
        cids = [cl['cid'] for cl in cls if mn.get(cl['cid']) == mm]
        if cids:
            wsums = set()
            for cid in cids:
                cl = cls[cid]
                if 'wsum' in cl:
                    wsums.add(cl['wsum'])
            mi[mm]['wsums'] = wsums
    
    # Check if wsum mod (n-1) is a valid coloring
    for modulus in [n-1, n-2, n]:
        valid = True
        for (a, b) in edges:
            wa = mi[a].get('wsums', set())
            wb = mi[b].get('wsums', set())
            if not wa or not wb: continue
            # All elements in wa should have same value mod modulus
            wa_mods = {w % modulus for w in wa}
            wb_mods = {w % modulus for w in wb}
            if len(wa_mods) > 1 or len(wb_mods) > 1:
                valid = False  # not well-defined on iso classes
                break
            if wa_mods & wb_mods:
                valid = False
                break
        
        # Also check well-definedness
        well_defined = True
        for mm in range(V):
            ws = mi[mm].get('wsums', set())
            if len({w % modulus for w in ws}) > 1:
                well_defined = False
                break
        
        colors_used = len({list(mi[mm].get('wsums', {0}))[0] % modulus for mm in range(V) if mi[mm].get('wsums')})
        print(f"    wsum mod {modulus}: well-defined={well_defined}, valid={valid}, colors={colors_used}")
    
    # TRY: H mod (n-1)
    h_mod_valid = True
    for (a, b) in edges:
        if mi[a]['H'] % (n-1) == mi[b]['H'] % (n-1):
            h_mod_valid = False
            break
    print(f"    H mod {n-1}: valid={h_mod_valid}")
    
    # TRY: c3 mod (n-1)
    c3_mod_valid = True
    for (a, b) in edges:
        if mi[a]['c3'] % (n-1) == mi[b]['c3'] % (n-1):
            c3_mod_valid = False
            break
    print(f"    c3 mod {n-1}: valid={c3_mod_valid}")
    
    # TRY: (H + c3) mod (n-1)
    hc3_mod_valid = True
    for (a, b) in edges:
        if (mi[a]['H'] + mi[a]['c3']) % (n-1) == (mi[b]['H'] + mi[b]['c3']) % (n-1):
            hc3_mod_valid = False
            break
    print(f"    (H+c3) mod {n-1}: valid={hc3_mod_valid}")
    
    # TRY: H mod 2
    h2_valid = True
    for (a, b) in edges:
        if mi[a]['H'] % 2 == mi[b]['H'] % 2:
            h2_valid = False
            break
    print(f"    H mod 2: valid={h2_valid}")
    
    # What about score sequence partitions?
    print(f"\n  SCORE SEQUENCE PARTITION:")
    score_types = defaultdict(list)
    for mm in range(V):
        score_types[mi[mm]['score']].append(mm)
    n_types = len(score_types)
    print(f"    {n_types} distinct score types")
    
    # Are all edges BETWEEN score types (not within)?
    within = 0; between = 0
    for (a, b) in edges:
        if mi[a]['score'] == mi[b]['score']:
            within += 1
        else:
            between += 1
    print(f"    Within same score type: {within}/{E}")
    print(f"    Between score types: {between}/{E}")
    
    # H parity on edges
    same_h_parity = sum(1 for (a,b) in edges if mi[a]['H'] % 2 == mi[b]['H'] % 2)
    print(f"\n  H parity: same={same_h_parity}/{E}, diff={E-same_h_parity}/{E}")

print("\nDONE.")
