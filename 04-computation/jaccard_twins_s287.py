#!/usr/bin/env python3
"""
jaccard_twins_s287.py — Find Jaccard=1 pairs (twin nodes) in the metagraph
opus-2026-03-24-S287

Jaccard(i,j) = |N(i) ∩ N(j)| / |N(i) ∪ N(j)|
When Jaccard = 1: nodes i,j have IDENTICAL neighbor sets = "twins."

If twins exist, merging them might reveal the Petersen graph underneath.
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
        cls.append({'cid':i, 'score':s, 'c3':cc, 'H':h, 'sc':sc,
                    'comp_canon':cc2, 'canon':canon})
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
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'score': cl0['score'], 'c3': cl0['c3']}

    me = set()
    for cl in cls:
        bits = int(lines_raw[cl['cid']], 2)
        for arc in range(m_total):
            fb = bits ^ (1 << (m_total-1-arc))
            fa = b2a(fb)
            nb = lk(fa)
            if nb >= 0 and nb != cl['cid']:
                a2, b2 = mn[cl['cid']], mn[nb]
                if a2 != b2: me.add((min(a2,b2), max(a2,b2)))

    return mid, mi, me

def analyze_jaccard(n):
    V, mi, me = build_metagraph(n)
    print(f"\n{'#'*72}")
    print(f"  n = {n}: JACCARD SIMILARITY (V={V}, E={len(me)})")
    print(f"{'#'*72}")

    # Build neighbor sets
    neighbors = defaultdict(set)
    for (a, b) in me:
        neighbors[a].add(b)
        neighbors[b].add(a)

    # Compute Jaccard for all pairs
    jaccard_pairs = []
    for i in range(V):
        for j in range(i+1, V):
            ni = neighbors[i]; nj = neighbors[j]
            intersection = len(ni & nj)
            union = len(ni | nj)
            jacc = intersection / union if union > 0 else 0
            jaccard_pairs.append((i, j, jacc, intersection, union))

    # Sort by Jaccard descending
    jaccard_pairs.sort(key=lambda x: -x[2])

    # Show top pairs
    print(f"\n  TOP JACCARD PAIRS:")
    print(f"  {'Pair':>10} {'Jacc':>6} {'|∩|':>4} {'|∪|':>4} {'H_i':>4} {'H_j':>4} "
          f"{'SC':>8} {'edge?':>6} {'deg_i':>6} {'deg_j':>6}")

    for i, j, jacc, inter, union in jaccard_pairs[:20]:
        is_edge = (min(i,j), max(i,j)) in me
        sc_type = ('SC' if mi[i]['sc'] else 'NS') + '-' + ('SC' if mi[j]['sc'] else 'NS')
        di = len(neighbors[i]); dj = len(neighbors[j])
        print(f"  ({i:>3},{j:>3}) {jacc:>6.3f} {inter:>4} {union:>4} "
              f"{mi[i]['H']:>4} {mi[j]['H']:>4} {sc_type:>8} "
              f"{'YES' if is_edge else 'no':>6} {di:>6} {dj:>6}")

    # Find Jaccard = 1 pairs
    twins = [(i, j) for i, j, jacc, _, _ in jaccard_pairs if jacc == 1.0]
    print(f"\n  JACCARD = 1 PAIRS (twins): {len(twins)}")
    for (i, j) in twins:
        ni = neighbors[i]; nj = neighbors[j]
        is_edge = (min(i,j), max(i,j)) in me
        print(f"    ({i},{j}): H=({mi[i]['H']},{mi[j]['H']}), "
              f"SC=({mi[i]['sc']},{mi[j]['sc']}), "
              f"neighbors={sorted(ni)}, edge={is_edge}, "
              f"score=({mi[i]['score']},{mi[j]['score']})")

    # If twins exist, what happens when we merge them?
    if twins:
        print(f"\n  MERGING TWINS → QUOTIENT GRAPH:")
        # Build equivalence classes from twins
        equiv = list(range(V))
        for (i, j) in twins:
            # Union-find
            ri, rj = i, j
            while equiv[ri] != ri: ri = equiv[ri]
            while equiv[rj] != rj: rj = equiv[rj]
            if ri != rj: equiv[max(ri, rj)] = min(ri, rj)

        # Flatten
        for i in range(V):
            r = i
            while equiv[r] != r: r = equiv[r]
            equiv[i] = r

        # Map to new node ids
        roots = sorted(set(equiv))
        new_id = {r: idx for idx, r in enumerate(roots)}

        # New edges
        new_edges = set()
        for (a, b) in me:
            na = new_id[equiv[a]]; nb = new_id[equiv[b]]
            if na != nb: new_edges.add((min(na, nb), max(na, nb)))

        V_new = len(roots)
        E_new = len(new_edges)

        # Degree sequence
        new_degs = [0] * V_new
        for (a, b) in new_edges:
            new_degs[a] += 1; new_degs[b] += 1

        print(f"    Quotient: V={V_new}, E={E_new}")
        print(f"    Degree sequence: {sorted(new_degs)}")
        print(f"    3-regular? {all(d == 3 for d in new_degs)}")

        # Check Petersen
        if V_new == 10 and E_new == 15:
            print(f"    *** COULD BE PETERSEN! ***")
        elif V_new == 5 and E_new == 5:
            print(f"    *** COULD BE C_5 (pentagon)! ***")

        # Node contents
        for r in roots:
            members = [i for i in range(V) if equiv[i] == r]
            h_vals = [mi[m]['H'] for m in members]
            print(f"    Supernode {new_id[r]}: nodes {members}, H={h_vals}, "
                  f"degree={new_degs[new_id[r]]}")

    # Also compute Jaccard for the WIGGLY-ONLY and COMPLEMENT-ONLY subgraphs
    # (if we have the data)

    # Find near-twins (Jaccard > 0.8)
    near_twins = [(i,j,jacc) for i,j,jacc,_,_ in jaccard_pairs if 0.8 <= jacc < 1.0]
    if near_twins:
        print(f"\n  NEAR-TWINS (Jaccard > 0.8): {len(near_twins)}")
        for i, j, jacc in near_twins[:10]:
            print(f"    ({i},{j}): Jaccard={jacc:.3f}, H=({mi[i]['H']},{mi[j]['H']})")

for n in [5, 6, 7]:
    analyze_jaccard(n)

print("\nDONE.")
