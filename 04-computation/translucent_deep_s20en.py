#!/usr/bin/env python3
"""
translucent_deep_s20en.py — Deep analysis of translucent subgraph
kind-pasteur-2026-03-23-S20en

Apply ALL metagraph analysis methods to the translucent (overlap-only) subgraph:
  - Edge decomposition (translucent-only / opaque-only / both)
  - SC-SC / SC-NS / NS-NS classification (spine/ribs/sea)
  - H-gradient (uphill/downhill/level)
  - Degree sequences
  - Lines per edge (thickness)
  - Self-loops per class
  - Connectivity, diameter
  - Spectral data (adjacency eigenvalues)
  - Blue/black tiling classification
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TRANSLUCENT DEEP ANALYSIS")
print("  kind-pasteur-2026-03-23-S20en")
print("=" * 80)

for n in range(4, 7):
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    perms = list(permutations(range(n)))

    # Classify arcs
    overlap_set = set()
    boundary_set = set()
    for k, (i,j) in enumerate(ALL_ARCS):
        if i >= 1 and j <= n-2:
            overlap_set.add(k)
        else:
            boundary_set.add(k)

    def canon(bits):
        best = None
        for p in perms:
            nb = 0
            for k, (i,j) in enumerate(ALL_ARCS):
                pi, pj = p[i], p[j]
                if pi < pj:
                    nk = ALL_ARCS.index((pi, pj))
                    if bits & (1 << k): nb |= (1 << nk)
                else:
                    nk = ALL_ARCS.index((pj, pi))
                    if not (bits & (1 << k)): nb |= (1 << nk)
            if best is None or nb < best: best = nb
        return best

    # Build canonical map and class data
    print(f"\n{'#'*70}")
    print(f"  n = {n}")
    print(f"  Overlap: {len(overlap_set)} arcs, Boundary: {len(boundary_set)} arcs")
    print(f"{'#'*70}")

    canon_map = {}
    classes = defaultdict(list)
    for bits in range(1 << m):
        cn = canon(bits)
        canon_map[bits] = cn
        classes[cn].append(bits)

    # Compute class properties: H, |Aut|, SC/NS, score
    class_info = {}
    for cn in classes:
        rep = cn
        # |Aut|
        aut = 0
        for p in perms:
            nb = 0
            for k, (i,j) in enumerate(ALL_ARCS):
                pi, pj = p[i], p[j]
                if pi < pj:
                    nk = ALL_ARCS.index((pi, pj))
                    if rep & (1 << k): nb |= (1 << nk)
                else:
                    nk = ALL_ARCS.index((pj, pi))
                    if not (rep & (1 << k)): nb |= (1 << nk)
            if nb == rep: aut += 1

        # Complement
        comp_bits = ((1 << m) - 1) ^ rep
        comp_cn = canon(comp_bits)
        is_sc = (comp_cn == cn)

        # Score sequence
        adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(ALL_ARCS):
            if rep & (1 << k): adj[i][j] = 1
            else: adj[j][i] = 1
        scores = tuple(sorted(sum(adj[i]) for i in range(n)))

        # H (number of Hamiltonian paths) - approximate via |class| * |Aut|
        H = len(classes[cn]) * aut  # This is actually tilings, not H

        class_info[cn] = {
            'aut': aut, 'sc': is_sc, 'comp': comp_cn,
            'scores': scores, 'size': len(classes[cn]),
            'tilings': len(classes[cn])
        }

    # Merge complement pairs for G_n/Z_2
    merged_classes = {}
    seen = set()
    for cn in classes:
        if cn in seen: continue
        comp = class_info[cn]['comp']
        merged_cn = min(cn, comp)
        merged_classes[merged_cn] = {
            'sc': class_info[cn]['sc'],
            'members': [cn] if cn == comp else [cn, comp],
            'scores': class_info[cn]['scores'],
            'aut': class_info[cn]['aut'],
            'tilings': class_info[cn]['tilings'] + (0 if cn == comp else class_info[comp]['tilings'])
        }
        seen.add(cn)
        seen.add(comp)

    V_merged = len(merged_classes)
    # Map each unmerged class to its merged class
    unmerged_to_merged = {}
    for mcn, minfo in merged_classes.items():
        for ucn in minfo['members']:
            unmerged_to_merged[ucn] = mcn

    print(f"  V = {len(classes)}, V_merged = {V_merged}, SC = {sum(1 for v in merged_classes.values() if v['sc'])}")

    # Build edge sets by type
    trans_edges = defaultdict(int)   # merged edge -> translucent line count
    opaque_edges = defaultdict(int)  # merged edge -> opaque line count
    trans_sl = defaultdict(int)      # merged class -> translucent self-loop lines
    opaque_sl = defaultdict(int)     # merged class -> opaque self-loop lines

    for bits in range(1 << m):
        cn1 = canon_map[bits]
        mcn1 = unmerged_to_merged[cn1]

        for k in range(m):
            flipped = bits ^ (1 << k)
            cn2 = canon_map[flipped]
            mcn2 = unmerged_to_merged[cn2]
            is_trans = k in overlap_set

            if mcn1 == mcn2:
                # Self-loop in merged graph
                if is_trans: trans_sl[mcn1] += 1
                else: opaque_sl[mcn1] += 1
            else:
                edge = (min(mcn1,mcn2), max(mcn1,mcn2))
                if is_trans: trans_edges[edge] += 1
                else: opaque_edges[edge] += 1

    # All edges
    all_edges = set(trans_edges.keys()) | set(opaque_edges.keys())
    trans_only = set(trans_edges.keys()) - set(opaque_edges.keys())
    opaque_only = set(opaque_edges.keys()) - set(trans_edges.keys())
    both_edges = set(trans_edges.keys()) & set(opaque_edges.keys())

    E_merged = len(all_edges)
    E_trans = len(trans_edges)
    E_opaque = len(opaque_edges)

    print(f"\n  EDGE DECOMPOSITION (merged metagraph):")
    print(f"    Total edges:         {E_merged}")
    print(f"    Translucent edges:   {E_trans} ({E_trans/E_merged*100:.1f}%)")
    print(f"    Opaque edges:        {E_opaque} ({E_opaque/E_merged*100:.1f}%)")
    print(f"    Trans-only:          {len(trans_only)}")
    print(f"    Opaque-only:         {len(opaque_only)}")
    print(f"    Both:                {len(both_edges)}")

    # SC-SC / SC-NS / NS-NS classification
    def edge_type(e):
        sc1 = merged_classes[e[0]]['sc']
        sc2 = merged_classes[e[1]]['sc']
        if sc1 and sc2: return 'SC-SC'
        elif sc1 or sc2: return 'SC-NS'
        else: return 'NS-NS'

    type_counts = {'SC-SC': {'trans': 0, 'opaque': 0, 'both': 0, 'total': 0},
                   'SC-NS': {'trans': 0, 'opaque': 0, 'both': 0, 'total': 0},
                   'NS-NS': {'trans': 0, 'opaque': 0, 'both': 0, 'total': 0}}

    for e in all_edges:
        et = edge_type(e)
        type_counts[et]['total'] += 1
        if e in trans_only: type_counts[et]['trans'] += 1
        elif e in opaque_only: type_counts[et]['opaque'] += 1
        else: type_counts[et]['both'] += 1

    print(f"\n  SPINE/RIBS/SEA DECOMPOSITION:")
    print(f"    {'Type':>8} {'Total':>6} {'Trans-only':>11} {'Opaque-only':>12} {'Both':>6}")
    for et in ['SC-SC', 'SC-NS', 'NS-NS']:
        tc = type_counts[et]
        print(f"    {et:>8} {tc['total']:6d} {tc['trans']:11d} {tc['opaque']:12d} {tc['both']:6d}")

    # Lines per edge (thickness)
    trans_thickness = [trans_edges[e] for e in trans_edges]
    opaque_thickness = [opaque_edges[e] for e in opaque_edges]
    total_thickness = [trans_edges.get(e,0) + opaque_edges.get(e,0) for e in all_edges]

    print(f"\n  LINES PER EDGE (thickness):")
    if trans_thickness:
        print(f"    Translucent: min={min(trans_thickness)}, max={max(trans_thickness)}, avg={sum(trans_thickness)/len(trans_thickness):.1f}")
    if opaque_thickness:
        print(f"    Opaque:      min={min(opaque_thickness)}, max={max(opaque_thickness)}, avg={sum(opaque_thickness)/len(opaque_thickness):.1f}")
    print(f"    Total:       min={min(total_thickness)}, max={max(total_thickness)}, avg={sum(total_thickness)/len(total_thickness):.1f}")

    # Translucent / opaque line ratio per edge
    ratios = []
    for e in both_edges:
        t = trans_edges[e]
        o = opaque_edges[e]
        ratios.append(t / (t + o))
    if ratios:
        print(f"    Trans/(Trans+Opaque) ratio per edge: min={min(ratios):.3f}, max={max(ratios):.3f}, avg={sum(ratios)/len(ratios):.3f}")
        expected_ratio = len(overlap_set) / m
        print(f"    Expected if uniform: {expected_ratio:.3f} = {len(overlap_set)}/{m}")

    # Degree sequences
    trans_degree = Counter()
    opaque_degree = Counter()
    total_degree = Counter()
    for e in all_edges:
        a, b = e
        if e in trans_edges or e in trans_only or e in both_edges:
            if e in set(trans_edges.keys()):
                trans_degree[a] += 1
                trans_degree[b] += 1
        if e in opaque_edges or e in opaque_only or e in both_edges:
            if e in set(opaque_edges.keys()):
                opaque_degree[a] += 1
                opaque_degree[b] += 1
        total_degree[a] += 1
        total_degree[b] += 1

    print(f"\n  DEGREE SEQUENCES (merged):")
    td = sorted(total_degree.values())
    trd = sorted(trans_degree.values()) if trans_degree else [0]
    opd = sorted(opaque_degree.values()) if opaque_degree else [0]
    print(f"    Total:       min={td[0]}, max={td[-1]}, avg={sum(td)/len(td):.1f}")
    print(f"    Translucent: min={trd[0]}, max={trd[-1]}, avg={sum(trd)/len(trd):.1f}")
    print(f"    Opaque:      min={opd[0]}, max={opd[-1]}, avg={sum(opd)/len(opd):.1f}")

    # Self-loops per class
    trans_sl_list = [trans_sl.get(mcn, 0) for mcn in merged_classes]
    opaque_sl_list = [opaque_sl.get(mcn, 0) for mcn in merged_classes]
    total_sl_list = [trans_sl.get(mcn,0) + opaque_sl.get(mcn,0) for mcn in merged_classes]

    print(f"\n  SELF-LOOPS PER CLASS:")
    print(f"    Classes with trans SL:  {sum(1 for x in trans_sl_list if x > 0)}/{V_merged}")
    print(f"    Classes with opaque SL: {sum(1 for x in opaque_sl_list if x > 0)}/{V_merged}")
    print(f"    Total trans SL lines:   {sum(trans_sl_list)}")
    print(f"    Total opaque SL lines:  {sum(opaque_sl_list)}")

    # Connectivity of translucent subgraph
    trans_adj = defaultdict(set)
    for e in trans_edges:
        trans_adj[e[0]].add(e[1])
        trans_adj[e[1]].add(e[0])

    # BFS for connected components
    visited = set()
    components = []
    for v in merged_classes:
        if v in visited: continue
        comp = set()
        queue = [v]
        while queue:
            u = queue.pop()
            if u in comp: continue
            comp.add(u)
            for w in trans_adj.get(u, set()):
                if w not in comp:
                    queue.append(w)
        visited |= comp
        components.append(comp)

    # Same for opaque
    opaque_adj = defaultdict(set)
    for e in opaque_edges:
        opaque_adj[e[0]].add(e[1])
        opaque_adj[e[1]].add(e[0])

    visited_o = set()
    components_o = []
    for v in merged_classes:
        if v in visited_o: continue
        comp = set()
        queue = [v]
        while queue:
            u = queue.pop()
            if u in comp: continue
            comp.add(u)
            for w in opaque_adj.get(u, set()):
                if w not in comp:
                    queue.append(w)
        visited_o |= comp
        components_o.append(comp)

    print(f"\n  CONNECTIVITY:")
    print(f"    Full graph:       {1} component (V={V_merged})")
    print(f"    Translucent:      {len(components)} components (sizes: {sorted([len(c) for c in components], reverse=True)[:5]})")
    print(f"    Opaque:           {len(components_o)} components (sizes: {sorted([len(c) for c in components_o], reverse=True)[:5]})")

    # Diameter (BFS from each vertex)
    def diameter(adj, verts):
        max_d = 0
        for start in verts:
            dist = {start: 0}
            queue = [start]
            while queue:
                u = queue.pop(0)
                for w in adj.get(u, set()):
                    if w not in dist:
                        dist[w] = dist[u] + 1
                        queue.append(w)
            if len(dist) == len(verts):
                max_d = max(max_d, max(dist.values()))
        return max_d if max_d > 0 else float('inf')

    full_adj = defaultdict(set)
    for e in all_edges:
        full_adj[e[0]].add(e[1])
        full_adj[e[1]].add(e[0])

    if V_merged <= 500:
        d_full = diameter(full_adj, set(merged_classes.keys()))
        # Only compute trans diameter if connected
        if len(components) == 1:
            d_trans = diameter(trans_adj, set(merged_classes.keys()))
        else:
            d_trans = float('inf')
        d_opaque = diameter(opaque_adj, set(merged_classes.keys())) if len(components_o) == 1 else float('inf')

        print(f"\n  DIAMETER:")
        print(f"    Full:         {d_full}")
        print(f"    Translucent:  {d_trans}")
        print(f"    Opaque:       {d_opaque}")

    # Spectral data for small n
    if V_merged <= 100:
        try:
            import numpy as np
            verts = sorted(merged_classes.keys())
            v_idx = {v: i for i, v in enumerate(verts)}
            N = len(verts)

            # Full adjacency
            A_full = np.zeros((N, N))
            A_trans = np.zeros((N, N))
            A_opaque = np.zeros((N, N))

            for e in all_edges:
                i, j = v_idx[e[0]], v_idx[e[1]]
                A_full[i,j] = A_full[j,i] = 1
            for e in trans_edges:
                i, j = v_idx[e[0]], v_idx[e[1]]
                A_trans[i,j] = A_trans[j,i] = 1
            for e in opaque_edges:
                i, j = v_idx[e[0]], v_idx[e[1]]
                A_opaque[i,j] = A_opaque[j,i] = 1

            eig_full = sorted(np.linalg.eigvalsh(A_full), reverse=True)
            eig_trans = sorted(np.linalg.eigvalsh(A_trans), reverse=True)
            eig_opaque = sorted(np.linalg.eigvalsh(A_opaque), reverse=True)

            print(f"\n  SPECTRAL DATA:")
            print(f"    Full spectrum (top 5):       {[f'{x:.2f}' for x in eig_full[:5]]}")
            print(f"    Translucent spectrum (top 5): {[f'{x:.2f}' for x in eig_trans[:5]]}")
            print(f"    Opaque spectrum (top 5):      {[f'{x:.2f}' for x in eig_opaque[:5]]}")
            print(f"    Full spectral gap:       {eig_full[0] - eig_full[1]:.4f}")
            if eig_trans[0] > 0:
                print(f"    Translucent spectral gap: {eig_trans[0] - eig_trans[1]:.4f}")
            if eig_opaque[0] > 0:
                print(f"    Opaque spectral gap:      {eig_opaque[0] - eig_opaque[1]:.4f}")

            # Check: A_full = A_trans + A_opaque?
            diff = np.max(np.abs(A_full - A_trans - A_opaque))
            print(f"\n    A_full = A_trans + A_opaque? max|diff| = {diff:.6f}")

        except ImportError:
            print("  (numpy not available for spectral analysis)")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
