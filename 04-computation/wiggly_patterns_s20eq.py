#!/usr/bin/env python3
"""
wiggly_patterns_s20eq.py — Hunt for patterns in wiggly line counts vs metagraph
kind-pasteur-2026-03-23-S20eq

For each merged class C and each of its metagraph neighbors C':
  - How many wiggly lines connect tilings in C to tilings in C'?
  - How many opaque lines?
  - What is the wiggly/opaque ratio per edge?
  - Does it correlate with H, |Aut|, SC/NS, spine/rib/sea type?

For each merged class C's self-loops:
  - How many wiggly self-loop lines vs opaque self-loop lines?
  - Does the wiggly self-loop count predict anything about C?

Also: the FIBER PROFILE of each class:
  - Across all fibers, how is class C distributed?
  - How many fibers contain at least one tiling in C?
  - What fraction of each fiber belongs to C?
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  WIGGLY PATTERN HUNTING")
print("  kind-pasteur-2026-03-23-S20eq")
print("=" * 80)

for n in range(4, 7):
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    perms = list(permutations(range(n)))

    overlap_indices = []
    boundary_indices = []
    for k, (i,j) in enumerate(ALL_ARCS):
        if i >= 1 and j <= n-2:
            overlap_indices.append(k)
        else:
            boundary_indices.append(k)

    n_overlap = len(overlap_indices)
    n_boundary = len(boundary_indices)

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

    canon_map = {}
    for bits in range(1 << m):
        canon_map[bits] = canon(bits)

    # Merged classes
    comp_map = {}
    for cn in set(canon_map.values()):
        comp_map[cn] = canon(((1 << m) - 1) ^ cn)

    merged_map = {}
    for cn in set(canon_map.values()):
        merged_map[cn] = min(cn, comp_map[cn])

    # Class info
    class_members = defaultdict(list)
    for bits in range(1 << m):
        cn = canon_map[bits]
        class_members[cn].append(bits)

    merged_classes = sorted(set(merged_map.values()))
    V_merged = len(merged_classes)

    # Compute |Aut| and score for each class
    class_info = {}
    for cn in set(canon_map.values()):
        aut = sum(1 for p in perms if canon(sum((((cn >> k) & 1) << k) for k in range(m))) == cn
                  ) if False else 0  # placeholder
        # Simpler: count members
        size = len(class_members[cn])
        aut = factorial(n) // size
        is_sc = (comp_map[cn] == cn)
        # Score sequence
        adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(ALL_ARCS):
            if cn & (1 << k): adj[i][j] = 1
            else: adj[j][i] = 1
        scores = tuple(sorted(sum(adj[i]) for i in range(n)))
        # H (tilings * |Aut|)
        tilings = size
        H = tilings * aut
        class_info[cn] = {'aut': aut, 'sc': is_sc, 'scores': scores,
                          'size': size, 'H': H, 'tilings': tilings}

    merged_info = {}
    for mcn in merged_classes:
        members = [cn for cn in set(canon_map.values()) if merged_map[cn] == mcn]
        ci = class_info[members[0]]
        total_tilings = sum(class_info[cn]['tilings'] for cn in members)
        merged_info[mcn] = {
            'aut': ci['aut'], 'sc': ci['sc'], 'scores': ci['scores'],
            'H': ci['H'], 'tilings': total_tilings
        }

    # Count wiggly and opaque lines per edge and self-loop
    edge_wiggly = defaultdict(int)  # (mcn1, mcn2) -> count
    edge_opaque = defaultdict(int)
    sl_wiggly = defaultdict(int)    # mcn -> count
    sl_opaque = defaultdict(int)

    # Also count per-class total wiggly/opaque TRANSITIONS (directed)
    class_wiggly_out = defaultdict(int)
    class_opaque_out = defaultdict(int)
    class_wiggly_sl = defaultdict(int)
    class_opaque_sl = defaultdict(int)

    # Fiber profile: for each class, in how many fibers does it appear?
    class_fiber_count = defaultdict(int)  # mcn -> number of fibers with >= 1 tiling in mcn
    class_fiber_tilings = defaultdict(list)  # mcn -> [count per fiber]

    for bb in range(1 << n_boundary):
        base_bits = 0
        for idx, k in enumerate(boundary_indices):
            if bb & (1 << idx):
                base_bits |= (1 << k)

        fiber_class_count = defaultdict(int)
        fiber_tilings = {}

        for ob in range(1 << n_overlap):
            full_bits = base_bits
            for idx, k in enumerate(overlap_indices):
                if ob & (1 << idx):
                    full_bits |= (1 << k)
            cn = canon_map[full_bits]
            mcn = merged_map[cn]
            fiber_class_count[mcn] += 1
            fiber_tilings[ob] = mcn

        for mcn, count in fiber_class_count.items():
            class_fiber_count[mcn] += 1
            class_fiber_tilings[mcn].append(count)

        # Count wiggly lines within this fiber
        for ob in range(1 << n_overlap):
            mcn1 = fiber_tilings[ob]
            for idx in range(n_overlap):
                ob2 = ob ^ (1 << idx)
                mcn2 = fiber_tilings[ob2]
                if mcn1 <= mcn2:  # avoid double counting
                    if mcn1 == mcn2:
                        sl_wiggly[mcn1] += 1
                        class_wiggly_sl[mcn1] += 1
                    else:
                        edge = (mcn1, mcn2)
                        edge_wiggly[edge] += 1
                        class_wiggly_out[mcn1] += 1
                        class_wiggly_out[mcn2] += 1

    # Opaque lines (boundary flips)
    for bits in range(1 << m):
        cn = canon_map[bits]
        mcn1 = merged_map[cn]
        for k in boundary_indices:
            flipped = bits ^ (1 << k)
            cn2 = canon_map[flipped]
            mcn2 = merged_map[cn2]
            if mcn1 <= mcn2:
                if mcn1 == mcn2:
                    sl_opaque[mcn1] += 1
                    class_opaque_sl[mcn1] += 1
                else:
                    edge = (mcn1, mcn2)
                    edge_opaque[edge] += 1
                    class_opaque_out[mcn1] += 1
                    class_opaque_out[mcn2] += 1

    all_edges = set(edge_wiggly.keys()) | set(edge_opaque.keys())

    print(f"\n{'#'*70}")
    print(f"  n = {n}, V_merged = {V_merged}, E = {len(all_edges)}")
    print(f"{'#'*70}")

    # ================================================================
    # PATTERN 1: Wiggly/Total ratio per edge — is it REALLY constant?
    # ================================================================
    ratios = []
    for e in all_edges:
        w = edge_wiggly.get(e, 0)
        o = edge_opaque.get(e, 0)
        total = w + o
        if total > 0:
            ratios.append(w / total)

    expected = n_overlap / m
    print(f"\n  PATTERN 1: Wiggly fraction per edge")
    print(f"    Expected (uniform): {expected:.6f} = {n_overlap}/{m}")
    if ratios:
        print(f"    Actual: min={min(ratios):.6f}, max={max(ratios):.6f}")
        print(f"    All equal? {len(set(f'{r:.10f}' for r in ratios)) == 1}")

    # ================================================================
    # PATTERN 2: Wiggly lines per edge vs H difference
    # ================================================================
    print(f"\n  PATTERN 2: Wiggly lines per edge, sorted by edge properties")
    print(f"    {'Edge':>20} {'Type':>6} {'Wiggly':>8} {'Opaque':>8} {'Total':>8} {'dH':>6} {'H_sum':>8}")

    edge_data = []
    for e in sorted(all_edges):
        w = edge_wiggly.get(e, 0)
        o = edge_opaque.get(e, 0)
        H1 = merged_info[e[0]]['H']
        H2 = merged_info[e[1]]['H']
        sc1 = merged_info[e[0]]['sc']
        sc2 = merged_info[e[1]]['sc']
        if sc1 and sc2: etype = 'SC-SC'
        elif sc1 or sc2: etype = 'SC-NS'
        else: etype = 'NS-NS'
        dH = abs(H1 - H2)
        H_sum = H1 + H2
        edge_data.append((e, etype, w, o, w+o, dH, H_sum))

    edge_data.sort(key=lambda x: x[4])  # sort by total lines
    for e, etype, w, o, total, dH, H_sum in edge_data[:20]:
        print(f"    {str(e):>20} {etype:>6} {w:8d} {o:8d} {total:8d} {dH:6d} {H_sum:8d}")
    if len(edge_data) > 20:
        print(f"    ... ({len(edge_data)} total edges)")

    # Check: is total lines proportional to H_sum or something?
    totals = [x[4] for x in edge_data]
    h_sums = [x[6] for x in edge_data]
    if len(set(totals)) > 1 and len(set(h_sums)) > 1:
        # Correlation
        mean_t = sum(totals) / len(totals)
        mean_h = sum(h_sums) / len(h_sums)
        cov = sum((t-mean_t)*(h-mean_h) for t,h in zip(totals, h_sums)) / len(totals)
        var_t = sum((t-mean_t)**2 for t in totals) / len(totals)
        var_h = sum((h-mean_h)**2 for h in h_sums) / len(h_sums)
        if var_t > 0 and var_h > 0:
            corr = cov / (var_t**0.5 * var_h**0.5)
            print(f"\n    Correlation(total_lines, H_sum): {corr:.4f}")

    # Also check total lines vs tilings product
    tilings_prods = [merged_info[e[0]]['tilings'] * merged_info[e[1]]['tilings'] for e, *_ in edge_data]
    if len(set(totals)) > 1 and len(set(tilings_prods)) > 1:
        mean_t2 = sum(totals) / len(totals)
        mean_tp = sum(tilings_prods) / len(tilings_prods)
        cov2 = sum((t-mean_t2)*(tp-mean_tp) for t,tp in zip(totals, tilings_prods)) / len(totals)
        var_t2 = sum((t-mean_t2)**2 for t in totals) / len(totals)
        var_tp = sum((tp-mean_tp)**2 for tp in tilings_prods) / len(tilings_prods)
        if var_t2 > 0 and var_tp > 0:
            corr2 = cov2 / (var_t2**0.5 * var_tp**0.5)
            print(f"    Correlation(total_lines, tilings_product): {corr2:.4f}")

    # Check exact proportionality: total_lines / (tilings1 * tilings2)
    ratios2 = [total / (merged_info[e[0]]['tilings'] * merged_info[e[1]]['tilings'])
               for (e, _, _, _, total, _, _) in edge_data
               if merged_info[e[0]]['tilings'] > 0 and merged_info[e[1]]['tilings'] > 0]
    if ratios2:
        print(f"    total_lines / (tilings1 * tilings2): min={min(ratios2):.4f}, max={max(ratios2):.4f}")
        print(f"    Constant? {max(ratios2) - min(ratios2) < 0.0001}")

    # ================================================================
    # PATTERN 3: Self-loop wiggly counts per class
    # ================================================================
    print(f"\n  PATTERN 3: Self-loop structure per class")
    print(f"    {'Class':>8} {'SC':>3} {'H':>6} {'Aut':>4} {'Tilings':>8} {'W_SL':>8} {'O_SL':>8} {'W_SL/T':>8} {'O_SL/T':>8} {'#Fibers':>8} {'AvgFib':>8}")

    for mcn in sorted(merged_classes, key=lambda c: merged_info[c]['H']):
        mi = merged_info[mcn]
        wsl = sl_wiggly.get(mcn, 0)
        osl = sl_opaque.get(mcn, 0)
        nfibers = class_fiber_count.get(mcn, 0)
        avg_fib = sum(class_fiber_tilings.get(mcn, [0])) / max(nfibers, 1)
        print(f"    {mcn:8d} {'Y' if mi['sc'] else 'N':>3} {mi['H']:6d} {mi['aut']:4d} {mi['tilings']:8d} {wsl:8d} {osl:8d} {wsl/max(mi['tilings'],1):8.1f} {osl/max(mi['tilings'],1):8.1f} {nfibers:8d} {avg_fib:8.2f}")

    # Check: W_SL / tilings = constant?
    wsl_per_tiling = [sl_wiggly.get(mcn, 0) / max(merged_info[mcn]['tilings'], 1) for mcn in merged_classes]
    osl_per_tiling = [sl_opaque.get(mcn, 0) / max(merged_info[mcn]['tilings'], 1) for mcn in merged_classes]
    print(f"\n    W_SL/tilings range: {min(wsl_per_tiling):.2f} to {max(wsl_per_tiling):.2f}")
    print(f"    O_SL/tilings range: {min(osl_per_tiling):.2f} to {max(osl_per_tiling):.2f}")

    # ================================================================
    # PATTERN 4: Fiber profile — how many fibers contain each class?
    # ================================================================
    print(f"\n  PATTERN 4: Fiber profile")
    total_fibers = 2 ** n_boundary
    for mcn in sorted(merged_classes, key=lambda c: class_fiber_count.get(c, 0)):
        mi = merged_info[mcn]
        nf = class_fiber_count.get(mcn, 0)
        fib_sizes = class_fiber_tilings.get(mcn, [])
        if fib_sizes:
            print(f"    {mcn:8d}: in {nf:4d}/{total_fibers} fibers ({nf/total_fibers*100:5.1f}%), avg {sum(fib_sizes)/len(fib_sizes):.2f} tilings/fiber, H={mi['H']}, tilings={mi['tilings']}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

print("\nDONE.")
print("=" * 80)
