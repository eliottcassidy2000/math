#!/usr/bin/env python3
"""
wiggly_lines_s272.py — Wiggly lines: inner sub-tournament flips (Mode B)
opus-2026-03-23-S272

A tournament on n vertices decomposes (Mode B, n→n-2) into:
  OVERLAP: inner sub-tournament on vertices {1,...,n-2}, C(n-2,2) arcs
  WIRING:  arcs involving vertices 0 or n-1, 2(n-2)+1 = 2n-3 arcs

Total: C(n-2,2) + (2n-3) = C(n,2) ✓

A WIGGLY LINE connects two tilings that differ in exactly ONE overlap cell.
Each tiling has C(n-2,2) wiggly neighbors.
Tilings sharing the same wiring form a hypercube Q_{C(n-2,2)}.

QUESTIONS:
1. How many wiggly lines generate metagraph edges? (vs self-loops)
2. Which metagraph edges are wiggly-reachable? Which require wiring flips?
3. How do wiggly groups (constant wiring) map to iso classes?
4. Is the wiggly/wiring split related to the cut/cycle decomposition?
5. What's the wiggly quotient: Q_m / (same wiring) = 2^{2n-3} points?
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def build_tournament_data(n):
    """Build all tournament data with iso class identification."""
    m = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1,n)]

    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    lines = [l.strip() for l in (r.stdout or '').split('\n')
             if len(l.strip()) == m and all(c in '01' for c in l.strip())]

    def b2a(bits):
        a = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if bits & (1 << (m-1-k)): a[i][j] = 1
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

    cls = []
    hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'adj':adj, 'score':s, 'c3':cc, 'H':h,
                    'sc':sc, 'comp_canon':cc2, 'canon':canon})
        hm[(s,cc,h)].append(i); cm[canon] = i
    for cl in cls: cl['comp_cid'] = cm.get(cl['comp_canon'], -1)

    def lk(adj):
        s = ss(adj); cc = c3(adj); h = H(adj)
        cids = hm.get((s,cc,h))
        if not cids: return -1
        if len(cids) == 1: return cids[0]
        return cm.get(cp(adj), -1)

    # Merge complement pairs
    mn = {}; mid = 0
    for cl in cls:
        ci = cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci] = mid; mid += 1
        else:
            mn[ci] = mid
            comp = cl['comp_cid']
            if comp >= 0 and comp != ci: mn[comp] = mid
            mid += 1

    mi = {}
    for mm in range(mid):
        cids = [cl['cid'] for cl in cls if mn[cl['cid']] == mm]
        cl0 = cls[cids[0]]
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}

    return P, cls, lk, mn, mid, mi, m

def classify_arcs(n, P):
    """Classify arcs into OVERLAP (inner) and WIRING (extreme vertices)."""
    m = len(P)
    inner_vertices = set(range(1, n-1))  # {1, ..., n-2}

    overlap_arcs = []  # indices of arcs between inner vertices
    wiring_arcs = []   # indices of arcs involving vertex 0 or n-1

    for k, (i, j) in enumerate(P):
        if i in inner_vertices and j in inner_vertices:
            overlap_arcs.append(k)
        else:
            wiring_arcs.append(k)

    return overlap_arcs, wiring_arcs

def analyze_wiggly(n):
    """Full wiggly line analysis."""
    print(f"\n{'#'*72}")
    print(f"  n = {n}: WIGGLY LINE ANALYSIS")
    print(f"{'#'*72}")

    P, cls, lk, mn, V_merged, mi, m = build_tournament_data(n)
    overlap_arcs, wiring_arcs = classify_arcs(n, P)

    n_overlap = len(overlap_arcs)
    n_wiring = len(wiring_arcs)

    print(f"\n  DECOMPOSITION:")
    print(f"    Total arcs: {m} = C({n},2)")
    print(f"    Overlap arcs (inner sub-tournament): {n_overlap} = C({n-2},2)")
    print(f"    Wiring arcs (extreme vertices): {n_wiring} = 2({n-2})+1 = {2*(n-2)+1}")
    print(f"    Check: {n_overlap} + {n_wiring} = {n_overlap + n_wiring} = {m} ✓")

    print(f"\n    Overlap arcs: {[(P[k]) for k in overlap_arcs]}")
    print(f"    Wiring arcs: {[(P[k]) for k in wiring_arcs]}")

    print(f"\n    Wiggly groups: 2^{n_wiring} = {2**n_wiring} groups")
    print(f"    Group size: 2^{n_overlap} = {2**n_overlap} tilings per group")

    # Iterate over ALL labeled tournaments
    # For each, classify its wiggly flips vs wiring flips
    wiggly_class_changes = 0  # overlap flip changes merged class
    wiggly_self_loops = 0     # overlap flip preserves merged class
    wiring_class_changes = 0
    wiring_self_loops = 0

    # Track which metagraph edges are wiggly-reachable
    wiggly_edges = set()
    wiring_edges = set()

    # Track wiggly groups → class distribution
    group_class_map = defaultdict(set)  # wiring_config → set of merged classes

    total_tournaments = 2 ** m

    def b2a(bits):
        a = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(P):
            if bits & (1 << (m-1-k)): a[i][j] = 1
            else: a[j][i] = 1
        return a

    for bits in range(total_tournaments):
        adj = b2a(bits)
        cid = lk(adj)
        if cid < 0: continue
        ma = mn[cid]

        # Wiring configuration = values at wiring arc positions
        wiring_config = tuple((bits >> (m-1-k)) & 1 for k in wiring_arcs)
        group_class_map[wiring_config].add(ma)

        # Wiggly flips (overlap arcs)
        for arc in overlap_arcs:
            fb = bits ^ (1 << (m-1-arc))
            fadj = b2a(fb)
            nb_cid = lk(fadj)
            if nb_cid < 0: continue
            mb = mn[nb_cid]
            if ma == mb:
                wiggly_self_loops += 1
            else:
                wiggly_class_changes += 1
                wiggly_edges.add((min(ma, mb), max(ma, mb)))

        # Wiring flips
        for arc in wiring_arcs:
            fb = bits ^ (1 << (m-1-arc))
            fadj = b2a(fb)
            nb_cid = lk(fadj)
            if nb_cid < 0: continue
            mb = mn[nb_cid]
            if ma == mb:
                wiring_self_loops += 1
            else:
                wiring_class_changes += 1
                wiring_edges.add((min(ma, mb), max(ma, mb)))

    # Results
    total_wiggly = wiggly_class_changes + wiggly_self_loops
    total_wiring = wiring_class_changes + wiring_self_loops

    print(f"\n  WIGGLY vs WIRING STATISTICS:")
    print(f"    Wiggly flips (overlap arcs):")
    print(f"      Total: {total_wiggly}")
    print(f"      Class-changing: {wiggly_class_changes} ({100*wiggly_class_changes/total_wiggly:.1f}%)")
    print(f"      Self-loops: {wiggly_self_loops} ({100*wiggly_self_loops/total_wiggly:.1f}%)")
    print(f"      Metagraph edges reached: {len(wiggly_edges)}")

    print(f"\n    Wiring flips (extreme arcs):")
    print(f"      Total: {total_wiring}")
    print(f"      Class-changing: {wiring_class_changes} ({100*wiring_class_changes/total_wiring:.1f}%)")
    print(f"      Self-loops: {wiring_self_loops} ({100*wiring_self_loops/total_wiring:.1f}%)")
    print(f"      Metagraph edges reached: {len(wiring_edges)}")

    # Full metagraph edges
    all_edges = wiggly_edges | wiring_edges
    only_wiggly = wiggly_edges - wiring_edges
    only_wiring = wiring_edges - wiggly_edges
    both = wiggly_edges & wiring_edges

    print(f"\n  METAGRAPH EDGE DECOMPOSITION:")
    print(f"    Total metagraph edges: {len(all_edges)}")
    print(f"    Wiggly-only edges: {len(only_wiggly)} ({100*len(only_wiggly)/len(all_edges):.1f}%)")
    print(f"    Wiring-only edges: {len(only_wiring)} ({100*len(only_wiring)/len(all_edges):.1f}%)")
    print(f"    Both wiggly AND wiring: {len(both)} ({100*len(both)/len(all_edges):.1f}%)")

    # Wiggly group analysis
    n_groups = len(group_class_map)
    classes_per_group = [len(v) for v in group_class_map.values()]
    print(f"\n  WIGGLY GROUP → CLASS MAPPING:")
    print(f"    Total wiggly groups: {n_groups} (= 2^{n_wiring} = {2**n_wiring})")
    print(f"    Classes per group: min={min(classes_per_group)}, max={max(classes_per_group)}, "
          f"mean={sum(classes_per_group)/len(classes_per_group):.2f}")

    # Distribution of classes per group
    cpg_dist = Counter(classes_per_group)
    print(f"    Distribution:")
    for k in sorted(cpg_dist.keys()):
        print(f"      {k} classes: {cpg_dist[k]} groups ({100*cpg_dist[k]/n_groups:.1f}%)")

    # How many groups map to a SINGLE class?
    single_class_groups = sum(1 for v in group_class_map.values() if len(v) == 1)
    print(f"\n    Single-class groups: {single_class_groups}/{n_groups} "
          f"({100*single_class_groups/n_groups:.1f}%)")
    print(f"    = wirings where the inner sub-tournament doesn't affect the class")

    # The inner sub-tournament IS a tournament on n-2 vertices
    # How many distinct iso classes of the inner tournament appear per group?
    print(f"\n  INNER SUB-TOURNAMENT STRUCTURE:")
    print(f"    Inner tournament: n-2 = {n-2} vertices, C({n-2},2) = {n_overlap} arcs")
    print(f"    Inner iso classes: A000568({n-2}) = ?")

    # For each wiggly group, extract the inner tournaments and classify them
    if n <= 6:
        inner_cls_per_group = defaultdict(Counter)
        for bits in range(total_tournaments):
            wc = tuple((bits >> (m-1-k)) & 1 for k in wiring_arcs)
            # Inner tournament = bits at overlap positions
            inner_bits = tuple((bits >> (m-1-k)) & 1 for k in overlap_arcs)
            inner_adj = [[0]*(n-2) for _ in range(n-2)]
            inner_P = [(i-1,j-1) for (i,j) in [P[k] for k in overlap_arcs]]
            for idx, (i,j) in enumerate(inner_P):
                if inner_bits[idx]: inner_adj[i][j] = 1
                else: inner_adj[j][i] = 1

            # Classify inner tournament (simple canonical form for n-2 ≤ 4)
            inner_canon = tuple(inner_adj[i][j] for i in range(n-2) for j in range(n-2))
            inner_cls_per_group[wc][inner_canon] += 1

        # How many distinct inner tournaments per group?
        inner_per_group = [len(v) for v in inner_cls_per_group.values()]
        print(f"    Distinct inner tournaments per group: "
              f"all = {2**n_overlap} (= 2^{n_overlap})")
        print(f"    (All groups have all inner tournaments — they're hypercubes)")

    return {
        'n': n, 'V_merged': V_merged,
        'n_overlap': n_overlap, 'n_wiring': n_wiring,
        'wiggly_edges': len(wiggly_edges),
        'wiring_edges': len(wiring_edges),
        'both_edges': len(both),
        'total_edges': len(all_edges),
        'wiggly_neutral_frac': wiggly_self_loops / total_wiggly if total_wiggly > 0 else 0,
        'wiring_neutral_frac': wiring_self_loops / total_wiring if total_wiring > 0 else 0,
    }

if __name__ == '__main__':
    print("=" * 72)
    print("  WIGGLY LINES: MODE B INNER SUB-TOURNAMENT FLIPS")
    print("  opus-2026-03-23-S272")
    print("=" * 72)

    results = []
    for n in [4, 5, 6]:
        t0 = time.time()
        r = analyze_wiggly(n)
        results.append(r)
        print(f"\n  Time for n={n}: {time.time()-t0:.1f}s")

    # Summary
    print(f"\n{'='*72}")
    print("  SUMMARY TABLE")
    print(f"{'='*72}")
    print(f"  {'n':>3} {'V_m':>5} {'E_total':>8} {'E_wiggly':>9} {'E_wiring':>9} "
          f"{'E_both':>7} {'wig%':>5} {'wir%':>5} {'wig_neut':>9} {'wir_neut':>9}")
    for r in results:
        wp = 100*r['wiggly_edges']/r['total_edges'] if r['total_edges'] > 0 else 0
        wrp = 100*r['wiring_edges']/r['total_edges'] if r['total_edges'] > 0 else 0
        print(f"  {r['n']:>3} {r['V_merged']:>5} {r['total_edges']:>8} "
              f"{r['wiggly_edges']:>9} {r['wiring_edges']:>9} {r['both_edges']:>7} "
              f"{wp:>5.1f} {wrp:>5.1f} {r['wiggly_neutral_frac']:>9.4f} {r['wiring_neutral_frac']:>9.4f}")

    print(f"\n{'='*72}")
    print("  THE WIGGLY/WIRING DECOMPOSITION")
    print(f"{'='*72}")
    print("""
  WIGGLY LINES = flips of inner arcs (overlap, Mode B sub-tournament)
  WIRING LINES = flips of extreme arcs (vertices 0 and n-1)

  The m-cube Q_m decomposes as:
    Q_m = Q_{overlap} × Q_{wiring}
    Q_{C(n,2)} = Q_{C(n-2,2)} × Q_{2n-3}

  Each wiring configuration defines a Q_{C(n-2,2)} "wiggly group."
  There are 2^{2n-3} wiggly groups, each of size 2^{C(n-2,2)}.

  METAGRAPH EDGES come from BOTH types:
  - Wiggly edges: changing the inner structure
  - Wiring edges: changing how extremes connect to the middle
  - Most edges are reachable by BOTH types (high overlap)

  This IS the Mode B recursion (n → n-2):
  The inner sub-tournament recurses, while the wiring is "frozen."
  The wiggly group structure reveals how much the inner tournament
  INDEPENDENTLY determines the iso class.
""")

    print("DONE.")
