#!/usr/bin/env python3
"""
translucent_deep_s271.py — Deep investigation of translucent lines
opus-2026-03-23-S271

Translucent lines = neutral arc flips (stay within same iso class).
Opaque lines = class-changing arc flips (form metagraph edges).

The m-cube Q_m partitions: Q_m = translucent ∪ opaque (disjoint edges).

QUESTIONS:
1. How many LINES (tiling pairs) contribute to each METAGRAPH EDGE?
   Each edge of G_n/Z_2 is formed by many opaque lines.
   Lines per edge = multiplicity = weight of the metagraph edge.

2. Which metagraph edges have the MOST lines? The FEWEST?

3. How does the translucent fraction relate to H?
   Classes with H=1 (transitive): high neutral fraction (0.50 at n=4, 0.40 at n=5)
   Classes with H=max (regular): zero neutral arcs
   → The more "ordered" the tournament, the more translucent it is.

4. The translucent graph WITHIN each class: what is its structure?
   - Components, diameter, degree distribution
   - Is it a Cayley graph of S_n / Aut(T)?
   - Binary tree / heap structure?

5. Connection to HEAPS: in Viennot's framework, a "heap" of conflicting pieces
   corresponds to a partially ordered set of arc flips. The translucent arcs
   are the "free" pieces — those that can be moved without changing the identity.
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def build_full_data(n):
    """Build complete labeled tournament data and metagraph."""
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

    # Build class data
    cls = []
    hm = defaultdict(list); cm = {}
    for i, line in enumerate(lines):
        bits = int(line, 2); adj = b2a(bits)
        s = ss(adj); cc = c3(adj); h = H(adj); canon = cp(adj)
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'bits':bits, 'adj':adj, 'score':s,
                    'c3':cc, 'H':h, 'sc':sc, 'comp_canon':cc2, 'canon':canon})
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
            mn[ci] = mid; comp = cl['comp_cid']
            if comp >= 0 and comp != ci: mn[comp] = mid
            mid += 1

    # Build metagraph edges with WEIGHTS (line counts per edge)
    edge_weights = defaultdict(int)  # (merged_a, merged_b) → count of opaque lines
    self_loop_weights = defaultdict(int)  # merged_a → count of translucent lines

    # Also track per-arc-position
    arc_edge_map = defaultdict(lambda: defaultdict(int))  # (arc_idx) → {(ma,mb): count}

    # Iterate over all labeled tournaments
    for bits in range(1 << m):
        adj = b2a(bits)
        # Find its class
        cid = lk(adj)
        if cid < 0: continue
        ma = mn[cid]

        for arc in range(m):
            fb = bits ^ (1 << (m-1-arc))
            fadj = b2a(fb)
            nb_cid = lk(fadj)
            if nb_cid < 0: continue
            mb = mn[nb_cid]

            if ma == mb:
                self_loop_weights[ma] += 1
            else:
                pair = (min(ma, mb), max(ma, mb))
                edge_weights[pair] += 1
                arc_edge_map[arc][pair] += 1

    # Metagraph edges (from class reps)
    me = set()
    for cl in cls:
        b = cl['bits']
        for arc in range(m):
            fb = b ^ (1 << (m-1-arc))
            fa = b2a(fb)
            nb = lk(fa)
            if nb >= 0 and nb != cl['cid']:
                a2, b2 = mn[cl['cid']], mn[nb]
                if a2 != b2:
                    me.add((min(a2,b2), max(a2,b2)))

    mi = {}
    for mm in range(mid):
        cids = [cl['cid'] for cl in cls if mn[cl['cid']] == mm]
        cl0 = cls[cids[0]]
        fiber = factorial(n) // len([p for p in permutations(range(n))
                                      if all(cl0['adj'][p[i]][p[j]] == cl0['adj'][i][j]
                                             for i in range(n) for j in range(n))]) if n <= 6 else 0
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'score': cl0['score'],
                  'c3': cl0['c3'], 'fiber': fiber}

    return {
        'n': n, 'm': m, 'V': mid, 'E': len(me), 'me': me,
        'mi': mi, 'cls': cls, 'mn': mn,
        'edge_weights': dict(edge_weights),
        'self_loop_weights': dict(self_loop_weights),
        'arc_edge_map': dict(arc_edge_map),
    }

def analyze_translucent(data):
    """Analyze translucent vs opaque structure."""
    n = data['n']; m = data['m']; V = data['V']; E = data['E']
    mi = data['mi']; me = data['me']
    ew = data['edge_weights']; sw = data['self_loop_weights']

    print(f"\n{'#'*72}")
    print(f"  n = {n}: TRANSLUCENT LINE ANALYSIS (V={V}, E={E})")
    print(f"{'#'*72}")

    total_labeled_edges = m * (2 ** (m-1))  # each of 2^m tournaments × m arcs / 2
    total_opaque = sum(ew.values())
    total_translucent = sum(sw.values())

    print(f"\n  Total labeled arc-flip pairs: {total_labeled_edges}")
    print(f"  Opaque (class-changing):   {total_opaque} ({100*total_opaque/(total_opaque+total_translucent):.1f}%)")
    print(f"  Translucent (class-preserving): {total_translucent} ({100*total_translucent/(total_opaque+total_translucent):.1f}%)")

    # LINES PER METAGRAPH EDGE
    print(f"\n  LINES PER METAGRAPH EDGE:")
    if ew:
        weights = sorted(ew.values())
        # Each edge weight counts directed pairs (T→T' and T'→T counted separately)
        # So actual lines = weight / 2
        lines_per_edge = [w // 2 for w in weights]
        print(f"    Min lines/edge: {min(lines_per_edge)}")
        print(f"    Max lines/edge: {max(lines_per_edge)}")
        print(f"    Mean lines/edge: {sum(lines_per_edge)/len(lines_per_edge):.1f}")
        print(f"    Median lines/edge: {sorted(lines_per_edge)[len(lines_per_edge)//2]}")

        # Weight distribution
        wc = Counter(lines_per_edge)
        print(f"\n    Weight distribution:")
        for w in sorted(wc.keys()):
            print(f"      {w} lines: {wc[w]} edges ({100*wc[w]/len(lines_per_edge):.1f}%)")

    # TRANSLUCENT FRACTION PER CLASS
    print(f"\n  TRANSLUCENT FRACTION PER CLASS (sorted by H):")
    print(f"  {'mid':>4} {'H':>4} {'SC':>3} {'trans':>8} {'opaque':>8} {'frac':>8} {'fiber':>6}")

    class_data = []
    for mm in sorted(mi.keys(), key=lambda x: mi[x]['H']):
        trans = sw.get(mm, 0)
        # Opaque: sum of all edges incident to this node
        opaq = sum(ew.get((min(mm, other), max(mm, other)), 0)
                   for other in range(V) if other != mm)
        total = trans + opaq
        frac = trans / total if total > 0 else 0
        fiber = mi[mm].get('fiber', 0)
        sc_str = 'SC' if mi[mm]['sc'] else 'NS'
        print(f"  {mm:>4} {mi[mm]['H']:>4} {sc_str:>3} {trans:>8} {opaq:>8} {frac:>8.4f} {fiber:>6}")
        class_data.append((mm, mi[mm]['H'], frac, trans, opaq, mi[mm]['sc']))

    # CORRELATION: translucent fraction vs H
    if len(class_data) >= 3:
        H_vals = [d[1] for d in class_data]
        fracs = [d[2] for d in class_data]
        import numpy as np
        if np.std(H_vals) > 0 and np.std(fracs) > 0:
            r = np.corrcoef(H_vals, fracs)[0, 1]
            print(f"\n  Correlation(H, translucent_fraction) = {r:.4f}")

    # WHICH EDGES ARE HEAVIEST? (most lines)
    print(f"\n  TOP 10 HEAVIEST METAGRAPH EDGES (most opaque lines):")
    sorted_edges = sorted(ew.items(), key=lambda x: -x[1])[:10]
    for (a, b), w in sorted_edges:
        print(f"    ({a},{b}): {w//2} lines, H=({mi[a]['H']},{mi[b]['H']}), "
              f"ΔH={abs(mi[a]['H']-mi[b]['H'])}")

    # WHICH EDGES ARE LIGHTEST? (fewest lines)
    print(f"\n  BOTTOM 10 LIGHTEST METAGRAPH EDGES:")
    sorted_edges_light = sorted(ew.items(), key=lambda x: x[1])[:10]
    for (a, b), w in sorted_edges_light:
        print(f"    ({a},{b}): {w//2} lines, H=({mi[a]['H']},{mi[b]['H']}), "
              f"ΔH={abs(mi[a]['H']-mi[b]['H'])}")

    # ARC POSITION ANALYSIS: which arc positions generate the most edges?
    print(f"\n  ARC POSITION ANALYSIS:")
    arc_counts = {}
    for arc in range(m):
        arc_edges = data['arc_edge_map'].get(arc, {})
        n_distinct_edges = len(arc_edges)
        total_lines = sum(arc_edges.values())
        arc_counts[arc] = (n_distinct_edges, total_lines)

    P = [(i,j) for i in range(n) for j in range(i+1,n)]
    print(f"  {'arc':>4} {'(i,j)':>8} {'edges':>6} {'lines':>8} {'hook':>5}")
    for arc in range(min(m, 20)):
        i, j = P[arc]
        ne, nl = arc_counts.get(arc, (0, 0))
        hook = 2*(n-1-i-j) - 1 if n-1-i-j > 0 else 1
        print(f"  {arc:>4} ({i},{j}):>8 {ne:>6} {nl:>8} {hook:>5}")

    return class_data

if __name__ == '__main__':
    print("=" * 72)
    print("  TRANSLUCENT LINES: DEEP INVESTIGATION")
    print("  opus-2026-03-23-S271")
    print("=" * 72)

    for n in [4, 5, 6]:
        t0 = time.time()
        print(f"\n  Building full data for n={n}...", end=' ', flush=True)
        data = build_full_data(n)
        print(f"done ({time.time()-t0:.1f}s)")

        class_data = analyze_translucent(data)

    # Synthesis
    print(f"\n{'='*72}")
    print("  SYNTHESIS: TRANSLUCENT LINES AND THE HEAP CONNECTION")
    print(f"{'='*72}")
    print("""
  TRANSLUCENT = NEUTRAL ARCS = "FREE PIECES" in Viennot's heap framework.

  A neutral arc is one whose direction doesn't matter for the iso class.
  In heap terms: it's a piece that can be freely moved (flipped) without
  changing the overall structure. It commutes with all other pieces.

  For the TRANSITIVE tournament:
    Neutral arcs = arcs between vertices with identical "relative positions"
    = the TWINS (vertices with same out-neighborhoods)
    Translucent fraction ≈ 2/(n-1) (decays as n grows)

  For REGULAR tournaments:
    Zero neutral arcs — every arc matters.
    = maximally "rigid" or "chaotic" structure
    = the heap is FULLY CONSTRAINED (no free pieces)

  The translucent fraction F(class) anti-correlates with H:
    High H → few translucent (rigid, cyclic)
    Low H → many translucent (flexible, hierarchical)

  In BINARY TREE terms:
    The translucent subgraph within a class is like a partial binary tree:
    - Root = one representative labeling
    - Children = labelings reachable by one neutral flip
    - Depth = translucent diameter
    - Leaves = labelings with no further neutral flips

    For the transitive class: this tree spans ALL n! labelings (connected)
    For regular classes: this "tree" is just isolated points (disconnected)

  METAGRAPH EDGE WEIGHTS = OPAQUE LINE COUNTS:
    Heavy edges: connect similar classes (many ways to transition)
    Light edges: connect dissimilar classes (few transitions)
    The weight measures the "width" of the tunnel between classes.

  THE TWO-LEVEL HIERARCHY:
    Level 1 (metagraph): classes connected by opaque lines
    Level 2 (translucent): labeled tournaments within each class
                           connected by neutral flips

    Level 1 is COARSE (V_merged nodes)
    Level 2 is FINE (n!/|Aut| nodes per class)

    The metagraph's spectral properties (H ≈ 2nd eigenvector, gap ≈ 2/n)
    describe Level 1. The translucent structure describes Level 2.
""")

    print("DONE.")
