#!/usr/bin/env python3
"""
wiggly_even_graph_s278.py — Wiggly tilings ARE even graphs
opus-2026-03-24-S278

KEY INSIGHT: Each tiling = binary filling of δ_{n-2} = cycle-space element.
The tile values (0/1 per tile) encode a degree-even graph on n vertices.

  tiling ↦ tournament (by adding base path P_0)
  tiling ↦ even graph (the tile values ARE the cycle-space vector)

So the 2^m tilings = 2^m labeled even graphs (m = C(n-1,2)).

QUESTIONS:
1. For each tiling, what is its tournament iso class AND its even graph iso class?
2. How do even graph classes partition the tournament classes?
3. Does each even graph class correspond to exactly one tournament class?
   (If so: this IS the Royle bijection through the tiling model!)
4. How does the wiggly structure (tile flips) interact with even graph isomorphism?
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def build_data(n):
    """Full tournament + even graph data from tiling model."""
    m = comb(n-1, 2)  # number of tiles
    m_total = comb(n, 2)  # total arcs

    # Tiles: (a,b) with a-b ≥ 2, a,b ∈ {1,...,n}
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))

    # All pairs for gentourng (0-indexed)
    P = [(i,j) for i in range(n) for j in range(i+1, n)]

    # Build tile-to-arc mapping
    # Tile (x,y) 1-indexed → arc between vertices x-1 and y-1 (0-indexed)
    # In P: pair (y-1, x-1) since y-1 < x-1
    tile_to_P = {}
    for t_idx, (x, y) in enumerate(tiles):
        pair = (y-1, x-1)
        p_idx = P.index(pair)
        tile_to_P[t_idx] = p_idx

    # Base path arcs: (k, k-1) 1-indexed → (k-2, k-1) 0-indexed
    base_P = set()
    for k in range(2, n+1):
        pair = (k-2, k-1)
        p_idx = P.index(pair)
        base_P.add(p_idx)

    # Get iso class reps from nauty
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    lines = [l.strip() for l in (r.stdout or '').split('\n')
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

    # Canonical form for undirected graphs (even graphs)
    def cp_undir(a):
        degs = [sum(a[i][j] for j in range(n)) for i in range(n)]
        sg = defaultdict(list)
        for v in range(n): sg[degs[v]].append(v)
        gs = [sg[d] for d in sorted(sg.keys())]
        best = None; cnt = 0
        def gp(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for rest in gp(gs2[1:]): yield list(p)+rest
        for p in gp(gs):
            f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(i+1, n))
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
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}

    return {
        'n': n, 'm': m, 'm_total': m_total, 'tiles': tiles,
        'tile_to_P': tile_to_P, 'base_P': base_P, 'P': P,
        'cls': cls, 'lk': lk, 'mn': mn, 'V_merged': mid, 'mi': mi,
        'b2a': b2a, 'cp_undir': cp_undir,
    }

def analyze(data):
    """Map each tiling to both its tournament class and even graph class."""
    n = data['n']; m = data['m']; tiles = data['tiles']
    P = data['P']; lk = data['lk']; mn = data['mn']; mi = data['mi']
    b2a = data['b2a']; cp_undir = data['cp_undir']
    tile_to_P = data['tile_to_P']; base_P = data['base_P']
    m_total = data['m_total']

    print(f"\n{'#'*72}")
    print(f"  n = {n}: TILING → TOURNAMENT CLASS × EVEN GRAPH CLASS")
    print(f"{'#'*72}")

    total_tilings = 2 ** m
    print(f"  Tilings: {total_tilings} = 2^{m}")
    print(f"  Tournament classes (merged): {data['V_merged']}")

    tourn_classes = []  # tournament merged class per tiling
    even_classes = []   # even graph class per tiling
    even_canon_map = {}  # canonical even graph → class id
    even_cid = 0

    for tiling_bits in range(total_tilings):
        # Build tournament from tiling
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k in base_P:
                # Base path: j→i (j = i+1, so j beats i → bit = 0)
                pass
            else:
                # Find tile index
                t_idx = None
                for ti in range(m):
                    if tile_to_P[ti] == k:
                        t_idx = ti
                        break
                if t_idx is not None:
                    bit = (tiling_bits >> (m - 1 - t_idx)) & 1
                    if bit:
                        full_bits |= (1 << (m_total - 1 - k))

        adj = b2a(full_bits)
        cid = lk(adj)
        if cid < 0:
            tourn_classes.append(-1)
            even_classes.append(-1)
            continue
        ma = mn[cid]
        tourn_classes.append(ma)

        # Build even graph from tile values
        # The even graph has edges at tile positions where bit = 1
        even_adj = [[0]*n for _ in range(n)]
        for ti in range(m):
            bit = (tiling_bits >> (m - 1 - ti)) & 1
            if bit:
                x, y = tiles[ti]
                even_adj[x-1][y-1] = 1
                even_adj[y-1][x-1] = 1

        # Canonical even graph
        ec = cp_undir(even_adj)
        if ec not in even_canon_map:
            even_canon_map[ec] = even_cid
            even_cid += 1
        even_classes.append(even_canon_map[ec])

    n_even_classes = len(even_canon_map)
    print(f"  Even graph classes: {n_even_classes}")

    # Cross-tabulation: tournament class × even graph class
    cross = defaultdict(int)
    for i in range(total_tilings):
        if tourn_classes[i] >= 0 and even_classes[i] >= 0:
            cross[(tourn_classes[i], even_classes[i])] += 1

    # How many tournament classes does each even graph class touch?
    even_to_tourn = defaultdict(set)
    tourn_to_even = defaultdict(set)
    for (tc, ec), cnt in cross.items():
        even_to_tourn[ec].add(tc)
        tourn_to_even[tc].add(ec)

    print(f"\n  CROSS-TABULATION:")
    print(f"  Tournament classes per even graph class:")
    for ec in sorted(even_to_tourn.keys()):
        tcs = sorted(even_to_tourn[ec])
        h_vals = [mi[tc]['H'] for tc in tcs]
        n_tilings = sum(cross[(tc, ec)] for tc in tcs)
        print(f"    Even class {ec}: {len(tcs)} tourn classes, "
              f"H={h_vals}, {n_tilings} tilings")

    print(f"\n  Even graph classes per tournament class:")
    for tc in sorted(tourn_to_even.keys()):
        ecs = sorted(tourn_to_even[tc])
        n_tilings = sum(cross[(tc, ec)] for ec in ecs)
        print(f"    Tourn class {tc} (H={mi[tc]['H']}): "
              f"{len(ecs)} even classes, {n_tilings} tilings")

    # Is the map injective? surjective?
    is_tourn_determines_even = all(len(ecs) == 1 for ecs in tourn_to_even.values())
    is_even_determines_tourn = all(len(tcs) == 1 for tcs in even_to_tourn.values())

    print(f"\n  INJECTIVITY CHECK:")
    print(f"  Tournament class determines even class: {is_tourn_determines_even}")
    print(f"  Even class determines tournament class: {is_even_determines_tourn}")

    if is_tourn_determines_even and is_even_determines_tourn:
        print(f"  → BIJECTION! Tournament classes = Even graph classes through tilings!")
    elif not is_tourn_determines_even and not is_even_determines_tourn:
        print(f"  → MANY-TO-MANY: both directions have multiplicity")
    elif is_tourn_determines_even:
        print(f"  → Tournament classes REFINE even graph classes (finer partition)")
    else:
        print(f"  → Even graph classes REFINE tournament classes (finer partition)")

    # Count: #tilings per tournament class
    tilings_per_tourn = Counter(tc for tc in tourn_classes if tc >= 0)
    print(f"\n  TILINGS PER TOURNAMENT CLASS:")
    for tc in sorted(tilings_per_tourn.keys()):
        print(f"    Class {tc} (H={mi[tc]['H']}): {tilings_per_tourn[tc]} tilings")

    # Is #tilings = H(T)? (This would be the key connection!)
    print(f"\n  TILINGS vs H:")
    print(f"  {'Class':>6} {'H':>5} {'#tilings':>9} {'H == #til?':>11}")
    all_match = True
    for tc in sorted(tilings_per_tourn.keys()):
        h = mi[tc]['H']
        nt = tilings_per_tourn[tc]
        match = "✓" if h == nt else "✗"
        if h != nt: all_match = False
        print(f"  {tc:>6} {h:>5} {nt:>9} {match:>11}")

    if all_match:
        print(f"\n  *** H(T) = #tilings for every class! ***")
        print(f"  This means: the number of tilings of class C = H(C).")
        print(f"  Each Hamiltonian path gives a different base-path embedding → tiling.")

    return {
        'n': n, 'n_even': n_even_classes, 'V_merged': data['V_merged'],
        'tilings_match_H': all_match,
    }

if __name__ == '__main__':
    print("=" * 72)
    print("  WIGGLY TILINGS = EVEN GRAPHS?")
    print("  opus-2026-03-24-S278")
    print("=" * 72)

    for n in [3, 4, 5, 6]:
        t0 = time.time()
        data = build_data(n)
        r = analyze(data)
        print(f"\n  Time: {time.time()-t0:.1f}s")

    print("\nDONE.")
