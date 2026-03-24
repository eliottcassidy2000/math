#!/usr/bin/env python3
"""
wiggly_shape_s293.py — The wiggly space as a standalone object
opus-2026-03-24-S293

THE OBJECT: The hypercube Q_m where m = C(n-1,2) (triangular number).
  - 2^m tilings (vertices)
  - Each tiling has exactly m neighbors (one per wiggly class)
  - Quotiented by S_n → the wiggly metagraph

WITHOUT blue/black: this is just "binary patterns under permutation."
Each tiling is a binary word of length m. S_n acts by permuting bits.
Orbits = iso classes. The quotient graph = orbits connected by single-bit flips.

WHAT VARIES across tilings (within the same class):
  - #distinct classes among neighbors
  - Average degree of neighboring classes
  - Local clustering coefficient

WHAT'S CONSTANT:
  - Total neighbors = m (always)
  - One neighbor per wiggly class (always)
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_data(n):
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

    m_tiles = comb(n-1, 2)
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2: tiles.append((x, y))
    tile_to_P = {}; base_P = set()
    for k, (i, j) in enumerate(P):
        if j == i + 1: base_P.add(k)
        else:
            for t_idx, (tx, ty) in enumerate(tiles):
                if (ty-1, tx-1) == (i, j):
                    tile_to_P[t_idx] = k; break

    total = 2 ** m_tiles; V = mid
    tiling_class = [0] * total
    for tb in range(total):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k not in base_P:
                for ti in range(m_tiles):
                    if tile_to_P.get(ti) == k:
                        bit = (tb >> (m_tiles - 1 - ti)) & 1
                        if bit: full_bits |= (1 << (m_total - 1 - k))
                        break
        adj = b2a(full_bits); cid = lk(adj)
        tiling_class[tb] = mn[cid] if cid >= 0 else -1

    return V, mi, tiling_class, m_tiles, total

def analyze_shape(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: THE WIGGLY SPACE Q_{comb(n-1,2)}")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, tc, m, total = build_data(n)
    print(f"  Built in {time.time()-t0:.1f}s")
    print(f"  m = C({n-1},2) = {m} (triangular number)")
    print(f"  2^m = {total} tilings, V = {V} merged classes")
    print(f"  Each tiling has exactly {m} wiggly neighbors")

    # For each tiling: count distinct classes among its m neighbors
    # And count self-loops (neighbor in same class)
    class_diversity = []  # per tiling: #distinct classes among neighbors
    self_loop_count = []  # per tiling: #neighbors in same class
    neighbor_class_dist = defaultdict(list)  # class → list of diversities

    fiber = Counter()
    for tb in range(total):
        if tc[tb] < 0: continue
        fiber[tc[tb]] += 1

    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue

        neighbor_classes = set()
        self_loops = 0
        for ti in range(m):
            flipped = tb ^ (1 << (m - 1 - ti))
            mb = tc[flipped]
            if mb >= 0:
                neighbor_classes.add(mb)
                if mb == ma: self_loops += 1

        diversity = len(neighbor_classes)
        class_diversity.append(diversity)
        self_loop_count.append(self_loops)
        neighbor_class_dist[ma].append(diversity)

    print(f"\n  PER-TILING: #DISTINCT CLASSES AMONG {m} NEIGHBORS")
    print(f"  Min: {min(class_diversity)}, Max: {max(class_diversity)}, "
          f"Mean: {np.mean(class_diversity):.2f}, Std: {np.std(class_diversity):.2f}")

    dist = Counter(class_diversity)
    print(f"  Distribution:")
    for k in sorted(dist.keys()):
        print(f"    {k} distinct classes: {dist[k]} tilings ({100*dist[k]/total:.1f}%)")

    # Per-class average diversity
    print(f"\n  PER-CLASS AVERAGE DIVERSITY:")
    print(f"  {'mid':>4} {'H':>4} {'SC':>3} {'fiber':>6} {'mean_div':>9} {'std_div':>8} "
          f"{'min_div':>8} {'max_div':>8}")

    class_stats = []
    for mm in sorted(mi.keys()):
        divs = neighbor_class_dist.get(mm, [])
        if not divs: continue
        mean_d = np.mean(divs)
        std_d = np.std(divs)
        class_stats.append((mm, mi[mm]['H'], mi[mm]['sc'], fiber[mm],
                           mean_d, std_d, min(divs), max(divs)))
        if n <= 6 or (mm < 5 or mm > V - 5):
            print(f"  {mm:>4} {mi[mm]['H']:>4} {'SC' if mi[mm]['sc'] else 'NS':>3} "
                  f"{fiber[mm]:>6} {mean_d:>9.2f} {std_d:>8.2f} "
                  f"{min(divs):>8} {max(divs):>8}")

    # SELF-LOOP ANALYSIS
    print(f"\n  SELF-LOOPS (neighbors in same class):")
    sl_dist = Counter(self_loop_count)
    print(f"  Min: {min(self_loop_count)}, Max: {max(self_loop_count)}, "
          f"Mean: {np.mean(self_loop_count):.2f}")
    for k in sorted(sl_dist.keys()):
        if sl_dist[k] >= total * 0.01 or k <= 3:
            print(f"    {k} self-loops: {sl_dist[k]} tilings ({100*sl_dist[k]/total:.1f}%)")

    # Correlation: diversity vs H
    if class_stats:
        H_vals = [s[1] for s in class_stats]
        mean_divs = [s[4] for s in class_stats]
        if np.std(H_vals) > 0 and np.std(mean_divs) > 0:
            r = np.corrcoef(H_vals, mean_divs)[0, 1]
            print(f"\n  Corr(H, mean_diversity) = {r:.4f}")

    # THE SHAPE: Hamming weight distribution per class
    print(f"\n  HAMMING WEIGHT (number of 1-bits = flipped tiles) PER CLASS:")
    hw_per_class = defaultdict(list)
    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue
        hw = bin(tb).count('1')
        hw_per_class[ma].append(hw)

    print(f"  {'mid':>4} {'H':>4} {'mean_hw':>8} {'std_hw':>7} {'min':>4} {'max':>4}")
    for mm in sorted(mi.keys()):
        hws = hw_per_class.get(mm, [])
        if not hws: continue
        if n <= 6 or (mm < 5 or mm > V - 5):
            print(f"  {mm:>4} {mi[mm]['H']:>4} {np.mean(hws):>8.2f} {np.std(hws):>7.2f} "
                  f"{min(hws):>4} {max(hws):>4}")

    # Correlation: H vs mean Hamming weight
    if class_stats:
        mean_hws = [np.mean(hw_per_class.get(s[0], [0])) for s in class_stats]
        if np.std(mean_hws) > 0:
            r_hw = np.corrcoef(H_vals, mean_hws)[0, 1]
            print(f"\n  Corr(H, mean_hamming_weight) = {r_hw:.4f}")

    # THE TRANSITIVE = all-zeros tiling (Hamming weight 0)
    all_zeros_class = tc[0]
    # THE ALL-ONES tiling (Hamming weight m) = "anti-transitive"
    all_ones_class = tc[(1 << m) - 1]
    print(f"\n  SPECIAL TILINGS:")
    print(f"  All-zeros (transitive): class {all_zeros_class}, "
          f"H={mi[all_zeros_class]['H']}")
    print(f"  All-ones (anti-trans):  class {all_ones_class}, "
          f"H={mi[all_ones_class]['H']}")

    # Distance from transitive (Hamming distance to all-zeros)
    # = Hamming weight of the tiling
    # So mean Hamming weight = mean distance from transitive

    return V, mi

if __name__ == '__main__':
    print("=" * 72)
    print("  THE WIGGLY SPACE: Q_m AS A STANDALONE OBJECT")
    print("  opus-2026-03-24-S293")
    print("=" * 72)

    for n in [4, 5, 6]:
        analyze_shape(n)

    print(f"\n{'='*72}")
    print("  WHAT IS THE WIGGLY SPACE?")
    print(f"{'='*72}")
    print(f"""
  WITHOUT blue/black lines, the wiggly space is:

  THE HYPERCUBE Q_m with m = C(n-1, 2) = TRIANGULAR NUMBER
    m = 1, 3, 6, 10, 15, 21, 28, 36 for n=3..10

  Each vertex = a binary word of length m (a tiling of δ_{{n-2}})
  Each edge = a single-bit flip (one wiggly line)
  Each vertex has EXACTLY m neighbors, one per wiggly class

  The GROUP S_n acts on Q_m by permuting bits (relabeling vertices).
  Orbits = tournament iso classes.
  The quotient Q_m / S_n = the wiggly metagraph.

  THE UNDERLYING SHAPE:
  - Q_m is vertex-transitive (uniform before quotienting)
  - After quotienting: diversity varies by class
  - Low-H classes: neighbors spread across FEW classes (concentrated)
  - High-H classes: neighbors spread across MANY classes (diverse)
  - Hamming weight ↔ distance from transitive ↔ "how flipped"
  - H correlates with Hamming weight (more flips → more H, roughly)

  THIS IS A CODING THEORY OBJECT:
  - Binary code of length m = C(n-1,2) = triangular number
  - Group action = S_n on the coordinate set
  - The "code" = all 2^m words (no redundancy)
  - The "quotient code" = the iso classes
  - Hamming distance = number of different tiles
  - The metagraph = the quotient Hamming graph
""")

    print("DONE.")
