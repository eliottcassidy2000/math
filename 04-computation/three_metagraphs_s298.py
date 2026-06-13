#!/usr/bin/env python3
"""
three_metagraphs_s298.py — The three metagraphs in tandem
opus-2026-03-24-S298

THREE METAGRAPHS on the same node set (merged iso classes):

1. WIGGLY (W): edges from single-tile flips (order 1)
2. WAGGLY-2 (W2): edges from 2-tile flips (order 2)
3. BLUE/BLACK (B): edges from tiling complement (order m = flip ALL tiles)

Each gives a weighted graph on V_merged nodes.
Study them together: overlaps, differences, spectra, correlations.
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_everything(n):
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
        cl0 = cls[cids[0]]; mi[mm] = {'H': cl0['H'], 'sc': cl0['sc']}
    m = comb(n-1, 2)
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
    total = 2 ** m; V = mid
    tc = [0] * total
    for tb in range(total):
        full_bits = 0
        for k, (i, j) in enumerate(P):
            if k not in base_P:
                for ti in range(m):
                    if tile_to_P.get(ti) == k:
                        bit = (tb >> (m - 1 - ti)) & 1
                        if bit: full_bits |= (1 << (m_total - 1 - k))
                        break
        adj = b2a(full_bits); cid = lk(adj)
        tc[tb] = mn[cid] if cid >= 0 else -1

    # Build all three weight matrices
    W1 = np.zeros((V, V), dtype=int)  # wiggly (order 1)
    W2 = np.zeros((V, V), dtype=int)  # waggly-2 (order 2)
    B = np.zeros((V, V), dtype=int)   # complement (order m)
    fiber = Counter()

    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue
        fiber[ma] += 1

        # Wiggly (order 1)
        for ti in range(m):
            mb = tc[tb ^ (1 << (m-1-ti))]
            if mb >= 0: W1[ma][mb] += 1

        # Waggly-2 (order 2)
        for ti in range(m):
            for tj in range(ti+1, m):
                mb = tc[tb ^ (1 << (m-1-ti)) ^ (1 << (m-1-tj))]
                if mb >= 0: W2[ma][mb] += 1

        # Complement (order m)
        mb = tc[((1 << m) - 1) ^ tb]
        if mb >= 0: B[ma][mb] += 1

    return V, mi, W1, W2, B, dict(fiber), m, total

def analyze_tandem(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: THREE METAGRAPHS IN TANDEM")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, W1, W2, B, fiber, m, total = build_everything(n)
    print(f"  Built in {time.time()-t0:.1f}s. V={V}, m={m}")

    # Edge sets
    E1 = {(i,j) for i in range(V) for j in range(i+1,V) if W1[i][j] > 0}
    E2 = {(i,j) for i in range(V) for j in range(i+1,V) if W2[i][j] > 0}
    EB = {(i,j) for i in range(V) for j in range(i+1,V) if B[i][j] > 0}
    total_pairs = V*(V-1)//2

    print(f"\n  EDGE COUNTS:")
    print(f"    Wiggly (E1):      {len(E1):>5} ({100*len(E1)/total_pairs:.1f}%)")
    print(f"    Waggly-2 (E2):    {len(E2):>5} ({100*len(E2)/total_pairs:.1f}%)")
    print(f"    Complement (EB):  {len(EB):>5} ({100*len(EB)/total_pairs:.1f}%)")

    # Venn diagram
    only1 = E1 - E2 - EB
    only2 = E2 - E1 - EB
    onlyB = EB - E1 - E2
    e12 = (E1 & E2) - EB
    e1B = (E1 & EB) - E2
    e2B = (E2 & EB) - E1
    e12B = E1 & E2 & EB
    none = total_pairs - len(E1 | E2 | EB)

    print(f"\n  VENN DIAGRAM:")
    print(f"    Only wiggly:        {len(only1):>5}")
    print(f"    Only waggly-2:      {len(only2):>5}")
    print(f"    Only complement:    {len(onlyB):>5}")
    print(f"    Wiggly ∩ waggly-2:  {len(e12):>5}")
    print(f"    Wiggly ∩ complement:{len(e1B):>5}")
    print(f"    Waggly-2 ∩ comp:    {len(e2B):>5}")
    print(f"    ALL THREE:          {len(e12B):>5}")
    print(f"    NONE (gap):         {none:>5}")

    # Symmetry checks
    print(f"\n  SYMMETRY:")
    print(f"    W1 symmetric: {np.allclose(W1, W1.T)}")
    print(f"    W2 symmetric: {np.allclose(W2, W2.T)}")
    print(f"    B symmetric:  {np.allclose(B, B.T)}")

    # Row sums
    print(f"\n  ROW SUMS:")
    print(f"    W1_row = fiber × {m} (check: {np.allclose([np.sum(W1[i]) for i in range(V)], [fiber.get(i,0)*m for i in range(V)])})")
    print(f"    W2_row = fiber × C({m},2) = fiber × {comb(m,2)} (check: {np.allclose([np.sum(W2[i]) for i in range(V)], [fiber.get(i,0)*comb(m,2) for i in range(V)])})")
    print(f"    B_row = fiber × 1 = fiber (check: {np.allclose([np.sum(B[i]) for i in range(V)], [fiber.get(i,0) for i in range(V)])})")

    # RATIO W2/W1 per edge
    both_12 = E1 & E2
    if both_12:
        ratios = []
        for (i,j) in both_12:
            r = W2[i][j] / W1[i][j] if W1[i][j] > 0 else 0
            ratios.append(r)
        print(f"\n  W2/W1 RATIO (for edges in both):")
        print(f"    Min: {min(ratios):.2f}, Max: {max(ratios):.2f}, Mean: {np.mean(ratios):.2f}")

    # Spectra comparison
    print(f"\n  SPECTRA (top 5 eigenvalues):")
    for name, M in [("W1", W1), ("W2", W2), ("B", B)]:
        evals = sorted(np.linalg.eigvalsh(M.astype(float)), reverse=True)
        print(f"    {name}: [{', '.join(f'{e:.1f}' for e in evals[:5])}]")

    # H correlation with each graph's 2nd eigenvector
    H_vals = np.array([mi[i]['H'] for i in range(V)], dtype=float)
    H_centered = H_vals - np.mean(H_vals)
    H_norm = np.linalg.norm(H_centered)

    if H_norm > 0:
        H_unit = H_centered / H_norm
        print(f"\n  H PROJECTION ONTO 2ND EIGENVECTOR:")
        for name, M in [("W1 (wiggly)", W1), ("W2 (waggly-2)", W2), ("B (complement)", B)]:
            evals, evecs = np.linalg.eigh(M.astype(float))
            idx = np.argsort(-evals)
            projs = np.abs(evecs[:, idx].T @ H_unit) ** 2
            best_i = np.argmax(projs[1:]) + 1  # skip trivial
            print(f"    {name}: v_1 captures {100*projs[1]:.1f}% of H, "
                  f"best non-trivial = v_{best_i} at {100*projs[best_i]:.1f}%")

    # Per-edge ΔH comparison
    print(f"\n  MEAN |ΔH| BY EDGE TYPE:")
    for name, edges in [("Only wiggly", only1), ("Only waggly-2", only2),
                         ("Only complement", onlyB), ("W1∩W2", e12),
                         ("W1∩B", e1B), ("W2∩B", e2B), ("All three", e12B)]:
        if edges:
            dh = [abs(mi[a]['H'] - mi[b]['H']) for (a,b) in edges]
            print(f"    {name:>20}: {np.mean(dh):>6.1f} ({len(edges)} edges)")

    # The COMBINED weight: W_total = W1 + W2 + B
    # Is there a formula relating them?
    print(f"\n  WEIGHT RELATIONSHIPS:")
    # For each node: W1_diag + W2_diag + B_diag = ?
    for i in range(min(V, 10)):
        w1d = W1[i][i]; w2d = W2[i][i]; bd = B[i][i]
        f = fiber.get(i, 0)
        total_self = w1d + w2d + bd
        total_possible = f * (m + comb(m,2) + 1)
        print(f"    Node {i} (H={mi[i]['H']}): W1_self={w1d}, W2_self={w2d}, B_self={bd}, "
              f"fiber={f}, total_self/fiber={total_self/f:.1f}" if f > 0 else "")

    return V, E1, E2, EB

if __name__ == '__main__':
    print("=" * 72)
    print("  THREE METAGRAPHS IN TANDEM")
    print("  opus-2026-03-24-S298")
    print("=" * 72)

    for n in [4, 5, 6]:
        analyze_tandem(n)

    print("\nDONE.")
