#!/usr/bin/env python3
"""
constrain_H_aut_s283.py — Using W and B to separate H and |Aut|
opus-2026-03-24-S283

We know: fiber(C) = H(C)/|Aut(C)| × mult(C)
         fiber = W_row / m (from wiggly weights)

This gives H/|Aut|. To get H and |Aut| SEPARATELY, we need more info.

INSIGHT: The wiggly self-loop W[i][i] counts how many (tiling, tile) pairs
in class i are neutral. The complement B connects tilings to their
cycle-space complements. Together, they give TWO independent equations
per class that can separate H and |Aut|.

Also: the TOTAL CONSTRAINT:
  Σ_{all unmerged C} H(C)/|Aut(C)| = 2^{C(n-1,2)}  [tiling total]
  Σ_{all unmerged C} n!/|Aut(C)| = 2^{C(n,2)}       [labeled total]

Dividing: Σ (n!/|Aut|) / Σ (H/|Aut|) = 2^{n-1}
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def build_and_compute(n):
    """Build everything needed."""
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
    def aut_size(a):
        return sum(1 for p in permutations(range(n))
                   if all(a[p[i]][p[j]] == a[i][j] for i in range(n) for j in range(n)))
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
        aut = aut_size(adj) if n <= 6 else 0
        comp = [[1-adj[a][b] if a!=b else 0 for b in range(n)] for a in range(n)]
        cc2 = cp(comp); sc = (canon == cc2)
        cls.append({'cid':i, 'score':s, 'c3':cc, 'H':h, 'sc':sc,
                    'comp_canon':cc2, 'canon':canon, 'aut':aut})
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
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'aut': cl0['aut']}

    # Tiling computation
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

    total = 2 ** m_tiles
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

    # W and B matrices
    V = mid
    W = np.zeros((V, V), dtype=int)
    B = np.zeros((V, V), dtype=int)
    fiber = Counter()

    for tb in range(total):
        ma = tiling_class[tb]
        if ma < 0: continue
        fiber[ma] += 1
        for ti in range(m_tiles):
            flipped = tb ^ (1 << (m_tiles - 1 - ti))
            mb = tiling_class[flipped]
            if mb >= 0: W[ma][mb] += 1
        comp = ((1 << m_tiles) - 1) ^ tb
        mb_c = tiling_class[comp]
        if mb_c >= 0: B[ma][mb_c] += 1

    return {
        'n': n, 'V': V, 'm_tiles': m_tiles, 'mi': mi,
        'W': W, 'B': B, 'fiber': dict(fiber),
        'nf': factorial(n),
    }

def analyze(data):
    n = data['n']; V = data['V']; m = data['m_tiles']
    mi = data['mi']; W = data['W']; B = data['B']
    fiber = data['fiber']; nf = data['nf']

    print(f"\n{'#'*72}")
    print(f"  n = {n}: CONSTRAINING H AND |Aut| FROM W AND B")
    print(f"{'#'*72}")

    # The sum constraints
    sum_fiber = sum(fiber.values())
    sum_H_over_aut = 0
    sum_nf_over_aut = 0
    for mm in mi:
        h = mi[mm]['H']; aut = mi[mm]['aut']; sc = mi[mm]['sc']
        if aut > 0:
            sum_H_over_aut += h / aut * (1 if sc else 2)
            sum_nf_over_aut += nf / aut * (1 if sc else 2)

    print(f"\n  SUM CONSTRAINTS:")
    print(f"  Σ fiber = {sum_fiber} = 2^{m} = {2**m}")
    if n <= 6:
        print(f"  Σ H/|Aut| × mult = {sum_H_over_aut:.0f} = 2^{m} = {2**m}")
        print(f"  Σ n!/|Aut| × mult = {sum_nf_over_aut:.0f} = 2^{comb(n,2)} = {2**comb(n,2)}")

    # For each class: fiber, W_self, B_self, W_degree, B_off_diag
    print(f"\n  PER-CLASS DATA:")
    print(f"  {'mid':>4} {'H':>4} {'|Aut|':>5} {'SC':>3} {'fiber':>6} "
          f"{'W_self':>7} {'B_self':>7} {'B_off':>6} {'W_deg':>6} "
          f"{'neut/m':>7}")

    for mm in sorted(mi.keys()):
        h = mi[mm]['H']; aut = mi[mm]['aut']; sc = mi[mm]['sc']
        f = fiber.get(mm, 0)
        ws = W[mm][mm]; bs = B[mm][mm]
        b_off = sum(B[mm]) - bs
        w_deg = sum(1 for j in range(V) if j != mm and W[mm][j] > 0)
        neut_per_tiling = ws / f if f > 0 else 0

        print(f"  {mm:>4} {h:>4} {aut:>5} {'SC' if sc else 'NS':>3} {f:>6} "
              f"{ws:>7} {bs:>7} {b_off:>6} {w_deg:>6} "
              f"{neut_per_tiling:>7.2f}")

    # THE KEY: can W_self and B_self separate H and |Aut|?
    # fiber = H/|Aut| × mult (known)
    # W_self = fiber × (avg neutral arcs)
    # B_self = #{tilings whose cycle-complement is in same class}

    # Check: does W_self correlate with |Aut|?
    if n <= 6:
        aut_vals = [mi[mm]['aut'] for mm in sorted(mi.keys())]
        ws_vals = [W[mm][mm] for mm in sorted(mi.keys())]
        fiber_vals = [fiber.get(mm, 0) for mm in sorted(mi.keys())]

        if np.std(aut_vals) > 0 and np.std(ws_vals) > 0:
            r_aut_ws = np.corrcoef(aut_vals, ws_vals)[0, 1]
            print(f"\n  Corr(|Aut|, W_self) = {r_aut_ws:.4f}")

        # Does W_self / fiber correlate with |Aut|?
        neut_frac = [W[mm][mm] / fiber.get(mm, 1) for mm in sorted(mi.keys())]
        if np.std(neut_frac) > 0:
            r_aut_nf = np.corrcoef(aut_vals, neut_frac)[0, 1]
            print(f"  Corr(|Aut|, neutral_frac) = {r_aut_nf:.4f}")

        # Does B_self / fiber tell us about SC structure?
        bs_frac = [B[mm][mm] / fiber.get(mm, 1) for mm in sorted(mi.keys())]
        if np.std(bs_frac) > 0:
            r_aut_bs = np.corrcoef(aut_vals, bs_frac)[0, 1]
            print(f"  Corr(|Aut|, B_self_frac) = {r_aut_bs:.4f}")

    # THE DOUBLE EQUATION SYSTEM:
    # Equation 1: fiber = H/|Aut| × mult   → known from W
    # Equation 2: n!/|Aut| = labeled class size  → known from Burnside
    # So: H = fiber × |Aut| / mult = fiber × (n! / class_size) / mult

    # But class_size = n!/|Aut| is NOT directly available from the tiling model...
    # Unless we can get it from W or B.

    # IDEA: The number of DISTINCT B-neighbors of class i.
    # B_off_diag(i) = fiber(i) - B_self(i) = #{tilings whose complement is elsewhere}
    # The number of DISTINCT classes j with B[i][j] > 0 tells us about
    # the "complement spread" of class i.

    print(f"\n  COMPLEMENT SPREAD:")
    for mm in sorted(mi.keys()):
        b_neighbors = sum(1 for j in range(V) if j != mm and B[mm][j] > 0)
        f = fiber.get(mm, 0)
        print(f"    Class {mm} (H={mi[mm]['H']}): {b_neighbors} complement neighbors, "
              f"B_self={B[mm][mm]}/{f}")

    # THE TILING IDENTITY:
    # Σ fiber(C) = 2^m
    # Σ fiber(C)^2 / total = <fiber> = E[fiber over random tiling]
    # This gives the "concentration" of tilings in classes.
    fiber_sq_sum = sum(f**2 for f in fiber.values())
    concentration = fiber_sq_sum / sum_fiber
    print(f"\n  CONCENTRATION: Σ fiber² / Σ fiber = {concentration:.2f}")
    print(f"  If uniform: concentration = {sum_fiber / V:.2f}")
    print(f"  Ratio: {concentration / (sum_fiber / V):.4f}")

    # The Cauchy-Schwarz bound: Σ H/|Aut| ≤ sqrt(Σ (H/|Aut|)²) × sqrt(V)
    # Since Σ H/|Aut| = 2^m, this gives a bound on the variance of H/|Aut|.

if __name__ == '__main__':
    print("=" * 72)
    print("  CONSTRAINING H AND |Aut| FROM DUAL METRICS")
    print("  opus-2026-03-24-S283")
    print("=" * 72)

    for n in [4, 5, 6]:
        t0 = time.time()
        data = build_and_compute(n)
        analyze(data)
        print(f"\n  Time: {time.time()-t0:.1f}s")

    # The grand identity
    print(f"\n{'='*72}")
    print("  THE GRAND IDENTITY")
    print(f"{'='*72}")
    print("""
  THREE COUNTING IDENTITIES on the same set of V_n nodes:

  1. TILING SUM:    Σ H(C)/|Aut(C)| × mult(C) = 2^{C(n-1,2)}
  2. LABELED SUM:   Σ n!/|Aut(C)| × mult(C) = 2^{C(n,2)}
  3. SZELE AVERAGE:  Σ H(C) × n!/|Aut(C)| × mult(C) = n! × 2^{C(n-1,2)}

  From 1 and 2:
     Σ n!/|Aut| / Σ H/|Aut| = 2^{n-1}
     → the weighted average of n!/H over classes is 2^{n-1}

  The WIGGLY MATRIX W determines identity 1 (via row sums).
  The COMPLEMENT MATRIX B provides a SECOND independent view.
  The BLUE/BLACK classification (SC vs NS) determines mult.

  Together: W + B + SC/NS → three equations → constrain H, |Aut|, mult.

  WHAT'S STILL NEEDED for full separation:
  |Aut| cannot be read from W or B alone (it requires knowing the
  full S_n action, not just the tiling model with fixed base path).
  The tiling model has Stab(P_0) = {id}, so it cannot see |Aut| > 1.

  BUT: H/|Aut| IS determined, and H can be recovered if we know
  the EIGENVALUES of W (since H ≈ 2nd eigenvector of the metagraph).
""")

    print("DONE.")
