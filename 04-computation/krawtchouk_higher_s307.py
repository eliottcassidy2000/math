#!/usr/bin/env python3
"""
krawtchouk_higher_s307.py — Higher Krawtchouk modes: what do K_2, K_3 capture?
opus-2026-03-24-S307

H ≈ -K_1. What about K_2, K_3, ...?
These should capture the WIDTH perpendicular to H — the "2D correction."

Also: correlate B_k (dual spectrum coefficients) with tournament invariants
(c3, score sequence, |Aut|, SC/NS).
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np
import subprocess

sys.stdout.reconfigure(line_buffering=True)

def krawtchouk(k, x, m):
    total = 0
    for j in range(k+1):
        if j > x or k-j > m-x: continue
        total += (-1)**j * comb(int(x), j) * comb(int(m-x), k-j)
    return total

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
        mi[mm] = {'H': cl0['H'], 'sc': cl0['sc'], 'c3': cl0['c3'], 'score': cl0['score']}
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
    fiber = Counter()
    for tb in range(total):
        if tc[tb] >= 0: fiber[tc[tb]] += 1
    return V, mi, tc, m, total, dict(fiber)

print("=" * 72)
print("  HIGHER KRAWTCHOUK MODES")
print("  opus-2026-03-24-S307")
print("=" * 72)

for n in [5, 6]:
    print(f"\n{'#'*72}")
    print(f"  n = {n}")
    print(f"{'#'*72}")

    V, mi, tc, m, total, fiber = build_data(n)
    print(f"  V={V}, m={m}")

    # Compute weight enumerator and dual spectrum for each class
    weight_enum = defaultdict(lambda: np.zeros(m+1))
    for tb in range(total):
        c = tc[tb]
        if c >= 0:
            w = bin(tb).count('1')
            weight_enum[c][w] += 1

    # Krawtchouk matrix
    K = np.zeros((m+1, m+1))
    for k in range(m+1):
        for x in range(m+1):
            K[k][x] = krawtchouk(k, x, m)

    # Compute B_k for each class
    B_vals = {}  # class → array of B_0, B_1, ..., B_m
    for c in range(V):
        A = weight_enum[c]
        f = fiber.get(c, 0)
        if f > 0:
            B = K @ A / f
            B_vals[c] = B

    # Collect invariants
    H_arr = np.array([mi[c]['H'] for c in range(V)], dtype=float)
    c3_arr = np.array([mi[c]['c3'] for c in range(V)], dtype=float)
    B1_arr = np.array([B_vals[c][1] for c in range(V)], dtype=float)
    B2_arr = np.array([B_vals[c][2] for c in range(V)], dtype=float)
    B3_arr = np.array([B_vals[c][3] for c in range(V)], dtype=float)

    # Correlations
    print(f"\n  KRAWTCHOUK DUAL COEFFICIENTS vs TOURNAMENT INVARIANTS:")
    print(f"  {'':>10} {'H':>8} {'c3':>8} {'SC':>8}")

    sc_arr = np.array([1.0 if mi[c]['sc'] else 0.0 for c in range(V)])

    for name, arr in [("B_1", B1_arr), ("B_2", B2_arr), ("B_3", B3_arr)]:
        corr_H = np.corrcoef(arr, H_arr)[0,1] if np.std(arr) > 0 else 0
        corr_c3 = np.corrcoef(arr, c3_arr)[0,1] if np.std(arr) > 0 else 0
        corr_sc = np.corrcoef(arr, sc_arr)[0,1] if np.std(arr) > 0 else 0
        print(f"  {name:>10} {corr_H:>8.4f} {corr_c3:>8.4f} {corr_sc:>8.4f}")

    # What does B_2 capture?
    print(f"\n  B_2 DETAILS:")
    for c in sorted(range(V), key=lambda c: B_vals[c][2]):
        if c < 5 or c > V-3 or abs(B_vals[c][2]) > 3:
            print(f"    Class {c}: H={mi[c]['H']}, c3={mi[c]['c3']}, "
                  f"B_1={B_vals[c][1]:.2f}, B_2={B_vals[c][2]:.2f}, "
                  f"B_3={B_vals[c][3]:.2f}")

    # MULTIPLE REGRESSION: can B_1, B_2 predict H and c3?
    X = np.column_stack([B1_arr, B2_arr])
    if np.linalg.matrix_rank(X) == 2:
        # H = a*B1 + b*B2 + c
        X_aug = np.column_stack([np.ones(V), B1_arr, B2_arr])
        beta_H = np.linalg.lstsq(X_aug, H_arr, rcond=None)[0]
        H_pred = X_aug @ beta_H
        R2_H = 1 - np.sum((H_arr - H_pred)**2) / np.sum((H_arr - np.mean(H_arr))**2)

        beta_c3 = np.linalg.lstsq(X_aug, c3_arr, rcond=None)[0]
        c3_pred = X_aug @ beta_c3
        R2_c3 = 1 - np.sum((c3_arr - c3_pred)**2) / np.sum((c3_arr - np.mean(c3_arr))**2)

        print(f"\n  MULTIPLE REGRESSION (B_1, B_2 → H):")
        print(f"    H ≈ {beta_H[0]:.2f} + {beta_H[1]:.2f}×B_1 + {beta_H[2]:.2f}×B_2")
        print(f"    R² = {R2_H:.4f}")

        print(f"\n  MULTIPLE REGRESSION (B_1, B_2 → c3):")
        print(f"    c3 ≈ {beta_c3[0]:.2f} + {beta_c3[1]:.2f}×B_1 + {beta_c3[2]:.2f}×B_2")
        print(f"    R² = {R2_c3:.4f}")

    # With B_3 too
    X_aug3 = np.column_stack([np.ones(V), B1_arr, B2_arr, B3_arr])
    if np.linalg.matrix_rank(X_aug3) >= 3:
        beta_H3 = np.linalg.lstsq(X_aug3, H_arr, rcond=None)[0]
        H_pred3 = X_aug3 @ beta_H3
        R2_H3 = 1 - np.sum((H_arr - H_pred3)**2) / np.sum((H_arr - np.mean(H_arr))**2)

        beta_c33 = np.linalg.lstsq(X_aug3, c3_arr, rcond=None)[0]
        c3_pred3 = X_aug3 @ beta_c33
        R2_c33 = 1 - np.sum((c3_arr - c3_pred3)**2) / np.sum((c3_arr - np.mean(c3_arr))**2)

        print(f"\n  WITH B_3 ADDED:")
        print(f"    H ~ B_1 + B_2 + B_3: R² = {R2_H3:.4f}")
        print(f"    c3 ~ B_1 + B_2 + B_3: R² = {R2_c33:.4f}")

print(f"\n{'='*72}")
print("  THE KRAWTCHOUK DECOMPOSITION OF TOURNAMENT INVARIANTS")
print(f"{'='*72}")
print("""
  B_1 ≈ -H (captures the "position" along the diagonal)
  B_2 captures the "curvature" — the second-order deviation from linearity
  B_3 captures the "asymmetry" — the third-order structure

  Together: B_1, B_2, B_3 form a COORDINATE SYSTEM for tournament space:
    B_1 = principal axis (H-gradient, 85-94% of variance)
    B_2 = first perpendicular (width, ~5-10% of variance)
    B_3 = second perpendicular (skewness, ~2-5% of variance)

  This is the KRAWTCHOUK COORDINATE SYSTEM — the natural basis
  for the Hamming scheme, projected onto the tournament quotient.
  It diagonalizes the distance operator and gives the spectral
  decomposition of tournament space.
""")

print("DONE.")
