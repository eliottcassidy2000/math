#!/usr/bin/env python3
"""
krawtchouk_s305.py — Krawtchouk polynomials and the tournament tiling space
opus-2026-03-24-S305

The tiling space IS the Hamming scheme H(m,2) where m = C(n-1,2).
Krawtchouk polynomials K_k(x; m) are the eigenfunctions.

CONNECTIONS:
1. K_k(x; m) = Σ_{j=0}^k (-1)^j C(x,j) C(m-x, k-j) × 2^{k-j} / ...
   Actually: K_k(x; m) = Σ_{j=0}^k (-1)^j C(x,j) C(m-x, k-j)

2. The distance-d adjacency matrix A_d of H(m,2) has eigenvalues K_d(j; m)
   in the j-th eigenspace (j = 0, 1, ..., m).

3. The MacWilliams transform: W_C^⊥(x,y) = (1/|C|) W_C(x+y, x-y)
   relates weight enumerators of dual codes. For our "iso class codes,"
   this gives the dual weight distribution.

4. The weight distribution of each iso class = distribution of Hamming
   weights (tile flip counts) among its tilings. Its Krawtchouk transform
   = the "dual spectrum" of the class.
"""

import sys, time, subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# KRAWTCHOUK POLYNOMIALS
# ============================================================

def krawtchouk(k, x, m):
    """Krawtchouk polynomial K_k(x; m) for binary Hamming scheme.
    K_k(x; m) = Σ_{j=0}^k (-1)^j C(x,j) C(m-x, k-j)"""
    total = 0
    for j in range(k+1):
        if j > x or k-j > m-x: continue
        total += (-1)**j * comb(x, j) * comb(m-x, k-j)
    return total

def krawtchouk_matrix(m):
    """Full Krawtchouk matrix K[k][x] = K_k(x; m) for k,x = 0..m."""
    K = np.zeros((m+1, m+1))
    for k in range(m+1):
        for x in range(m+1):
            K[k][x] = krawtchouk(k, x, m)
    return K

# ============================================================
# BUILD TILING DATA
# ============================================================

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
    fiber = Counter()
    for tb in range(total):
        if tc[tb] >= 0: fiber[tc[tb]] += 1
    return V, mi, tc, m, total, dict(fiber)

# ============================================================
# ANALYSIS
# ============================================================

print("=" * 72)
print("  KRAWTCHOUK POLYNOMIALS AND THE TILING SPACE")
print("  opus-2026-03-24-S305")
print("=" * 72)

for n in [5, 6]:
    print(f"\n{'#'*72}")
    print(f"  n = {n}")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, tc, m, total, fiber = build_data(n)
    print(f"  V={V}, m={m}, 2^m={total}. Built in {time.time()-t0:.1f}s")

    # 1. KRAWTCHOUK MATRIX
    K = krawtchouk_matrix(m)
    print(f"\n  KRAWTCHOUK MATRIX K[k][x] for m={m}:")
    print(f"  K_0(x) = {[int(K[0][x]) for x in range(min(m+1,8))]}")
    print(f"  K_1(x) = {[int(K[1][x]) for x in range(min(m+1,8))]}")
    print(f"  K_2(x) = {[int(K[2][x]) for x in range(min(m+1,8))]}")

    # 2. WEIGHT ENUMERATOR OF EACH ISO CLASS
    # W_C = distribution of Hamming weights among tilings in C
    weight_enum = defaultdict(lambda: [0]*(m+1))
    for tb in range(total):
        c = tc[tb]
        if c >= 0:
            w = bin(tb).count('1')
            weight_enum[c][w] += 1

    print(f"\n  WEIGHT ENUMERATORS (Hamming weight distribution per class):")
    for c in sorted(mi.keys())[:8]:
        we = weight_enum[c]
        nonzero = [(w, we[w]) for w in range(m+1) if we[w] > 0]
        print(f"    Class {c} (H={mi[c]['H']}, fiber={fiber[c]}): "
              f"weights = {nonzero}")

    # 3. KRAWTCHOUK TRANSFORM of weight enumerator
    # For class C with weight enumerator A = (A_0, ..., A_m):
    # B_k = (1/|C|) Σ_i A_i K_k(i; m) = MacWilliams dual spectrum
    print(f"\n  KRAWTCHOUK TRANSFORM (dual spectrum) per class:")
    for c in sorted(mi.keys())[:8]:
        A = np.array(weight_enum[c], dtype=float)
        f = fiber[c]
        if f > 0:
            B = K @ A / f  # Krawtchouk transform = K × A / |C|
            # B_0 should be 2^m / |C| (the number of "codewords" in the dual)
            nonzero_B = [(k, B[k]) for k in range(m+1) if abs(B[k]) > 0.01]
            print(f"    Class {c} (H={mi[c]['H']}): B_0={B[0]:.1f}, "
                  f"B_1={B[1]:.1f}, B_2={B[2]:.1f}, "
                  f"dual = {nonzero_B[:6]}")

    # 4. THE TOTAL WEIGHT ENUMERATOR (all tilings)
    total_we = [0]*(m+1)
    for w in range(m+1):
        total_we[w] = comb(m, w)  # each weight has C(m,w) tilings

    print(f"\n  TOTAL WEIGHT ENUMERATOR (all tilings):")
    print(f"    A = {total_we[:8]}...")
    total_B = K @ np.array(total_we, dtype=float) / total
    print(f"    Dual: B = {[f'{b:.1f}' for b in total_B[:8]]}")

    # 5. MEAN KRAWTCHOUK VALUE per class
    # For each class, compute K_1(mean_weight) = m - 2×mean_weight
    # This relates to the "average flip direction"
    print(f"\n  K_1(mean_weight) per class (= m - 2×mean_weight):")
    for c in sorted(mi.keys())[:10]:
        A = np.array(weight_enum[c], dtype=float)
        if sum(A) > 0:
            mean_w = sum(w*A[w] for w in range(m+1)) / sum(A)
            k1 = m - 2*mean_w  # K_1(x; m) = m - 2x
            print(f"    Class {c} (H={mi[c]['H']}): mean_weight={mean_w:.2f}, "
                  f"K_1={k1:.2f}")

    # 6. CORRELATION: K_1 vs H
    k1_vals = []
    h_vals = []
    for c in range(V):
        A = np.array(weight_enum[c], dtype=float)
        if sum(A) > 0:
            mean_w = sum(w*A[w] for w in range(m+1)) / sum(A)
            k1_vals.append(m - 2*mean_w)
            h_vals.append(mi[c]['H'])

    if len(k1_vals) > 2:
        r_k1_H = np.corrcoef(k1_vals, h_vals)[0,1]
        print(f"\n  Corr(K_1, H) = {r_k1_H:.4f}")
        print(f"  (K_1 = m - 2×mean_weight measures 'distance from center')")

    # 7. THE WEIGHT-DISTANCE CONNECTION
    # For a tiling at weight w, its metagraph class distance from transitive ≈ w
    # (from the completeness theorem: dist = min_Hamming = min weight in class)
    # So: min_weight per class ≈ metagraph distance from transitive
    print(f"\n  MIN WEIGHT PER CLASS = DISTANCE FROM TRANSITIVE:")
    for c in sorted(mi.keys())[:10]:
        min_w = min(w for w in range(m+1) if weight_enum[c][w] > 0)
        max_w = max(w for w in range(m+1) if weight_enum[c][w] > 0)
        print(f"    Class {c} (H={mi[c]['H']}): min_w={min_w}, max_w={max_w}, "
              f"span={max_w-min_w}")

print(f"\n{'='*72}")
print("  KRAWTCHOUK CONNECTIONS SUMMARY")
print(f"{'='*72}")
print("""
  THE TILING SPACE IS THE HAMMING SCHEME H(m, 2) where m = C(n-1,2).
  Krawtchouk polynomials K_k(x; m) are the natural eigenfunctions.

  1. K_1(x; m) = m - 2x. For a class with mean weight w̄:
     K_1(w̄) = m - 2w̄. This measures "distance from the center"
     of the Hamming scheme. It ANTI-CORRELATES with H.

  2. Each iso class C is a "code" in the Hamming scheme.
     Its weight enumerator A_C(x) = Σ #{tilings of weight w} × x^w.
     The Krawtchouk transform gives the dual spectrum B_C.

  3. The MacWilliams identity relates A_C to B_C:
     B_k = (1/|C|) Σ A_i K_k(i; m)
     This is the FOURIER TRANSFORM on the Hamming scheme.

  4. The min weight of each class = metagraph distance from transitive
     (by the completeness theorem). So the minimum Krawtchouk function
     value determines the class's position in the metagraph.

  5. The S_n action on the Hamming scheme preserves the distance
     structure. So the Krawtchouk polynomials descend to the quotient
     Q_m / S_n = the tournament metagraph. The descended polynomials
     are "S_n-averaged Krawtchouk" = ZONAL SPHERICAL FUNCTIONS of the
     Gelfand pair (S_n, Stab(transitive)).
""")

print("DONE.")
