#!/usr/bin/env python3
"""
twist_helix_s308.py — The twist structure: is (B1, B2, B3) a helix?
opus-2026-03-24-S308

Plot all iso classes in (B1, B2, B3) Krawtchouk coordinates.
Look for helical/spiral structure as we walk along the H-gradient.
Check if SC and NS form two interleaved strands (double helix).
"""

import sys, subprocess
from math import comb
from itertools import permutations
from collections import defaultdict, Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def krawtchouk(k, x, m):
    total = 0
    for j in range(k+1):
        if j > x or k-j > m-x: continue
        total += (-1)**j * comb(int(x), j) * comb(int(m-x), k-j)
    return total

def build_full(n):
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
        adj = b2a(full_bits); cid = cls[0]['cid'] if False else 0
        from collections import defaultdict as dd2
        def lk2(adj2):
            s2 = ss(adj2); cc2 = c3(adj2); h2 = H(adj2)
            cids2 = hm.get((s2,cc2,h2))
            if not cids2: return -1
            if len(cids2) == 1: return cids2[0]
            return cm.get(cp(adj2), -1)
        cid = lk2(b2a(full_bits))
        tc[tb] = mn[cid] if cid >= 0 else -1
    fiber = Counter()
    for tb in range(total):
        if tc[tb] >= 0: fiber[tc[tb]] += 1
    # Weight enumerators
    weight_enum = defaultdict(lambda: np.zeros(m+1))
    for tb in range(total):
        c = tc[tb]
        if c >= 0:
            w = bin(tb).count('1')
            weight_enum[c][w] += 1
    # Krawtchouk transform
    K = np.zeros((m+1, m+1))
    for kk in range(m+1):
        for x in range(m+1):
            K[kk][x] = krawtchouk(kk, x, m)
    B_vals = {}
    for c in range(V):
        A = weight_enum[c]
        f = fiber.get(c, 0)
        if f > 0:
            B_vals[c] = K @ A / f
    return V, mi, B_vals, m, dict(fiber)

print("=" * 72)
print("  THE TWIST: (B1, B2, B3) STRUCTURE")
print("  opus-2026-03-24-S308")
print("=" * 72)

for n in [5, 6, 7]:
    print(f"\n{'#'*72}")
    print(f"  n = {n}")
    print(f"{'#'*72}")

    V, mi, B_vals, m, fiber = build_full(n)
    print(f"  V={V}, m={m}")

    # Sort classes by H (= -B1)
    sorted_classes = sorted(range(V), key=lambda c: mi[c]['H'])

    # Extract B1, B2, B3 along the H-gradient
    B1 = [B_vals[c][1] for c in sorted_classes]
    B2 = [B_vals[c][2] for c in sorted_classes]
    B3 = [B_vals[c][3] for c in sorted_classes]
    H_sorted = [mi[c]['H'] for c in sorted_classes]
    c3_sorted = [mi[c]['c3'] for c in sorted_classes]
    sc_sorted = [mi[c]['sc'] for c in sorted_classes]

    # Normalize B2 and B3 relative to B1
    # Remove the linear trend (B2 ≈ a*B1 + b) to see the RESIDUAL
    X = np.column_stack([np.ones(V), np.array(B1)])
    beta2 = np.linalg.lstsq(X, np.array(B2), rcond=None)[0]
    B2_residual = np.array(B2) - X @ beta2

    beta3 = np.linalg.lstsq(X, np.array(B3), rcond=None)[0]
    B3_residual = np.array(B3) - X @ beta3

    print(f"\n  B2 = {beta2[0]:.3f} + {beta2[1]:.3f} × B1 + residual")
    print(f"  B3 = {beta3[0]:.3f} + {beta3[1]:.3f} × B1 + residual")
    print(f"  B2 residual std: {np.std(B2_residual):.4f}")
    print(f"  B3 residual std: {np.std(B3_residual):.4f}")

    # THE TWIST: does (B2_residual, B3_residual) rotate as B1 changes?
    # Check if the residual traces a circle/helix
    print(f"\n  TWIST ANALYSIS (B2_resid, B3_resid along H-gradient):")

    # Compute the "angle" of (B2_resid, B3_resid) for each class
    angles = []
    for i in range(V):
        r2 = B2_residual[i]; r3 = B3_residual[i]
        angle = np.arctan2(r3, r2) if (abs(r2) > 0.01 or abs(r3) > 0.01) else 0
        angles.append(angle)

    # Check if angles increase/decrease monotonically with H (= helical rotation)
    angle_diffs = [angles[i+1] - angles[i] for i in range(V-1)]

    positive_diffs = sum(1 for d in angle_diffs if d > 0)
    negative_diffs = sum(1 for d in angle_diffs if d < 0)

    print(f"    Angle increases: {positive_diffs}/{V-1} ({100*positive_diffs/(V-1):.0f}%)")
    print(f"    Angle decreases: {negative_diffs}/{V-1} ({100*negative_diffs/(V-1):.0f}%)")

    # If ~50/50: random, no helix. If >70%: consistent rotation = helix!
    if positive_diffs > 0.7 * (V-1) or negative_diffs > 0.7 * (V-1):
        direction = "CLOCKWISE" if negative_diffs > positive_diffs else "COUNTER-CLOCKWISE"
        print(f"    *** HELICAL ROTATION DETECTED: {direction} ***")
    else:
        print(f"    No consistent rotation (random-like)")

    # The AMPLITUDE of the twist
    amplitudes = [np.sqrt(B2_residual[i]**2 + B3_residual[i]**2) for i in range(V)]
    print(f"\n    Twist amplitude: min={min(amplitudes):.3f}, max={max(amplitudes):.3f}, "
          f"mean={np.mean(amplitudes):.3f}")

    # SC vs NS in the twist
    sc_amp = [amplitudes[i] for i in range(V) if sc_sorted[i]]
    ns_amp = [amplitudes[i] for i in range(V) if not sc_sorted[i]]
    if sc_amp and ns_amp:
        print(f"    SC twist amplitude: {np.mean(sc_amp):.3f}")
        print(f"    NS twist amplitude: {np.mean(ns_amp):.3f}")

    # SC and NS strands: do they separate in the twist plane?
    sc_angles = [angles[i] for i in range(V) if sc_sorted[i]]
    ns_angles = [angles[i] for i in range(V) if not sc_sorted[i]]
    if sc_angles and ns_angles:
        # Circular mean of SC and NS angles
        sc_mean = np.arctan2(np.mean(np.sin(sc_angles)), np.mean(np.cos(sc_angles)))
        ns_mean = np.arctan2(np.mean(np.sin(ns_angles)), np.mean(np.cos(ns_angles)))
        angle_sep = abs(sc_mean - ns_mean)
        if angle_sep > np.pi: angle_sep = 2*np.pi - angle_sep
        print(f"\n    SC mean angle: {sc_mean:.3f} rad ({np.degrees(sc_mean):.1f}°)")
        print(f"    NS mean angle: {ns_mean:.3f} rad ({np.degrees(ns_mean):.1f}°)")
        print(f"    Angular separation: {np.degrees(angle_sep):.1f}°")

        if angle_sep > np.pi/3:  # > 60°
            print(f"    *** SC AND NS ARE ANGULARLY SEPARATED! ***")
            print(f"    → TWO STRANDS in the twist plane")

    # Display (B1, B2_resid, B3_resid) for sorted classes
    if V <= 40:
        print(f"\n    SORTED BY H:")
        print(f"    {'H':>4} {'c3':>3} {'SC':>3} {'B1':>7} {'B2r':>7} {'B3r':>7} {'amp':>7} {'angle':>7}")
        for i, c in enumerate(sorted_classes):
            print(f"    {mi[c]['H']:>4} {mi[c]['c3']:>3} {'SC' if mi[c]['sc'] else 'NS':>3} "
                  f"{B1[i]:>7.2f} {B2_residual[i]:>7.3f} {B3_residual[i]:>7.3f} "
                  f"{amplitudes[i]:>7.3f} {np.degrees(angles[i]):>7.1f}°")

print(f"\n{'='*72}")
print("  THE TWIST STRUCTURE")
print(f"{'='*72}")
print("""
  The TWIST = the residual of (B2, B3) after removing the linear
  dependence on B1 (= H). It captures what H DOESN'T explain.

  If the twist traces a HELIX as H increases:
    → tournament space is a SPIRAL STAIRCASE
    → each H-level is rotated relative to the previous
    → SC and NS form the two rails of the spiral

  If the twist is RANDOM:
    → no coherent rotation
    → the perpendicular structure is unstructured noise

  If SC and NS are ANGULARLY SEPARATED:
    → the two types occupy different sectors of the twist plane
    → this is a DOUBLE HELIX with SC and NS as the two strands
    → like DNA: the two strands carry different information
       (SC = structural, NS = generic) wrapped around the H-axis
""")

print("DONE.")
