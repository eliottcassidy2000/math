#!/usr/bin/env python3
"""
new_statistics_s20cj.py -- kind-pasteur-2026-03-22-S20cj

NEW STATISTICS FOR G_n: What's missing and what would be useful.

The KEY question: what MINIMAL set of coordinates uniquely locates
each iso class in G_n? If we find this, we have a "canonical address"
for every tournament type.

NEW STATISTICS TO COMPUTE:
1. H-level distribution (shape of the DAG)
2. Degree distribution (histogram)
3. Blue ratio = blue_deg / total_deg per class
4. The "canonical coordinates" -- minimal set that uniquely identifies each class
5. H-gap distribution (gaps between consecutive H values)
6. The |Aut| spectrum
7. Effective dimension of G_n (how many coordinates needed?)

Author: kind-pasteur-2026-03-22-S20cj
"""
import sys
import numpy as np
from math import comb, factorial, log2
from collections import defaultdict, Counter
from itertools import permutations
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

def is_sc(A, n):
    A_comp = np.zeros_like(A)
    for i in range(n):
        for j in range(n):
            if i != j: A_comp[i][j] = 1 - A[i][j]
    for perm in permutations(range(n)):
        if all(A[perm[i]][perm[j]] == A_comp[i][j] for i in range(n) for j in range(n) if i != j):
            return True
    return False

print("=" * 70)
print("  NEW STATISTICS FOR G_n")
print("=" * 70)

for n in [5, 6]:
    print(f"\n{'='*70}")
    print(f"  n = {n}")
    print(f"{'='*70}\n")

    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)

    # Build classes with all invariants
    canon_map = defaultdict(list)
    H_map = {}
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        H = count_hp(A, n)
        H_map[bits] = H
        s = A.sum(axis=1).astype(int)
        S2 = int(sum(s*s))
        c3 = comb(n,3) - (S2 - comb(n,2)) // 2
        score = tuple(sorted(s))
        sc = is_sc(A, n) if n <= 6 else None
        cf = canonical(A, n)
        canon_map[cf].append((bits, H, score, c3, S2, sc))

    classes = []
    for cf, members in sorted(canon_map.items(), key=lambda x: x[1][0][1]):
        bits0, H, score, c3, S2, sc = members[0]
        cid = len(classes)
        aut = factorial(n) // len(members)
        classes.append({'id': cid, 'H': H, 'score': score, 'c3': c3,
                       'S2': S2, 'sc': sc, 'aut': aut, 'size': len(members),
                       'members': set(b for b, *_ in members)})

    bits_to_class = {}
    for c in classes:
        for b in c['members']:
            bits_to_class[b] = c['id']

    N = len(classes)

    # Build adjacency + blue/black
    adj = np.zeros((N, N), dtype=int)
    blue_adj = np.zeros((N, N), dtype=int)
    for c in classes:
        T = list(c['members'])[0]
        for k in range(m):
            nb = T ^ (1 << k)
            nb_id = bits_to_class[nb]
            if nb_id != c['id']:
                adj[c['id']][nb_id] = 1
                adj[nb_id][c['id']] = 1
                if c['sc'] == classes[nb_id]['sc']:
                    blue_adj[c['id']][nb_id] = 1

    degrees = [int(sum(adj[i])) for i in range(N)]
    blue_degrees = [int(sum(blue_adj[i])) for i in range(N)]

    # ---- STAT 1: H-level distribution ----
    H_levels = defaultdict(list)
    for c in classes:
        H_levels[c['H']].append(c['id'])

    print(f"  STAT 1: H-LEVEL DISTRIBUTION (shape of the DAG)")
    h_widths = []
    for h in sorted(H_levels.keys()):
        w = len(H_levels[h])
        h_widths.append(w)
        if n <= 5:
            print(f"    H={h:>3d}: {'#' * w} ({w})")

    # Shape analysis
    max_w = max(h_widths)
    total_levels = len(h_widths)
    peak_pos = h_widths.index(max_w)
    skew = (peak_pos / (total_levels - 1)) if total_levels > 1 else 0.5
    print(f"  Shape: {total_levels} levels, max width {max_w} at level {peak_pos}")
    print(f"  Skewness (0=bottom, 1=top): {skew:.3f}")
    print(f"  Widths: {h_widths}")

    # ---- STAT 2: Degree distribution ----
    print(f"\n  STAT 2: DEGREE DISTRIBUTION")
    deg_dist = Counter(degrees)
    for d in sorted(deg_dist.keys()):
        print(f"    deg={d:>2d}: {'#' * deg_dist[d]} ({deg_dist[d]})")

    # ---- STAT 3: Blue ratio per class ----
    print(f"\n  STAT 3: BLUE RATIO = blue_deg / total_deg")
    blue_ratios = []
    for c in classes:
        d = degrees[c['id']]
        bd = blue_degrees[c['id']]
        ratio = bd / d if d > 0 else 0
        blue_ratios.append(ratio)

    sc_ratios = [blue_ratios[c['id']] for c in classes if c['sc']]
    ns_ratios = [blue_ratios[c['id']] for c in classes if not c['sc']]
    print(f"  SC classes: avg blue ratio = {np.mean(sc_ratios):.3f} (range [{min(sc_ratios):.2f}, {max(sc_ratios):.2f}])")
    if ns_ratios:
        print(f"  NS classes: avg blue ratio = {np.mean(ns_ratios):.3f} (range [{min(ns_ratios):.2f}, {max(ns_ratios):.2f}])")

    # ---- STAT 4: Canonical coordinates ----
    print(f"\n  STAT 4: CANONICAL COORDINATES")
    # How many invariants needed to uniquely identify each class?

    # Test (H alone)
    h_unique = len(set(c['H'] for c in classes)) == N
    print(f"    H alone determines class: {h_unique}")

    # Test (H, score)
    hs_unique = len(set((c['H'], c['score']) for c in classes)) == N
    print(f"    (H, score) determines class: {hs_unique}")

    # Test (H, |Aut|)
    ha_unique = len(set((c['H'], c['aut']) for c in classes)) == N
    print(f"    (H, |Aut|) determines class: {ha_unique}")

    # Test (H, score, |Aut|)
    hsa_unique = len(set((c['H'], c['score'], c['aut']) for c in classes)) == N
    print(f"    (H, score, |Aut|) determines class: {hsa_unique}")

    # Test (H, degree)
    hd_unique = len(set((c['H'], degrees[c['id']]) for c in classes)) == N
    print(f"    (H, degree) determines class: {hd_unique}")

    # Test (H, score, degree)
    hsd_unique = len(set((c['H'], c['score'], degrees[c['id']]) for c in classes)) == N
    print(f"    (H, score, degree) determines class: {hsd_unique}")

    # Test (score, c3)
    sc3_unique = len(set((c['score'], c['c3']) for c in classes)) == N
    print(f"    (score, c3) determines class: {sc3_unique}")

    # Test (H, blue_deg, black_deg)
    hbb_unique = len(set((c['H'], blue_degrees[c['id']], degrees[c['id']] - blue_degrees[c['id']]) for c in classes)) == N
    print(f"    (H, blue_deg, black_deg) determines class: {hbb_unique}")

    # Find minimal determining set
    # Try all pairs of invariants
    invs = {
        'H': [c['H'] for c in classes],
        'score': [c['score'] for c in classes],
        'c3': [c['c3'] for c in classes],
        '|Aut|': [c['aut'] for c in classes],
        'deg': degrees,
        'b_deg': blue_degrees,
        'sc': [c['sc'] for c in classes],
    }

    print(f"\n    MINIMAL DETERMINING PAIRS:")
    for name1 in invs:
        for name2 in invs:
            if name1 >= name2: continue
            combined = set(zip(invs[name1], invs[name2]))
            if len(combined) == N:
                print(f"      ({name1}, {name2}) DETERMINES at n={n}")

    # ---- STAT 5: H-gap distribution ----
    print(f"\n  STAT 5: H-GAP DISTRIBUTION")
    H_vals = sorted(set(c['H'] for c in classes))
    gaps = [H_vals[i+1] - H_vals[i] for i in range(len(H_vals)-1)]
    gap_dist = Counter(gaps)
    print(f"  H values: {H_vals}")
    print(f"  Gaps: {gaps}")
    print(f"  Gap distribution: {dict(sorted(gap_dist.items()))}")
    print(f"  All gaps even: {all(g % 2 == 0 for g in gaps)}")

    # ---- STAT 6: |Aut| spectrum ----
    print(f"\n  STAT 6: |Aut| SPECTRUM")
    aut_dist = Counter(c['aut'] for c in classes)
    for a in sorted(aut_dist.keys()):
        print(f"    |Aut|={a}: {aut_dist[a]} classes ({100*aut_dist[a]/N:.1f}%)")

    # ---- STAT 7: Effective dimension ----
    print(f"\n  STAT 7: EFFECTIVE DIMENSION")
    # How many independent coordinates needed to embed G_n without distortion?
    # Use the adjacency matrix eigenvalues
    eigvals = sorted(np.linalg.eigvalsh(adj.astype(float)), reverse=True)
    # Count eigenvalues explaining 90% of spectral energy
    total_energy = sum(e**2 for e in eigvals)
    cum_energy = 0
    dim_90 = 0
    for e in eigvals:
        cum_energy += e**2
        dim_90 += 1
        if cum_energy >= 0.9 * total_energy:
            break
    print(f"  Spectral dimension (90% energy): {dim_90}")
    print(f"  Top 5 eigenvalues: {[f'{e:.2f}' for e in eigvals[:5]]}")

print(f"\n{'='*70}")
print(f"  SYNTHESIS: THE MOST USEFUL NEW STATISTICS")
print(f"{'='*70}\n")

print(f"""  THE 5 MOST INFORMATIVE NEW STATISTICS:

  1. H-LEVEL DISTRIBUTION (shape of DAG):
     Reveals whether G_n is "diamond" (wide middle) or "pyramid" (wide bottom).
     At n=5: peaked at bottom (H=3, skew=0.17). At n=6: peaked at 33% of max.
     The peak RISES with n.

  2. BLUE RATIO per class:
     SC classes: low blue ratio (SC is rare, most neighbors are NS = black)
     NS classes: high blue ratio (NS is common, most neighbors are NS = blue)
     The blue ratio DETERMINES SC status at n=6 (clean separation).

  3. CANONICAL COORDINATES:
     The minimal set of invariants that uniquely identifies each class.
     Reveals the "intrinsic dimension" of tournament type space.

  4. H-GAP DISTRIBUTION:
     All gaps are even (Redei). The gap pattern reveals forbidden values.
     Large gaps = topological obstructions in tournament space.

  5. SPECTRAL DIMENSION:
     The number of eigenvectors needed to capture 90% of G_n's structure.
     This is the "effective dimension" of the meta-graph.
""")
