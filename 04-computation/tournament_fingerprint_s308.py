#!/usr/bin/env python3
"""
tournament_fingerprint_s308.py -- opus-2026-03-24-S308

TOURNAMENT FINGERPRINT: A COMPACT SPECTRAL HASH

IDEA: Since H is band-limited (hat{H}_k = 0 for k >= m/2),
the iso class of a tournament is determined by its low-frequency
Krawtchouk spectrum. Use the first ceil(m/2) spectral coefficients
as a FINGERPRINT.

PROPERTIES:
  - Size: ceil(m/2) real numbers (vs m bits for the full tiling)
  - Invariant: same fingerprint ↔ same iso class (if band-limitedness is exact)
  - Compact: 50% size reduction
  - Searchable: can use as coordinates for nearest-neighbor search

THIS SCRIPT: implements and verifies the fingerprint at n=5,6,7.
Tests: do fingerprints uniquely identify iso classes?
How many coefficients are actually needed?

Author: opus-2026-03-24-S308
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(n):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(n))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def krawtchouk(k, x, m):
    val = 0
    for j in range(k+1):
        if j > x or k-j > m-x: continue
        val += ((-1)**j) * comb(x, j) * comb(m-x, k-j)
    return val

def tiling_fingerprint(bits, m, n_coeffs):
    """Compute the Krawtchouk fingerprint of a single tiling."""
    wt = bin(bits).count('1')
    return tuple(krawtchouk(k, wt, m) for k in range(n_coeffs))

def class_fingerprint(tilings, m, n_coeffs):
    """Compute the AVERAGE Krawtchouk fingerprint over all tilings in a class."""
    if not tilings:
        return tuple(0.0 for _ in range(n_coeffs))
    fps = [tiling_fingerprint(b, m, n_coeffs) for b in tilings]
    return tuple(np.mean([fp[k] for fp in fps]) for k in range(n_coeffs))

# ============================================================
banner("APPLICATION 1: TOURNAMENT FINGERPRINTING")
# ============================================================

for n in [5, 6, 7]:
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs = [(i,j) for i,j in all_p if (i,j) not in base_set]
    tile_pairs.sort()

    # Build classes
    cmap = {}; tc = {}; members = defaultdict(list)
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs):
            if (bits >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap)
        tc[bits] = cmap[c]; members[cmap[c]].append(bits)
    nc = len(cmap)

    # Merge complement pairs
    info = {}
    for cid in range(nc):
        b = members[cid][0]
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k2, (lo, hi) in enumerate(tile_pairs):
            if (b >> k2) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        comp_c = canon(complement(A, n), n)
        info[cid] = {'comp': cmap.get(comp_c, -1), 'H': Hcount(A, n)}

    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        if comp >= 0 and comp != cid:
            merged[mid] = (cid, comp); mid_map[cid] = mid; mid_map[comp] = mid
        else:
            merged[mid] = (cid,); mid_map[cid] = mid
        mid += 1
    nm = mid

    print("  n=%d: m=%d, %d merged classes\n" % (n, m, nm))

    # Test: how many Krawtchouk coefficients are needed to uniquely identify classes?
    for n_coeffs in range(1, min(m+1, 8)):
        # Compute class fingerprint for each merged class
        class_fps = {}
        for mi in range(nm):
            all_tilings = []
            for c in merged[mi]:
                all_tilings.extend(members[c])
            fp = class_fingerprint(all_tilings, m, n_coeffs)
            # Round to avoid floating point issues
            fp_rounded = tuple(round(x, 6) for x in fp)
            class_fps[mi] = fp_rounded

        # How many UNIQUE fingerprints?
        unique_fps = len(set(class_fps.values()))

        # Are all classes distinguished?
        all_distinct = (unique_fps == nm)

        # Collisions?
        fp_to_classes = defaultdict(list)
        for mi, fp in class_fps.items():
            fp_to_classes[fp].append(mi)
        n_collisions = sum(1 for fp, cls in fp_to_classes.items() if len(cls) > 1)

        status = "✓ PERFECT" if all_distinct else "✗ %d collisions" % n_collisions
        print("    %d coefficients: %d unique / %d classes  %s" %
              (n_coeffs, unique_fps, nm, status))

    # Show the fingerprints for each class (first few coefficients)
    if n <= 6:
        print("\n    CLASS FINGERPRINTS (first 5 coefficients):\n")
        n_show = min(5, m)
        print("    %-4s %-4s %-8s" % ("mid", "H", "fiber"), end="")
        for k in range(n_show):
            print(" K_%d      " % k, end="")
        print()
        print("    " + "-"*(30 + 10*n_show))

        order = sorted(range(nm), key=lambda x: info[merged[x][0]]['H'])
        for mi in order:
            H = info[merged[mi][0]]['H']
            all_tilings = []
            for c in merged[mi]:
                all_tilings.extend(members[c])
            fiber = len(all_tilings)
            fp = class_fingerprint(all_tilings, m, n_show)
            print("    %-4d %-4d %-8d" % (mi, H, fiber), end="")
            for k in range(n_show):
                print(" %8.3f " % fp[k], end="")
            print()

    print()

# ============================================================
banner("APPLICATION 2: ERROR DETECTION")
# ============================================================

# A valid tournament tiling has zero high-frequency energy.
# Inject random bit errors and check if the spectrum detects them.

n = 5; m = comb(n-1, 2); total = 2**m
base_set = set((i, i+1) for i in range(n-1))
all_p2 = [(i,j) for i in range(n) for j in range(i+1,n)]
tile_pairs2 = [(i,j) for i,j in all_p2 if (i,j) not in base_set]
tile_pairs2.sort()

import random
random.seed(42)

print("  ERROR DETECTION via high-frequency spectrum (n=5, m=6)\n")

# For each valid tiling, inject 1-3 bit errors, check if detected
n_trials = 200
detect_counts = {1: 0, 2: 0, 3: 0}
total_counts = {1: 0, 2: 0, 3: 0}

for _ in range(n_trials):
    # Pick a random valid tiling
    bits = random.randint(0, total-1)

    for n_errors in [1, 2, 3]:
        # Inject errors at random positions
        error_positions = random.sample(range(m), n_errors)
        corrupted = bits
        for p in error_positions:
            corrupted ^= (1 << p)

        # Compute high-frequency energy of corrupted tiling
        wt_orig = bin(bits).count('1')
        wt_corr = bin(corrupted).count('1')

        # High-freq energy: sum of hat{H}_k^2 for k > m//2
        high_energy = 0
        for k in range(m//2 + 1, m+1):
            # The tiling-level Krawtchouk value at the corrupted weight
            Kk = krawtchouk(k, wt_corr, m)
            high_energy += Kk ** 2

        # Also check: is the corrupted word in a DIFFERENT iso class?
        A_orig = [[0]*n for _ in range(n)]
        for i in range(1, n): A_orig[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs2):
            if (bits >> k) & 1: A_orig[lo][hi] = 1
            else: A_orig[hi][lo] = 1

        A_corr = [[0]*n for _ in range(n)]
        for i in range(1, n): A_corr[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs2):
            if (corrupted >> k) & 1: A_corr[lo][hi] = 1
            else: A_corr[hi][lo] = 1

        same_class = (canon(A_orig, n) == canon(A_corr, n))
        class_changed = not same_class

        # Detection: high-frequency energy differs from original
        high_orig = sum(krawtchouk(k, wt_orig, m)**2 for k in range(m//2+1, m+1))

        detected = (abs(high_energy - high_orig) > 0.1) or (wt_corr != wt_orig)

        total_counts[n_errors] += 1
        if detected and class_changed:
            detect_counts[n_errors] += 1

print("  %-10s %-15s %-15s %-15s" %
      ("Errors", "Class changed", "Detected", "Detection rate"))
print("  " + "-"*55)

for ne in [1, 2, 3]:
    # Recount properly
    class_changed_count = 0
    detected_count = 0
    for _ in range(n_trials):
        bits = random.randint(0, total-1)
        error_positions = random.sample(range(m), ne)
        corrupted = bits
        for p in error_positions:
            corrupted ^= (1 << p)

        A_orig = [[0]*n for _ in range(n)]
        for i in range(1, n): A_orig[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs2):
            if (bits >> k) & 1: A_orig[lo][hi] = 1
            else: A_orig[hi][lo] = 1

        A_corr = [[0]*n for _ in range(n)]
        for i in range(1, n): A_corr[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs2):
            if (corrupted >> k) & 1: A_corr[lo][hi] = 1
            else: A_corr[hi][lo] = 1

        changed = (canon(A_orig, n) != canon(A_corr, n))
        if changed:
            class_changed_count += 1

        # Simple detection: weight changed
        wt_changed = (bin(bits).count('1') != bin(corrupted).count('1'))
        if wt_changed and changed:
            detected_count += 1

    rate = detected_count / class_changed_count if class_changed_count > 0 else 0
    print("  %-10d %-15d %-15d %-15.1f%%" %
          (ne, class_changed_count, detected_count, 100*rate))

# ============================================================
banner("APPLICATION 3: COMPRESSION RATIO MEASUREMENT")
# ============================================================

# Measure ACTUAL compression: how many bits to identify a class
# using Krawtchouk fingerprint vs full tiling?

for n in [5, 6]:
    m = comb(n-1, 2)
    total = 2**m
    base_set = set((i, i+1) for i in range(n-1))
    all_p3 = [(i,j) for i in range(n) for j in range(i+1,n)]
    tile_pairs3 = [(i,j) for i,j in all_p3 if (i,j) not in base_set]
    tile_pairs3.sort()

    cmap = {}; tc = {}; members = defaultdict(list)
    for bits in range(total):
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs3):
            if (bits >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap)
        tc[bits] = cmap[c]; members[cmap[c]].append(bits)
    nc = len(cmap)

    # Merge
    inf2 = {}
    for cid in range(nc):
        b = members[cid][0]
        A = [[0]*n for _ in range(n)]
        for i in range(1, n): A[i][i-1] = 1
        for k, (lo, hi) in enumerate(tile_pairs3):
            if (b >> k) & 1: A[lo][hi] = 1
            else: A[hi][lo] = 1
        comp_c = canon(complement(A, n), n)
        inf2[cid] = {'comp': cmap.get(comp_c, -1)}
    merged2 = {}; mid_map2 = {}; mid2 = 0
    for cid in range(nc):
        if cid in mid_map2: continue
        comp = inf2[cid]['comp']
        if comp >= 0 and comp != cid:
            merged2[mid2] = (cid, comp); mid_map2[cid] = mid2; mid_map2[comp] = mid2
        else:
            merged2[mid2] = (cid,); mid_map2[cid] = mid2
        mid2 += 1
    nm2 = mid2

    # Minimum coefficients needed for unique identification
    for n_coeffs in range(1, m+1):
        class_fps = {}
        for mi in range(nm2):
            all_tilings = []
            for c in merged2[mi]:
                all_tilings.extend(members[c])
            fp = class_fingerprint(all_tilings, m, n_coeffs)
            fp_rounded = tuple(round(x, 4) for x in fp)
            class_fps[mi] = fp_rounded
        if len(set(class_fps.values())) == nm2:
            min_coeffs = n_coeffs
            break
    else:
        min_coeffs = m

    from math import log2, ceil
    bits_naive = m  # full tiling
    bits_class = ceil(log2(nm2)) if nm2 > 1 else 1  # class index
    bits_fingerprint = min_coeffs * 8  # 8 bits per coefficient (quantized)

    print("  n=%d (m=%d, %d classes):" % (n, m, nm2))
    print("    Full tiling: %d bits" % bits_naive)
    print("    Class index: %d bits (log₂(%d))" % (bits_class, nm2))
    print("    Fingerprint: %d coefficients needed, ≈%d bits (quantized)" %
          (min_coeffs, bits_fingerprint))
    print("    Compression: %d → %d coefficients = %.1f:1" %
          (m, min_coeffs, m/min_coeffs))
    print()

print("Done. opus-2026-03-24-S308")
