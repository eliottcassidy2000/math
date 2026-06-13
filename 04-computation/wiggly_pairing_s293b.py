#!/usr/bin/env python3
"""
wiggly_pairing_s293b.py -- opus-2026-03-24-S293

THE 2-FOLD PAIRING IN WIGGLY PROFILES

From S293: profile orbits are almost always size 2.
HYPOTHESIS: the pairing comes from PATH REVERSAL.

If HP = (v₁, v₂, ..., vₙ), the reversed HP = (vₙ, ..., v₂, v₁).
The reversed HP views the SAME tournament through the "opposite" end.

For the STANDARD base path (n-1, n-2, ..., 1, 0):
The reversed path is (0, 1, ..., n-2, n-1).
These give DIFFERENT tilings of the SAME tournament.

If T has the standard path as an HP, does it have the reversed path too?
YES iff the standard path visits vertices in BOTH directions — i.e., iff
there's also an HP going 0 → 1 → ... → n-1.

Actually: given HP P = (v₁,...,vₙ), the reversed HP is P^rev = (vₙ,...,v₁).
P^rev is an HP of T iff arc(vₙ, vₙ₋₁), arc(vₙ₋₁, vₙ₋₂), ..., arc(v₂, v₁) all exist.
But P is an HP means arc(v₁, v₂), ..., arc(vₙ₋₁, vₙ) exist.
P^rev is an HP means arc(vₙ, vₙ₋₁), ..., arc(v₂, v₁) exist.
These are the OPPOSITE directions of the same arcs!

So P^rev is an HP of T iff the base-path arcs of P are ALL reversed in T.
That means P^rev is an HP of T^op (the complement), NOT of T.

Wait: P^rev visits the same vertices but traverses the arcs backwards.
If arc(v₁→v₂) is in T, then arc(v₂→v₁) is in T^op.
So P^rev is an HP of T^op, not T.

THE PAIRING IS NOT PATH REVERSAL (that would give T^op).

Let me think again about what creates the 2-fold pairing...

THE STAIRCASE REFLECTION: the staircase has a symmetry
(x,y) → (n-y+1, n-x+1) = the anti-diagonal reflection.
This maps tile (x,y) to tile (n-y+1, n-x+1).

In the base-path model: if we relabel vertex v → (n-1-v),
we get a different HP of the SAME tournament (if T is invariant under
this relabeling... which it isn't in general).

Hmm. Let me just LOOK at the data to find the pairing mechanism.

Author: opus-2026-03-24-S293
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

def find_all_HPs(A, n):
    paths = []
    def dfs(path, visited):
        if len(path) == n:
            paths.append(tuple(path))
            return
        v = path[-1]
        for w in range(n):
            if w not in visited and A[v][w]:
                dfs(path + [w], visited | {w})
    for start in range(n):
        dfs([start], {start})
    return paths

# ============================================================
banner("THE 2-FOLD PAIRING: WHAT CREATES IT? (n=5)")
# ============================================================

n = 5
m = comb(n-1, 2)
total = 2**m

tile_pairs = []
for hi in range(n):
    for lo in range(hi):
        if hi - lo >= 2:
            tile_pairs.append((hi, lo))
tile_pairs.sort()

# Build classes
cmap = {}; sizes = {}; tc = {}; members = defaultdict(list)
for bits in range(total):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    for k, (hi, lo) in enumerate(tile_pairs):
        if (bits >> k) & 1:
            A[lo][hi] = 1
        else:
            A[hi][lo] = 1
    c = canon(A, n)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
    cid = cmap[c]
    sizes[cid] += 1; tc[bits] = cid
    members[cid].append(bits)

nc = len(cmap)
info = {}
for cid in range(nc):
    bits0 = members[cid][0]
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    for k, (hi, lo) in enumerate(tile_pairs):
        if (bits0 >> k) & 1:
            A[lo][hi] = 1
        else:
            A[hi][lo] = 1
    comp_c = canon(complement(A, n), n)
    info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c], 'SC': cmap[comp_c] == cid, 'size': sizes[cid]}

merged = {}; mid_map = {}; mid = 0
for cid in range(nc):
    if cid in mid_map: continue
    comp = info[cid]['comp']
    merged[mid] = (cid,) if comp == cid else (cid, comp)
    mid_map[cid] = mid
    if comp != cid: mid_map[comp] = mid
    mid += 1
nm = mid

# For each class with fiber > 1: find the paired tilings
print("  For each class, find what PAIRS tilings with the same class-count profile.\n")

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid0 = merged[mi][0]
    H = info[cid0]['H']
    all_tilings = []
    for c in merged[mi]:
        all_tilings.extend(members[c])
    fiber = len(all_tilings)
    if fiber <= 2: continue

    # Get class-count profile for each tiling
    tiling_profiles = {}
    for bits in all_tilings:
        profile = Counter()
        for k in range(m):
            mj = mid_map[tc[bits ^ (1 << k)]]
            profile[mj] += 1
        tiling_profiles[bits] = tuple(sorted(profile.items()))

    # Find pairs with same profile
    profile_groups = defaultdict(list)
    for bits, prof in tiling_profiles.items():
        profile_groups[prof].append(bits)

    print("  mid=%d (H=%d, fiber=%d):" % (mi, H, fiber))

    for prof, group in sorted(profile_groups.items(), key=lambda x: -len(x[1])):
        if len(group) == 1:
            continue  # Singleton, no pairing

        # For each pair: what is the RELATIONSHIP between the two tilings?
        if len(group) == 2:
            b1, b2 = group
            xor = b1 ^ b2
            hamming = bin(xor).count('1')
            diff_positions = [k for k in range(m) if (xor >> k) & 1]
            diff_tiles = [tile_pairs[k] for k in diff_positions]

            # Which tiles are shared neutral?
            n1 = set(k for k in range(m) if mid_map[tc[b1 ^ (1<<k)]] == mi)
            n2 = set(k for k in range(m) if mid_map[tc[b2 ^ (1<<k)]] == mi)
            shared_neutral = n1 & n2
            unique_n1 = n1 - n2
            unique_n2 = n2 - n1

            # What's the complement relationship?
            comp_bits = b1 ^ ((1 << m) - 1)
            is_complement = (comp_bits == b2)

            # What about the LABELED tournament relationship?
            # Build tournament for b1
            A1 = [[0]*n for _ in range(n)]
            for i in range(1, n): A1[i][i-1] = 1
            for k, (hi, lo) in enumerate(tile_pairs):
                if (b1 >> k) & 1: A1[lo][hi] = 1
                else: A1[hi][lo] = 1

            # Build tournament for b2
            A2 = [[0]*n for _ in range(n)]
            for i in range(1, n): A2[i][i-1] = 1
            for k, (hi, lo) in enumerate(tile_pairs):
                if (b2 >> k) & 1: A2[lo][hi] = 1
                else: A2[hi][lo] = 1

            # Check if A2 = σ(A1) for some permutation σ
            relabeling = None
            for perm in permutations(range(n)):
                ok = True
                for i in range(n):
                    for j in range(n):
                        if i != j and A1[i][j] != A2[perm[i]][perm[j]]:
                            ok = False; break
                    if not ok: break
                if ok:
                    relabeling = perm
                    break

            # Check: does σ map the base path to another HP?
            bp1 = list(range(n-1, -1, -1))  # standard: 4,3,2,1,0
            if relabeling:
                bp2 = [relabeling[v] for v in bp1]
                # Is bp2 a valid HP of A1?
                is_hp = all(A1[bp2[i]][bp2[i+1]] for i in range(n-1))
            else:
                bp2 = None; is_hp = False

            print("    PAIR: bits=%s vs %s (hamming=%d, complement=%s)" %
                  (format(b1, '06b'), format(b2, '06b'), hamming, is_complement))
            if relabeling:
                print("      relabeling σ=%s maps base path [4,3,2,1,0] → %s (valid HP: %s)" %
                      (relabeling, bp2, is_hp))
                # What does σ look like?
                cycles = []
                visited = [False]*n
                for i in range(n):
                    if visited[i]: continue
                    cycle = []
                    j = i
                    while not visited[j]:
                        visited[j] = True
                        cycle.append(j)
                        j = relabeling[j]
                    if len(cycle) > 1:
                        cycles.append(tuple(cycle))
                if cycles:
                    print("      σ cycles: %s" % cycles)
                else:
                    print("      σ = identity (?)")

        elif len(group) > 2:
            # Larger orbit — compute all pairwise relationships
            hamming_dists = []
            for i in range(len(group)):
                for j in range(i+1, len(group)):
                    hamming_dists.append(bin(group[i] ^ group[j]).count('1'))
            print("    GROUP of %d: hamming distances = %s" %
                  (len(group), sorted(hamming_dists)))

# ============================================================
banner("THE REGULAR TOURNAMENT (mid=9, H=15, |Aut|=5)")
# ============================================================

# This class has fiber=3, and ALL 3 profiles are distinct.
# The 15 HPs / 5 Aut = 3 orbits.
# Let's trace exactly what happens.

mi = 9
cid0 = merged[mi][0]
bits0 = members[cid0][0]
A = [[0]*n for _ in range(n)]
for i in range(1, n): A[i][i-1] = 1
for k, (hi, lo) in enumerate(tile_pairs):
    if (bits0 >> k) & 1: A[lo][hi] = 1
    else: A[hi][lo] = 1

print("  THE REGULAR TOURNAMENT (unique at n=5, scores=(2,2,2,2,2)):\n")
print("  Adjacency matrix:")
for i in range(n):
    print("    %s" % A[i])

# All HPs
hps = find_all_HPs(A, n)
print("\n  %d Hamiltonian paths:" % len(hps))

# Automorphisms
auts = []
for perm in permutations(range(n)):
    ok = True
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j] != A[perm[i]][perm[j]]:
                ok = False; break
        if not ok: break
    if ok:
        auts.append(tuple(perm))

print("  |Aut| = %d" % len(auts))
print("  Automorphisms:")
for aut in auts:
    cycles = []
    visited = [False]*n
    for i in range(n):
        if visited[i]: continue
        cycle = []
        j = i
        while not visited[j]:
            visited[j] = True
            cycle.append(j)
            j = aut[j]
        if len(cycle) > 1:
            cycles.append(tuple(cycle))
    print("    %s → cycles %s" % (aut, cycles if cycles else "identity"))

# Group HPs by Aut orbit
hp_set = set(hps)
hp_orbits = []
remaining = set(hps)
while remaining:
    hp = next(iter(remaining))
    orbit = set()
    for aut in auts:
        mapped = tuple(aut[v] for v in hp)
        orbit.add(mapped)
    hp_orbits.append(sorted(orbit))
    remaining -= orbit

print("\n  HP orbits (each = one tiling):")
for idx, orbit in enumerate(hp_orbits):
    print("    Orbit %d (size %d):" % (idx, len(orbit)))
    for hp in orbit[:3]:
        print("      HP %s" % (hp,))

    # What tiling does this orbit correspond to?
    hp0 = orbit[0]
    sigma = {}
    for i, v in enumerate(hp0):
        sigma[v] = n - 1 - i
    A_rel = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                A_rel[sigma[i]][sigma[j]] = A[i][j]
    bits = 0
    for k, (hi, lo) in enumerate(tile_pairs):
        if A_rel[lo][hi]:
            bits |= (1 << k)

    # Neighbor profile
    profile = []
    for k in range(m):
        mj = mid_map[tc[bits ^ (1 << k)]]
        profile.append(mj)
    class_counts = Counter(profile)
    n_neutral = class_counts.get(mi, 0)

    print("    → tiling %s, profile: %s, neutral=%d" %
          (format(bits, '06b'), dict(class_counts), n_neutral))

# ============================================================
banner("THE PATTERN: COMPLEMENT PAIRING IN TILINGS")
# ============================================================

# Let me check: for paired tilings (same class-count profile),
# is one the COMPLEMENT of the other (flip ALL tiles)?
# The complement in tiling space: bits → bits XOR (2^m - 1).

print("  For each pair: is one the tiling-complement of the other?\n")

complement_count = 0; non_complement_count = 0
relabeling_count = 0

for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid0 = merged[mi][0]
    H = info[cid0]['H']
    all_tilings = []
    for c in merged[mi]:
        all_tilings.extend(members[c])
    fiber = len(all_tilings)
    if fiber <= 2: continue

    tiling_profiles = {}
    for bits in all_tilings:
        profile = Counter()
        for k in range(m):
            mj = mid_map[tc[bits ^ (1 << k)]]
            profile[mj] += 1
        tiling_profiles[bits] = tuple(sorted(profile.items()))

    profile_groups = defaultdict(list)
    for bits, prof in tiling_profiles.items():
        profile_groups[prof].append(bits)

    for prof, group in profile_groups.items():
        if len(group) == 2:
            b1, b2 = group
            comp_b1 = b1 ^ ((1 << m) - 1)
            if comp_b1 == b2:
                complement_count += 1
            else:
                non_complement_count += 1
                # Check if they're related by a relabeling
                A1 = [[0]*n for _ in range(n)]
                for i in range(1, n): A1[i][i-1] = 1
                for k, (hi, lo) in enumerate(tile_pairs):
                    if (b1 >> k) & 1: A1[lo][hi] = 1
                    else: A1[hi][lo] = 1
                A2 = [[0]*n for _ in range(n)]
                for i in range(1, n): A2[i][i-1] = 1
                for k, (hi, lo) in enumerate(tile_pairs):
                    if (b2 >> k) & 1: A2[lo][hi] = 1
                    else: A2[hi][lo] = 1
                for perm in permutations(range(n)):
                    ok = True
                    for i in range(n):
                        for j in range(n):
                            if i != j and A1[i][j] != A2[perm[i]][perm[j]]:
                                ok = False; break
                        if not ok: break
                    if ok:
                        relabeling_count += 1
                        break

print("  Complement pairs: %d" % complement_count)
print("  Non-complement pairs: %d (of which %d related by relabeling)" %
      (non_complement_count, relabeling_count))
print("  Total pairs: %d" % (complement_count + non_complement_count))

# ============================================================
banner("WHAT IS THE PAIRING? DETAILED EXAMINATION")
# ============================================================

# For each non-complement pair, find the actual relationship
for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
    cid0 = merged[mi][0]
    H = info[cid0]['H']
    all_tilings = []
    for c in merged[mi]:
        all_tilings.extend(members[c])
    fiber = len(all_tilings)
    if fiber <= 2: continue

    tiling_profiles = {}
    for bits in all_tilings:
        profile = Counter()
        for k in range(m):
            mj = mid_map[tc[bits ^ (1 << k)]]
            profile[mj] += 1
        tiling_profiles[bits] = tuple(sorted(profile.items()))

    profile_groups = defaultdict(list)
    for bits, prof in tiling_profiles.items():
        profile_groups[prof].append(bits)

    pairs = [(prof, grp) for prof, grp in profile_groups.items() if len(grp) == 2]
    if not pairs: continue

    print("\n  mid=%d (H=%d, fiber=%d):" % (mi, H, fiber))

    for prof, (b1, b2) in pairs[:3]:
        # Build tournaments
        A1 = [[0]*n for _ in range(n)]
        for i in range(1, n): A1[i][i-1] = 1
        for k, (hi, lo) in enumerate(tile_pairs):
            if (b1 >> k) & 1: A1[lo][hi] = 1
            else: A1[hi][lo] = 1

        A2 = [[0]*n for _ in range(n)]
        for i in range(1, n): A2[i][i-1] = 1
        for k, (hi, lo) in enumerate(tile_pairs):
            if (b2 >> k) & 1: A2[lo][hi] = 1
            else: A2[hi][lo] = 1

        # Find ALL permutations mapping A1 to A2
        relabelings = []
        for perm in permutations(range(n)):
            ok = True
            for i in range(n):
                for j in range(n):
                    if i != j and A1[i][j] != A2[perm[i]][perm[j]]:
                        ok = False; break
                if not ok: break
            if ok:
                relabelings.append(perm)

        # What do these relabelings do to the base path?
        bp = list(range(n-1, -1, -1))
        for sigma in relabelings[:3]:
            mapped_bp = [sigma[v] for v in bp]
            # Analyze the permutation σ
            cycles = []
            visited = [False]*n
            for i in range(n):
                if visited[i]: continue
                cycle = []
                j = i
                while not visited[j]:
                    visited[j] = True
                    cycle.append(j)
                    j = sigma[j]
                if len(cycle) > 1:
                    cycles.append(tuple(cycle))

            print("    σ=%s, cycles=%s, maps BP=%s→%s" %
                  (sigma, cycles if cycles else "id", bp, mapped_bp))

            # Does σ reverse the base path?
            is_reversal = (mapped_bp == list(range(n)))
            print("      BP reversed: %s" % is_reversal)

            # Is σ an involution (σ² = id)?
            sigma2 = tuple(sigma[sigma[i]] for i in range(n))
            is_involution = (sigma2 == tuple(range(n)))
            print("      Involution: %s" % is_involution)

            # Hamming distance in tiling space
            hamming = bin(b1 ^ b2).count('1')
            print("      Tiling hamming distance: %d / %d" % (hamming, m))

print("\nDone. opus-2026-03-24-S293")
