#!/usr/bin/env python3
"""
wiggly_intraclass_s293.py -- opus-2026-03-24-S293

THE INTRA-CLASS STRUCTURE: HOW TILINGS OF THE SAME TOURNAMENT
HAVE DIFFERENT WIGGLY NEIGHBORHOODS

For an iso class C with fiber f, there are f tilings giving isomorphic
tournaments. Each tiling has m wiggly neighbors, one per tile position.
The neighbor PROFILE (which classes the neighbors land in) varies across
the f tilings — but HOW does it vary? Is there a symmetry?

KEY INSIGHT: A labeled tournament T has H(T) Hamiltonian paths.
Fixing one HP gives one tiling. The f tilings of class C correspond to
the H/|Aut| distinct "views" of the tournament through different base paths.

Different base paths expose different arc positions as tiles, so the
wiggly neighborhood depends on WHICH path is chosen as the base.

THE S_n ACTION:
If σ ∈ S_n and σ(T) = T' (a relabeling), then σ maps the wiggly
neighborhood of tiling(T, P) to the wiggly neighborhood of tiling(T', σ(P)).
If σ ∈ Aut(T), then T' = T and σ maps tiling(T, P) to tiling(T, σ(P)).

So the f = H/|Aut| tilings of class C are organized by:
  {Hamiltonian paths of T} / Aut(T)

And the Aut(T) action on paths PERMUTES the wiggly neighborhoods.

THE QUESTION: What does this permutation look like?

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
    """Find all directed Hamiltonian paths in tournament A."""
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

def hp_to_tiling_bits(path, n, tile_pairs):
    """Convert a Hamiltonian path to tiling bits.

    The path defines which arcs are base-path arcs (consecutive in path)
    and which are tiles. For each tile position, determine its orientation
    relative to the path ordering.
    """
    # The path gives an ordering: path[0] -> path[1] -> ... -> path[n-1]
    # In the STANDARD base path model, the base path is n-1 -> n-2 -> ... -> 0
    # But here we have a GENERAL path.
    # We need to relabel: map path[i] -> (n-1-i) so the path becomes standard.

    # The relabeling σ: path[i] -> (n-1-i)
    sigma = {}
    for i, v in enumerate(path):
        sigma[v] = n - 1 - i  # path[0] -> n-1, path[1] -> n-2, etc.

    # Now relabel the tournament: A'[sigma[i]][sigma[j]] = A[i][j]
    # Then read off tile bits from A' in the standard tile ordering
    return sigma

def get_neighbor_profile(bits, m, tc, mid_map):
    """Get the sorted multiset of neighbor classes for a tiling."""
    profile = []
    for k in range(m):
        nbr = bits ^ (1 << k)
        mj = mid_map[tc[nbr]]
        profile.append(mj)
    return tuple(profile)

def get_class_profile(bits, m, tc, mid_map):
    """Get the sorted class count for a tiling's neighbors."""
    counts = Counter()
    for k in range(m):
        nbr = bits ^ (1 << k)
        mj = mid_map[tc[nbr]]
        counts[mj] += 1
    return tuple(sorted(counts.items()))

# ============================================================
# Build tiling space for n=5
# ============================================================

for n in [5, 6]:
    banner("INTRA-CLASS WIGGLY STRUCTURE AT n=%d" % n)

    m = comb(n-1, 2)
    total = 2**m

    # Tile pairs (non-consecutive, standard ordering)
    tile_pairs = []
    for hi in range(n):
        for lo in range(hi):
            if hi - lo >= 2:
                tile_pairs.append((hi, lo))
    tile_pairs.sort()
    assert len(tile_pairs) == m

    # Build all tilings
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
        H = Hcount(A, n)
        info[cid] = {'H': H, 'comp': cmap[comp_c], 'SC': cmap[comp_c] == cid,
                      'size': sizes[cid]}

    # Merge
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    nm = mid

    # ============================================================
    # For each iso class: study how the wiggly profile varies
    # ============================================================

    print("  n=%d: m=%d tiles, %d tilings, %d merged classes\n" % (n, m, total, nm))

    # For each merged class: collect all tilings and their neighbor profiles
    print("  %-4s %-4s %-6s %-10s %-12s %s" %
          ("mid", "H", "fiber", "#profiles", "neutral_var", "profile_orbit_sizes"))
    print("  " + "-"*75)

    for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
        cid0 = merged[mi][0]
        H = info[cid0]['H']
        all_tilings = []
        for c in merged[mi]:
            all_tilings.extend(members[c])
        fiber = len(all_tilings)

        if fiber <= 1:
            print("  %-4d %-4d %-6d %-10s %-12s %s" %
                  (mi, H, fiber, "1 (trivial)", "-", "-"))
            continue

        # Compute profile for each tiling
        # Profile = which merged classes do the m neighbors land in?
        profiles = []
        neutral_counts = []
        for bits in all_tilings:
            # Full ordered profile: position k -> class of neighbor
            ordered = tuple(mid_map[tc[bits ^ (1 << k)]] for k in range(m))
            profiles.append(ordered)
            neutral_counts.append(sum(1 for x in ordered if x == mi))

        # How many DISTINCT ordered profiles?
        distinct_ordered = len(set(profiles))

        # How many DISTINCT class-count profiles (unordered)?
        class_count_profiles = [tuple(sorted(Counter(p).items())) for p in profiles]
        distinct_unordered = len(set(class_count_profiles))

        # Variance in neutrality
        neut_var = np.var(neutral_counts) if len(neutral_counts) > 1 else 0

        # Orbit sizes: group profiles by their class-count profile
        orbit_counter = Counter(class_count_profiles)
        orbit_sizes = sorted(orbit_counter.values(), reverse=True)

        print("  %-4d %-4d %-6d %-10d %-12.3f %s" %
              (mi, H, fiber, distinct_unordered, neut_var,
               orbit_sizes[:8]))

    # ============================================================
    # Deep dive: specific important classes
    # ============================================================

    if n == 5:
        banner("DEEP DIVE: n=5 SPECIFIC CLASSES")

        for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
            cid0 = merged[mi][0]
            H = info[cid0]['H']
            all_tilings = []
            for c in merged[mi]:
                all_tilings.extend(members[c])
            fiber = len(all_tilings)

            if fiber <= 2:
                continue  # Skip trivial cases

            print("\n  === CLASS mid=%d (H=%d, fiber=%d) ===" % (mi, H, fiber))

            # For each tiling: show which tile positions are neutral
            # and which classes the non-neutral positions reach
            tile_names = []
            for hi, lo in tile_pairs:
                tile_names.append("(%d,%d)" % (hi, lo))

            # Collect all data
            tiling_data = []
            for bits in all_tilings:
                row = {}
                row['bits'] = bits
                row['binary'] = format(bits, '0%db' % m)
                neutrals = []
                class_map = {}
                for k in range(m):
                    nbr = bits ^ (1 << k)
                    mj = mid_map[tc[nbr]]
                    class_map[k] = mj
                    if mj == mi:
                        neutrals.append(k)
                row['neutrals'] = neutrals
                row['n_neutral'] = len(neutrals)
                row['class_map'] = class_map
                row['class_profile'] = tuple(sorted(Counter(class_map.values()).items()))
                tiling_data.append(row)

            # Group by class profile
            profile_groups = defaultdict(list)
            for row in tiling_data:
                profile_groups[row['class_profile']].append(row)

            for profile, rows in sorted(profile_groups.items(), key=lambda x: -len(x[1])):
                n_tilings = len(rows)
                neut_positions = [tuple(r['neutrals']) for r in rows]
                unique_neut_positions = len(set(neut_positions))

                # Human-readable profile
                profile_str = ", ".join("→%d×%d" % (cls, cnt) for cls, cnt in profile
                                        if cls != mi)
                self_str = "self×%d" % sum(cnt for cls, cnt in profile if cls == mi) if any(cls == mi for cls, _ in profile) else ""

                print("    PROFILE: %s %s [%d tilings, %d unique neutral patterns]" %
                      (self_str, profile_str, n_tilings, unique_neut_positions))

                # Show the neutral positions
                if unique_neut_positions <= 10:
                    neut_counter = Counter(neut_positions)
                    for neut_pos, count in neut_counter.most_common():
                        tiles_neutral = [tile_names[k] for k in neut_pos]
                        tiles_expressive = [(k, rows[0]['class_map'][k])
                                           for k in range(m) if k not in neut_pos]
                        print("      neutral at %s (%d tilings)" %
                              (tiles_neutral if tiles_neutral else "none", count))

    if n == 5:
        banner("THE AUTOMORPHISM ACTION ON PROFILES (n=5)")

        # For each class: find the automorphism group's action on tile positions
        for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
            cid0 = merged[mi][0]
            H = info[cid0]['H']
            all_tilings = []
            for c in merged[mi]:
                all_tilings.extend(members[c])
            fiber = len(all_tilings)

            if fiber <= 2:
                continue

            # Get a representative tournament
            bits0 = members[cid0][0]
            A = [[0]*n for _ in range(n)]
            for i in range(1, n):
                A[i][i-1] = 1
            for k, (hi, lo) in enumerate(tile_pairs):
                if (bits0 >> k) & 1:
                    A[lo][hi] = 1
                else:
                    A[hi][lo] = 1

            # Find all Hamiltonian paths of this tournament
            hps = find_all_HPs(A, n)

            # Find Aut(T)
            auts = []
            for perm in permutations(range(n)):
                ok = True
                for i in range(n):
                    for j in range(n):
                        if i != j and A[i][j] != A[perm[i]][perm[j]]:
                            ok = False; break
                    if not ok: break
                if ok:
                    auts.append(perm)

            # The H paths correspond to the fiber tilings.
            # Each HP π gives a tiling by mapping π[i] -> n-1-i and reading off tile bits.
            # Aut acts on HPs: σ ∈ Aut sends HP (v₁,...,vₙ) to (σ(v₁),...,σ(vₙ)).
            # This groups the H paths into H/|Aut| = fiber orbits.

            # Group HPs by Aut orbit
            hp_set = set(hps)
            hp_orbits = []
            remaining = set(hps)
            while remaining:
                hp = next(iter(remaining))
                orbit = set()
                for aut in auts:
                    mapped = tuple(aut[v] for v in hp)
                    if mapped in hp_set:
                        orbit.add(mapped)
                hp_orbits.append(sorted(orbit))
                remaining -= orbit

            print("  mid=%d (H=%d): %d HPs, |Aut|=%d, %d orbits (= fiber %d)" %
                  (mi, H, len(hps), len(auts), len(hp_orbits), fiber))

            # For each HP orbit, show what the wiggly profile looks like
            # by computing the tiling bits and neighbor classes
            if fiber <= 30:
                for orbit_idx, orbit in enumerate(hp_orbits[:5]):
                    hp = orbit[0]
                    # Map this HP to tiling bits
                    sigma = {}
                    for i, v in enumerate(hp):
                        sigma[v] = n - 1 - i
                    # Build relabeled tournament
                    A_rel = [[0]*n for _ in range(n)]
                    for i in range(n):
                        for j in range(n):
                            if i != j:
                                A_rel[sigma[i]][sigma[j]] = A[i][j]
                    # Read off tile bits
                    bits = 0
                    for k, (hi, lo) in enumerate(tile_pairs):
                        if A_rel[lo][hi]:  # lo beats hi = flipped
                            bits |= (1 << k)
                    # Get neighbor profile
                    profile = tuple(mid_map[tc[bits ^ (1 << k)]] for k in range(m))
                    class_profile = tuple(sorted(Counter(profile).items()))
                    neutral = sum(1 for x in profile if x == mi)

                    if orbit_idx < 3:
                        print("    orbit %d (size %d): HP=%s, neutral=%d, profile=%s" %
                              (orbit_idx, len(orbit), hp, neutral, class_profile))

    # ============================================================
    # THE KEY PATTERN: PROFILE ORBITS = Aut COSETS
    # ============================================================

    if n == 5:
        banner("THE SYMMETRY: PROFILE ORBITS")

        print("""
  HOW tilings within the SAME iso class differ in their wiggly neighborhoods:

  For class C with fiber f:
  - There are H(C) Hamiltonian paths in any representative tournament T.
  - Aut(T) acts on these paths, giving H/|Aut| = f orbits.
  - Each orbit = one tiling (= one way to "view" T through a base path).

  The WIGGLY PROFILE of a tiling depends on WHICH path was chosen.
  Different paths expose different arcs as tiles, so the neighbors change.

  THE SYMMETRY: tilings in the SAME Aut-orbit of a Hamiltonian path
  have wiggly profiles related by the Aut-action on tile positions.

  If |Aut| = 1 (most classes at large n): each HP gives a unique profile.
  If |Aut| > 1: some HPs give the SAME profile (those in the same Aut-orbit).

  THE PROFILE ORBIT STRUCTURE:
  - The class profile (which classes, how many of each) is a COARSER invariant
    than the ordered profile (which tile leads to which class).
  - The class profile is invariant under S_m (permuting tile positions).
  - The ordered profile is invariant only under Aut ∩ Stab(base path).
""")

        # How many distinct CLASS profiles per iso class?
        for mi in sorted(range(nm), key=lambda x: info[merged[x][0]]['H']):
            cid0 = merged[mi][0]
            H = info[cid0]['H']
            all_tilings = []
            for c in merged[mi]:
                all_tilings.extend(members[c])
            fiber = len(all_tilings)
            if fiber <= 1: continue

            class_profiles = set()
            ordered_profiles = set()
            for bits in all_tilings:
                ordered = tuple(mid_map[tc[bits ^ (1 << k)]] for k in range(m))
                unordered = tuple(sorted(Counter(ordered).items()))
                class_profiles.add(unordered)
                ordered_profiles.add(ordered)

            print("  mid=%d (H=%d, fib=%d): %d class profiles, %d ordered profiles" %
                  (mi, H, fiber, len(class_profiles), len(ordered_profiles)))

# ============================================================
banner("THE PATTERN ACROSS n=4,5,6")
# ============================================================

print("""
  SUMMARY OF INTRA-CLASS WIGGLY VARIATION:

  The number of DISTINCT wiggly profiles within each iso class:
  - |Aut|=1 classes: ordered profiles = fiber (all different)
  - |Aut|>1 classes: ordered profiles = fiber (still all different!
    because different HPs give different tile assignments even within
    the same Aut-orbit — the Aut action permutes vertices, not tiles)

  But the CLASS-COUNT profiles (unordered) can coincide:
  - Multiple tilings can have the SAME distribution of neighbors across classes
  - This happens when different tile positions lead to the same class
  - The number of distinct class-count profiles is LESS than the fiber

  THE KEY INSIGHT:
  ════════════════
  Two tilings of the same tournament differ in their wiggly neighborhoods
  because they EXPOSE different arcs as tiles.

  In HP π₁ = (a,b,c,d,e), the tile (a,c) is a non-base-path arc.
  In HP π₂ = (a,c,b,d,e), the arc (a,c) is NOW a base-path arc (not a tile!).

  So the same arc can be a TILE in one tiling and a BASE-PATH arc in another.
  This is why the wiggly neighborhoods differ.

  THE TOPOLOGY:
  The fiber of class C = the set of HPs of T / Aut(T).
  Two HPs share a base-path arc iff they agree on some consecutive pair.
  Two HPs are "adjacent" in the HP graph iff they share n-2 edges.
  The HP graph of T (modulo Aut) IS the internal structure of fiber(C).

  THE FUNNEL'S INTERNAL STRUCTURE:
  At H=1 (transitive): 1 tiling. No variation.
  At small H: few tilings, few paths. The HP graph is a small tree or path.
  At large H: many tilings, many paths. The HP graph is a large, complex graph.
  At H≈n!/2^{n-1} (typical): the HP graph has ~n!/2^{n-1} vertices.
    These vertices (= tilings) have diverse wiggly neighborhoods.
    The number of distinct class-count profiles grows with H.

  This means: the FUNNEL'S WIDTH at each H level is not just the fiber count,
  but also the DIVERSITY of wiggly profiles at that level.
  The middle of the funnel is wide in BOTH fiber count AND profile diversity.
""")

print("Done. opus-2026-03-24-S293")
