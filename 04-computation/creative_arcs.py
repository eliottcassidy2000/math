"""
creative_arcs.py — opus-2026-04-04-S15

THE USER'S KEY INSIGHT: From the transitive, only the source-to-sink arc
(the longest-range arc, skip = n-1) creates the cyclic at n=3.

Extend: for EACH arc (tile) flipped from the transitive, what is the
resulting tournament's "complexity"? Is there a universal law relating
arc range (skip) to the creativity of the flip?

Also: from the H-MAXIMIZER, which arcs are most destructive?
And: which arcs change the iso class most vs least?
"""
import numpy as np
from itertools import permutations
from collections import defaultdict, Counter
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, make_tiles, tiling_to_tournament, count_hamiltonian_paths

def tournament_canonical(T, n):
    best = None
    for perm in permutations(range(n)):
        key = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or key < best:
            best = key
    return best

def analyze_creative_arcs(n):
    """For the transitive tournament, analyze how each arc flip creates structure."""
    print(f"\n{'='*60}")
    print(f"CREATIVE ARCS AT n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)

    # The transitive tournament = all tiles forward (tiling = 0)
    T_trans = tiling_to_tournament(n, tiles, 0)

    print(f"\n  FROM TRANSITIVE (H=1):")
    print(f"  Each tile flip creates a new tournament. How much structure?")
    print(f"  {'tile':>8s} {'skip':>5s} {'H_new':>6s} {'ΔH':>5s} {'c3':>4s} "
          f"{'effect':>10s}")

    for k in range(m):
        tile = tiles[k]
        skip = tile[0] - tile[1]
        tiling_flip = 1 << k
        T_new = tiling_to_tournament(n, tiles, tiling_flip)
        h_new = count_hamiltonian_paths(T_new, n)
        delta_h = h_new - 1

        # Count 3-cycles created
        c3 = 0
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if T_new[a][b]+T_new[b][c]+T_new[c][a]==3 or \
                       T_new[a][c]+T_new[c][b]+T_new[b][a]==3:
                        c3 += 1

        # The linear coefficient is 2^{skip-1}
        expected = 2**(skip-1)
        print(f"  {str(tile):>8s} {skip:5d} {h_new:6d} {delta_h:+5d} {c3:4d} "
              f"{'c_k='+str(expected):>10s}")

    # From the H-MAXIMIZER
    H_all, _, _ = compute_all_H(n)
    max_h = max(H_all)
    max_tilings = [t for t in range(1 << m) if H_all[t] == max_h]

    print(f"\n  FROM H-MAXIMIZER (H={max_h}), {len(max_tilings)} tilings:")
    if max_tilings:
        t_max = max_tilings[0]
        for k in range(m):
            tile = tiles[k]
            skip = tile[0] - tile[1]
            t_flip = t_max ^ (1 << k)
            h_flip = int(H_all[t_flip])
            delta_h = h_flip - max_h
            was_backward = (t_max >> k) & 1
            direction = "←" if was_backward else "→"

            print(f"  {str(tile):>8s} skip={skip} {direction}: H={h_flip:4d} "
                  f"(ΔH={delta_h:+5d})")

    # THE KEY TABLE: for each skip value, what's the AVERAGE effect?
    print(f"\n  AVERAGE EFFECT BY SKIP:")
    by_skip = defaultdict(list)
    for k in range(m):
        skip = tiles[k][0] - tiles[k][1]
        h_from_trans = int(H_all[1 << k])
        by_skip[skip].append(h_from_trans - 1)

    for skip in sorted(by_skip.keys()):
        vals = by_skip[skip]
        print(f"    skip={skip}: mean ΔH from transitive = {np.mean(vals):.1f}, "
              f"values = {vals}")

def analyze_class_change_probability(n):
    """For each tile, what fraction of flips change the iso class?"""
    print(f"\n{'='*60}")
    print(f"CLASS CHANGE PROBABILITY BY TILE n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m
    H_all, _, _ = compute_all_H(n)

    # Compute iso class for each tiling
    tiling_class = {}
    for t in range(N):
        T = tiling_to_tournament(n, tiles, t)
        tiling_class[t] = tournament_canonical(T, n)

    # For each tile k: what fraction of flips of that tile change the class?
    print(f"  {'tile':>8s} {'skip':>5s} {'changes':>8s} {'total':>6s} {'frac':>8s} "
          f"{'mean_|ΔH|':>10s}")

    for k in range(m):
        tile = tiles[k]
        skip = tile[0] - tile[1]
        changes = 0
        total = 0
        delta_h_vals = []

        for t in range(N):
            t_flip = t ^ (1 << k)
            total += 1
            h_diff = abs(int(H_all[t_flip]) - int(H_all[t]))
            delta_h_vals.append(h_diff)
            if tiling_class[t] != tiling_class[t_flip]:
                changes += 1

        frac = changes / total
        mean_abs_dh = np.mean(delta_h_vals)
        print(f"  {str(tile):>8s} {skip:5d} {changes:8d} {total:6d} {frac:8.4f} "
              f"{mean_abs_dh:10.2f}")

    # The SKIP DEPENDENCE
    print(f"\n  SUMMARY BY SKIP:")
    by_skip_change = defaultdict(list)
    by_skip_dh = defaultdict(list)

    for k in range(m):
        skip = tiles[k][0] - tiles[k][1]
        changes = sum(1 for t in range(N)
                     if tiling_class[t] != tiling_class[t ^ (1 << k)])
        frac = changes / N
        mean_dh = np.mean([abs(int(H_all[t ^ (1<<k)]) - int(H_all[t]))
                          for t in range(N)])
        by_skip_change[skip].append(frac)
        by_skip_dh[skip].append(mean_dh)

    for skip in sorted(by_skip_change.keys()):
        fracs = by_skip_change[skip]
        dhs = by_skip_dh[skip]
        print(f"    skip={skip}: mean class-change frac = {np.mean(fracs):.4f}, "
              f"mean |ΔH| = {np.mean(dhs):.2f}")

def the_abstract_principle(n):
    """
    THE DEEP PATTERN from n=3:
    From order (transitive), the longest arc creates the most disorder.
    From disorder (H-max), every arc destroys structure.

    At n=3: 1 tile, 1 flip. From transitive: skip-2 arc creates the 3-cycle.
    At n=4: 3 tiles. From transitive: skip-3 → H=5 (max), skip-2 → H=3.
    At n=5: 6 tiles. From transitive: skip-4 → H=9, skip-3 → H=5, skip-2 → H=3.
    At n=6: 10 tiles. From transitive: skip-5 → H=17, skip-4 → H=9, etc.

    PATTERN: flipping skip-s tile from transitive gives H = 1 + 2^{s-1}.
    (This IS THM-284: c_k = 2^{skip-1}.)

    The user's observation at n=3 is the SIMPLEST CASE of THM-284:
    flipping the only tile (skip 2) adds 2^1 = 2 HPs, giving H = 3.

    ABSTRACT EXTENSION:
    "The transitive tournament is the ground state. Each tile is an excitation.
    Skip-s tiles carry energy 2^{s-1}. The first excitation is the apex (highest
    energy). The ground state has H=1, the first excited state has H=1+2^{n-2}."

    "From the H-maximizer (the 'fully excited' state), every de-excitation
    (flip) reduces H. The H-maximizer is the unique maximum-energy state."

    This is the ANTIFERROMAGNETIC PICTURE from THM-290:
    - The field h_k = 2^{skip-1} is the excitation energy
    - The coupling J_{ij} < 0 means excitations COMPETE
    - The H-maximizer is the state that resolves competition optimally
    """
    print(f"\n{'='*60}")
    print(f"THE ABSTRACT PRINCIPLE")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H_all, _, _ = compute_all_H(n)

    print(f"""
  THE TOURNAMENT AS A QUANTUM SYSTEM:

  Ground state:  Transitive (H=1, all arcs "forward")
  Excitations:   Each backward tile is an excitation with energy 2^(skip-1)
  Interactions:  Same-end tiles COMPETE (antiferromagnetic coupling)
                 Cross-end/disjoint tiles COOPERATE (but multiplicatively weaker)

  FROM GROUND STATE (transitive):
    Flipping tile k adds energy c_k = 2^(skip-1) = the LINEAR effect
    H(transitive + e_k) = 1 + 2^(skip-1)

    The USER'S n=3 OBSERVATION: the ONLY excitation (skip-2) creates
    the 3-cycle (the fundamental curvature quantum). H goes 1→3.

    At larger n: the APEX excitation (skip = n-1) creates the LARGEST
    jump: H → 1 + 2^(n-2). This is the "hardest" excitation.

  FROM EXCITED STATE (H-maximizer):
    Every flip REDUCES H (the maximum is a peak).
    The H-maximizer is FRAGILE in the sense that MANY flips reduce H.
    But it's NOT fragile in the iso-class sense — it may have self-loops
    (flips that stay in the same class but a different tournament).
    """)

    # Verify: from transitive, H(e_k) = 1 + 2^{skip-1}
    print(f"  Verification at n={n}:")
    for k in range(m):
        tile = tiles[k]
        skip = tile[0] - tile[1]
        h = int(H_all[1 << k])
        expected = 1 + 2**(skip - 1)
        match = "✓" if h == expected else "✗"
        print(f"    tile {tile} skip={skip}: H = {h}, 1+2^(skip-1) = {expected} {match}")

    # From H-maximizer: what's the mean ΔH per skip?
    max_h = max(H_all)
    max_tilings = [t for t in range(1 << m) if H_all[t] == max_h]
    if max_tilings:
        t_max = max_tilings[0]
        print(f"\n  From H-maximizer (H={max_h}):")
        by_skip = defaultdict(list)
        for k in range(m):
            skip = tiles[k][0] - tiles[k][1]
            h_flip = int(H_all[t_max ^ (1 << k)])
            by_skip[skip].append(h_flip - max_h)

        for skip in sorted(by_skip.keys()):
            vals = by_skip[skip]
            print(f"    skip={skip}: mean ΔH = {np.mean(vals):+.1f}, "
                  f"values = {vals[:6]}")

if __name__ == "__main__":
    for n in [3, 4, 5, 6]:
        analyze_creative_arcs(n)

    for n in [5, 6]:
        analyze_class_change_probability(n)

    for n in [5, 6]:
        the_abstract_principle(n)
