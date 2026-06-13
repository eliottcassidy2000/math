#!/usr/bin/env python3
"""The apex tile: the one arc that closes the polygon.

opus-2026-06-01-S529

The outside of the polygon = base path (n-1 edges) + apex tile (1 edge).
The apex tile (n,1) connects SOURCE to SINK of the base path.

In the tiling model: tile A = (n,1) = the first tile in the explorer.
It's the ONLY tile on the polygon boundary. All other tiles are interior.

ALIGNED apex (source → sink):  the ranking stays open, transitive tendency
ANTI-ALIGNED apex (sink → source): the ranking CLOSES into a Hamiltonian CYCLE

For LRC: if the observer is at the top of the ranking (source vertex),
the apex connects the observer to the FARTHEST runner. The lonely
condition requires this tile to be aligned (observer beats the farthest).

KEY QUESTIONS:
1. How often is the apex aligned vs anti-aligned along the LRC walk?
2. When the apex flips, what happens to the tournament's H and cycle structure?
3. How does the apex interact with the cascade proof?
4. For n=14: the apex connects runner 13 to runner 1 — what's its role?
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, comb, factorial
from functools import reduce
from collections import Counter
from itertools import permutations


ONE = Fraction(1)
ZERO = Fraction(0)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    f = frac(x)
    return min(f, ONE - f)


# ═══════════════════════════════════════════════════════════════
# PART 1: The apex tile in the tiling model
# ═══════════════════════════════════════════════════════════════

def apex_tile_basics(n_values=[4, 5, 6, 7, 8, 14]):
    """The apex tile (n,1) in the staircase."""
    print("=" * 70)
    print("PART 1: The apex tile — closing the polygon")
    print("=" * 70)
    print()

    for n in n_values:
        m = comb(n - 1, 2)  # number of tiles
        print(f"n={n}:")
        print(f"  base path: {n} → {n-1} → ... → 1  ({n-1} edges)")
        print(f"  apex tile: ({n}, 1) = tile A")
        print(f"  total tiles: C({n-1},2) = {m}")
        print(f"  interior tiles (non-apex): {m - 1}")
        print(f"  polygon boundary = base path + apex = {n} edges")
        print()

        # The apex's position in the staircase:
        # The staircase is a right triangle with rows y=1,...,n-2
        # and columns x=y+2,...,n.
        # Tile (n,1) is at row y=1, column x=n — the TOP-LEFT corner.
        # It's the LONGEST skip (distance n-1 in the base path).

        print(f"  apex position: row 1 (bottom), column {n} (rightmost)")
        print(f"  skip distance: {n-1} (longest possible)")
        print(f"  in the staircase, the apex is at the HYPOTENUSE vertex")
        print()

        # Aligned vs anti-aligned:
        # Aligned (n→1): the base-path source beats the sink.
        # The tournament is "more transitive" — source dominates.
        # Anti-aligned (1→n): sink beats source. A Hamiltonian cycle
        # n→n-1→...→1→n exists if all base-path arcs are consistent.

        print(f"  ALIGNED: {n}→1 (source beats sink, no boundary cycle)")
        print(f"  ANTI-ALIGNED: 1→{n} (sink beats source, Hamiltonian cycle)")
        print()

    print("THE APEX IN THE 'EVERYTHING IS THE TRIANGLE' FRAMEWORK:")
    print("  The staircase is a right triangle. The three sides are:")
    print("    - VERTICAL (source column): scores, hierarchy, cut space")
    print("    - HORIZONTAL (sink row): complement, anti-hierarchy")
    print("    - HYPOTENUSE (anti-diagonal): the fiber fraction, Walsh order-2")
    print("  The APEX TILE sits at the VERTEX where the hypotenuse meets")
    print("  the horizontal leg — the corner of the triangle.")
    print("  It connects the vertical (source) to the horizontal (sink).")
    print("  Flipping the apex rotates between transitive and cyclic global structure.")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 2: The apex along the LRC walk
# ═══════════════════════════════════════════════════════════════

def apex_on_lrc_walk(n_values=[4, 5, 6, 7]):
    """Track the apex tile's orientation along the LRC walk.

    For LRC with observer = vertex n (base path source):
    Apex = arc between observer and runner 1 (base path sink).
    Apex aligned (observer→runner1) iff ||v_1 t|| ≥ 1/n.
    (Runner 1 = the SLOWEST runner, speed v_1.)

    Wait — in the tiling model, vertex n is the source and vertex 1 is
    the sink. For LRC, the observer's position in the ranking depends
    on the speed ordering. Let me use the STANDARD LRC labeling:
    observer = vertex 0, runners = vertices 1,...,n-1 (by speed).

    Then the base path for the runner sub-tournament is:
    (n-1) → (n-2) → ... → 1 (runners in decreasing speed order).
    The apex tile of the RUNNER sub-tournament is (n-1, 1):
    speed n-1 (fastest) vs speed 1 (slowest).

    The OBSERVER is a separate vertex. The "full" tournament on all n
    vertices has a different base path. But for the runner sub-tournament
    (which determines the metagraph class), the apex is (n-1, 1).

    For LRC: the observer's arcs are separate from the runner tiles.
    The runner apex (n-1, 1) connects the fastest to the slowest runner.
    It has nothing directly to do with the observer's loneliness.

    BUT: the FULL n-vertex tournament's apex (if we include the observer
    as vertex n in the base path n→n-1→...→1→0) is (n, 0) = the arc
    between the observer and the SLOWEST runner.

    Actually, let me reconsider. If the base path is:
    observer(n) → runner_{n-1} → runner_{n-2} → ... → runner_1
    then the apex tile is (observer, runner_1) = (n, 1).

    The apex is aligned iff observer→runner_1, which is the LRC condition
    for the SLOWEST runner to be safe: ||v_1 t|| ≥ 1/n.

    So the apex tile's status = "is the slowest runner safe?"
    """
    print("=" * 70)
    print("PART 2: The apex tile on the LRC walk")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))  # runner speeds 1,...,n-1
        thr = Fraction(1, n)

        # In the full tournament with observer as vertex n (source):
        # base path: observer → runner_{n-1} → ... → runner_1
        # apex tile: (observer, runner_1) = is runner_1 (speed 1) safe?

        # Track apex status over the walk
        walls = set()
        walls.add(ZERO)
        for v in speeds:
            for a in [1, n - 1]:
                for k in range(v):
                    walls.add(Fraction(k * n + a, v * n))
            for k in range(v):
                walls.add(Fraction(k, v))

        walls = sorted(walls)
        walls_ext = walls + [ONE]

        apex_aligned = 0
        apex_anti = 0
        apex_aligned_lonely = 0
        total_cells = 0

        for idx in range(len(walls)):
            if walls_ext[idx + 1] <= walls_ext[idx]:
                continue
            t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
            total_cells += 1

            # Apex: is runner 1 (speed 1) safe?
            runner1_safe = dist0(Fraction(1) * t_mid) >= thr

            # All runners safe?
            all_safe = all(dist0(Fraction(v) * t_mid) >= thr for v in speeds)

            if runner1_safe:
                apex_aligned += 1
                if all_safe:
                    apex_aligned_lonely += 1
            else:
                apex_anti += 1

        # Also: what fraction of time is the APEX the last blocker?
        # (i.e., all runners safe EXCEPT runner 1)
        apex_is_last_blocker = 0
        for idx in range(len(walls)):
            if walls_ext[idx + 1] <= walls_ext[idx]:
                continue
            t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2

            runner1_safe = dist0(Fraction(1) * t_mid) >= thr
            others_safe = all(dist0(Fraction(v) * t_mid) >= thr
                              for v in speeds if v != 1)

            if not runner1_safe and others_safe:
                apex_is_last_blocker += 1

        print(f"n={n} (initial segment):")
        print(f"  apex = arc (observer, runner_1=speed 1)")
        print(f"  aligned (runner 1 safe): {apex_aligned}/{total_cells} "
              f"({100*apex_aligned/total_cells:.1f}%)")
        print(f"  anti-aligned (runner 1 blocks): {apex_anti}/{total_cells} "
              f"({100*apex_anti/total_cells:.1f}%)")
        print(f"  lonely (all safe): {apex_aligned_lonely}/{total_cells}")
        print(f"  apex is LAST BLOCKER: {apex_is_last_blocker}/{total_cells} "
              f"({100*apex_is_last_blocker/total_cells:.1f}%)")
        print()

        # Compare with other runners as "last blocker"
        print(f"  last-blocker distribution:")
        for v in speeds:
            count = 0
            for idx in range(len(walls)):
                if walls_ext[idx + 1] <= walls_ext[idx]:
                    continue
                t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
                v_safe = dist0(Fraction(v) * t_mid) >= thr
                others = all(dist0(Fraction(w) * t_mid) >= thr
                             for w in speeds if w != v)
                if not v_safe and others:
                    count += 1
            if count > 0:
                print(f"    speed {v}: last blocker in {count} cells")
        print()


# ═══════════════════════════════════════════════════════════════
# PART 3: The apex's H-impact
# ═══════════════════════════════════════════════════════════════

def apex_H_impact(n_values=[4, 5, 6]):
    """How much does flipping the apex change H(T)?

    The apex connects the extremes of the ranking. Flipping it
    should have the LARGEST effect on H (by the deletion-contraction
    interpretation: it affects the most paths).
    """
    print("=" * 70)
    print("PART 3: The apex's impact on H(T)")
    print("=" * 70)
    print()

    for n in n_values:
        if n > 7:
            print(f"n={n}: too large for H computation")
            continue

        # Generate all tilings (2^m), compute H for each
        m = comb(n - 1, 2)

        # Tile positions (same as CLAUDE.md)
        tiles = []
        for y in range(1, n - 1):
            for x in range(n, y + 1, -1):
                if x - y >= 2:
                    tiles.append((x, y))

        # The apex is tile 0 = (n, 1)
        apex_idx = tiles.index((n, 1))

        def build_tournament(tiling_bits):
            adj = [[0] * n for _ in range(n)]
            # Base path arcs
            for k in range(n - 1):
                adj[k + 1][k] = 1  # (k+1) → k  (base path: n→n-1→...→1, 1-indexed)
            # Wait, vertices are 1-indexed. adj[i][j]=1 means i→j.
            # Base path n→n-1: adj[n-1][n-2] = 1 (0-indexed: vertex n-1 → vertex n-2)
            # Hmm, let me use 0-indexed vertices 0,...,n-1.
            # Base path: (n-1)→(n-2)→...→0
            adj = [[0] * n for _ in range(n)]
            for k in range(n - 1):
                adj[k + 1][k] = 1  # vertex (k+1) → vertex k

            # Tiles
            for t_idx, (x, y) in enumerate(tiles):
                # x, y are 1-indexed in the staircase
                # Convert to 0-indexed: vertex x-1 and vertex y-1
                vx, vy = x - 1, y - 1
                if (tiling_bits >> t_idx) & 1:
                    adj[vx][vy] = 1  # aligned: x→y (higher→lower)
                else:
                    adj[vy][vx] = 1  # anti-aligned: y→x

            return adj

        def count_ham_paths(adj, n):
            """Count Hamiltonian paths by DP."""
            dp = [[0] * n for _ in range(1 << n)]
            for v in range(n):
                dp[1 << v][v] = 1
            for mask in range(1 << n):
                for v in range(n):
                    if dp[mask][v] == 0:
                        continue
                    for u in range(n):
                        if mask & (1 << u):
                            continue
                        if adj[v][u]:
                            dp[mask | (1 << u)][u] += dp[mask][v]
            full = (1 << n) - 1
            return sum(dp[full][v] for v in range(n))

        # For each tiling with apex aligned vs anti-aligned:
        # compute H and compare
        H_apex_aligned = []
        H_apex_anti = []

        for bits in range(1 << m):
            adj = build_tournament(bits)
            H = count_ham_paths(adj, n)

            if (bits >> apex_idx) & 1:
                H_apex_aligned.append(H)
            else:
                H_apex_anti.append(H)

        avg_aligned = sum(H_apex_aligned) / len(H_apex_aligned)
        avg_anti = sum(H_apex_anti) / len(H_apex_anti)

        print(f"n={n} (tiles: {m}, apex at index {apex_idx}):")
        print(f"  apex ALIGNED: avg H = {avg_aligned:.2f}, "
              f"range [{min(H_apex_aligned)}, {max(H_apex_aligned)}]")
        print(f"  apex ANTI-ALIGNED: avg H = {avg_anti:.2f}, "
              f"range [{min(H_apex_anti)}, {max(H_apex_anti)}]")
        print(f"  H ratio (anti/aligned): {avg_anti/avg_aligned:.4f}")

        # Compare with OTHER tiles: which tile has the largest H-impact?
        tile_impacts = []
        for t_idx, (x, y) in enumerate(tiles):
            H_on = []
            H_off = []
            for bits in range(1 << m):
                adj = build_tournament(bits)
                H = count_ham_paths(adj, n)
                if (bits >> t_idx) & 1:
                    H_on.append(H)
                else:
                    H_off.append(H)
            impact = abs(sum(H_on) / len(H_on) - sum(H_off) / len(H_off))
            tile_impacts.append(((x, y), impact))

        tile_impacts.sort(key=lambda x: -x[1])
        print(f"  tile H-impact ranking:")
        for (x, y), impact in tile_impacts[:5]:
            is_apex = "*** APEX ***" if (x, y) == (n, 1) else ""
            print(f"    ({x},{y}): impact={impact:.2f} {is_apex}")
        print()


# ═══════════════════════════════════════════════════════════════
# PART 4: The apex in the cascade
# ═══════════════════════════════════════════════════════════════

def apex_in_cascade(n_values=[5, 7, 14]):
    """In the cascade proof, the apex is processed at a specific step.

    The cascade processes runners from slowest to fastest.
    The APEX connects the observer to the slowest runner (speed 1).
    So the apex is processed FIRST in the cascade!

    But wait — in the LRC cascade (S527), we process runners by their
    speed: v_1=1 first, then v_2=2, etc. The FIRST runner constrains
    t ∈ [1/n, (n-1)/n]. This IS the apex tile's constraint!

    The apex tile determines the FIRST constraint in the cascade.
    All subsequent tiles (interior diagonals) operate WITHIN this
    constraint.

    The cascade says: the apex sets the stage (the initial feasible set),
    and the interior tiles refine it.
    """
    print("=" * 70)
    print("PART 4: The apex tile in the cascade")
    print("=" * 70)
    print()

    for n in n_values:
        box_width = (n - 2) / n

        print(f"n={n}:")
        print(f"  APEX tile = (observer, runner_1) where runner_1 has speed 1")
        print(f"  apex aligned iff ||t|| ≥ 1/{n}, i.e., t ∈ [1/{n}, {n-1}/{n}]")
        print(f"  initial feasible set: width = {n-2}/{n} = {box_width:.4f}")
        print()
        print(f"  CASCADE STEP 1 (apex): constrains t to [{1/n:.4f}, {(n-1)/n:.4f}]")
        print(f"    This is {box_width*100:.1f}% of the circle.")

        # The apex is the COARSEST constraint.
        # Runner 2 (speed 2) adds a FINER constraint within this.
        # Runner 3 adds an even finer one.
        # The LAST runner (speed n-1) adds the finest constraint.

        # The OUTSIDE (apex + remaining boundary) gives the coarse structure.
        # The INSIDE (diagonals) gives the fine structure.

        # At each cascade step k:
        #   the constraint from speed k has granularity 1/(n*k)
        #   the feasible set's measure ≈ ((n-2)/n)^k

        print(f"  CASCADE STEPS:")
        mu = box_width
        for k in range(1, min(n, 10)):
            v = k + 1  # speed of the (k+1)-th runner
            granularity = 1.0 / (n * v)
            image_len = v * mu
            wraps = image_len >= 1

            print(f"    step {k+1}: speed {v}, granularity 1/{n*v}={granularity:.4f}, "
                  f"image_len={image_len:.3f} {'✓' if wraps else '✗'}")
            if wraps:
                mu *= box_width
            else:
                break
        print()

    print("THE APEX'S ROLE:")
    print("  The apex tile is processed FIRST in the cascade.")
    print("  It sets the COARSEST constraint: t ∈ [1/n, (n-1)/n].")
    print("  This is the polygon's OUTSIDE — the boundary minus one edge.")
    print("  All subsequent constraints operate INSIDE this boundary.")
    print()
    print("  The cascade is a REFINEMENT from outside to inside:")
    print("  Step 1 (apex, speed 1): the coarsest constraint, width (n-2)/n")
    print("  Step 2 (speed 2): finer, granularity 1/(2n)")
    print("  Step k (speed k): finest, granularity 1/(kn)")
    print()
    print("  The LAST step (speed n-1) is the finest constraint.")
    print("  It corresponds to the SHORTEST diagonal in the polygon")
    print("  (connecting near-adjacent vertices).")
    print()
    print("  OUTSIDE → INSIDE = COARSE → FINE = SLOW → FAST")
    print("  The cascade refines from the polygon boundary to the interior.")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 5: The Hamiltonian cycle and the apex
# ═══════════════════════════════════════════════════════════════

def hamiltonian_cycle_and_apex():
    """When the apex is ANTI-ALIGNED, the polygon boundary is a
    Hamiltonian CYCLE. This cycle is the strongest possible constraint:
    it means every vertex is reachable from every other.

    For LRC: when 1→n (sink beats source), the tournament has a
    Ham cycle through the boundary. The observer is NOT a source
    (the sink beats it). But the observer is ON the cycle.

    The LRC walk must FLIP the apex from anti-aligned to aligned
    to make the observer a source. This means: the slowest runner
    must cross from "close" to "safe."

    The TIME to flip the apex: runner 1 (speed 1) has close-intervals
    of width 2/(n·1) = 2/n centered at t=0, 1. The first crossing
    from close to safe happens at t = 1/n.

    After the apex flips: the feasible set is [1/n, (n-1)/n], and
    the cascade must check all remaining runners are safe within this.
    """
    print("=" * 70)
    print("PART 5: The Hamiltonian cycle — apex anti-aligned")
    print("=" * 70)
    print()

    print("When the APEX is anti-aligned (1 → n, sink beats source):")
    print("  The polygon boundary forms a HAMILTONIAN CYCLE:")
    print("  n → n-1 → ... → 2 → 1 → n")
    print("  (base path: n→...→1, plus apex: 1→n)")
    print()
    print("  Every vertex reaches every other via the cycle.")
    print("  The tournament is 'maximally cyclic' at the boundary level.")
    print("  H is MAXIMIZED for this boundary (the cycle contributes")
    print("  many paths that wrap around).")
    print()
    print("When the APEX is aligned (n → 1, source beats sink):")
    print("  The boundary is a DIRECTED PATH (no cycle):")
    print("  n → n-1 → ... → 1, and n → 1")
    print("  The source n beats ALL others via the path + direct arc.")
    print("  The tournament is 'maximally transitive' at the boundary level.")
    print("  H is MINIMIZED for this boundary (only one path: the base path).")
    print()
    print("THE LRC TRANSITION:")
    print("  At t < 1/n: runner 1 is close → apex ANTI-ALIGNED → cyclic boundary")
    print("  At t = 1/n: runner 1 crosses threshold → apex FLIPS → transitive boundary")
    print("  At t > 1/n: runner 1 is safe → apex ALIGNED → observer can be source")
    print()
    print("  The FIRST event in the cascade (apex flipping at t=1/n) is the")
    print("  TRANSITION from cyclic to transitive boundary structure.")
    print("  All subsequent cascade steps happen WITHIN the transitive boundary.")
    print()
    print("  This is why the initial segment has its lonely time at t=1/n:")
    print("  it's exactly when the apex flips, opening the transitive window.")
    print()


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    print("The Apex Tile — opus-S529")
    print()

    apex_tile_basics()
    apex_on_lrc_walk()
    apex_H_impact()
    apex_in_cascade()
    hamiltonian_cycle_and_apex()

    print("=" * 70)
    print("SYNTHESIS: The apex closes the polygon")
    print("=" * 70)
    print()
    print("The regular polygon has n boundary edges:")
    print("  n-1 BASE PATH edges (the ranking)")
    print("  1 APEX TILE (closing the cycle)")
    print()
    print("The apex tile is the MOST IMPORTANT tile because:")
    print("  1. It determines the global topology (path vs cycle)")
    print("  2. It connects the two extremes (source to sink)")
    print("  3. It's the FIRST constraint in the cascade (coarsest)")
    print("  4. Its flip at t=1/n opens the transitive window")
    print("  5. Its H-impact is among the largest of all tiles")
    print()
    print("For LRC: the apex = 'is the slowest runner safe?'")
    print("  When the apex flips from anti-aligned (cyclic) to aligned")
    print("  (transitive), the observer's source potential opens up.")
    print("  The remaining cascade steps check the interior tiles.")
    print()
    print("THE PICTURE:")
    print("  t < 1/n: apex CLOSED → cyclic boundary → observer blocked")
    print("  t = 1/n: apex OPENS → transition to transitive boundary")
    print("  t > 1/n: apex OPEN → cascade checks interior → lonely or not")
    print()
    print("  The polygon's OUTSIDE (base path + apex) sets the STAGE.")
    print("  The polygon's INSIDE (diagonals) determines the OUTCOME.")
    print("  The apex is the GATE between these two worlds.")
    print()


if __name__ == "__main__":
    main()
