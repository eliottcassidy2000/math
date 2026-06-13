#!/usr/bin/env python3
"""Section Tournament: nodes = arcs of the circle, edges = occupancy relations.

opus-2026-06-01-S539

THE IDEA: divide the circle into n equal SECTIONS of width 1/n.
Section k = [k/n, (k+1)/n). The observer is at position 0, straddling
the boundary between sections 0 and n-1.

At time t, each runner occupies exactly one section. The OCCUPANCY VECTOR
c = (c_0, ..., c_{n-1}) counts runners per section.

TOURNAMENT ON SECTIONS: orient section_i → section_j iff c_i > c_j
(the more crowded section "beats" the less crowded).

LONELY CONDITION: no runner in the observer's neighborhood.
The observer's danger zone is [0, 1/n) ∪ ((n-1)/n, 1] = sections 0 and n-1.
Lonely iff c_0 = 0 AND c_{n-1} = 0.

THE DYNAMICS: runner at speed v crosses from section k to section k+1
at time (k+1)/(nv). This is a CHIP MOVE on the circular section graph.

KEY QUESTIONS:
1. How many distinct section tournament iso classes are realizable?
2. How does this compare with the standard half-turn mapping?
3. What structure do the lonely section tournaments have?
4. Does the chip-firing dynamics constrain the section tournament walk?
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, floor
from functools import reduce
from itertools import combinations, permutations
from collections import Counter, defaultdict


ONE = Fraction(1)
ZERO = Fraction(0)


def frac(x):
    return x - Fraction(x.numerator // x.denominator)


def dist0(x):
    f = frac(x)
    return min(f, ONE - f)


def section_of(pos, n):
    """Which section does position pos ∈ [0,1) belong to?
    Section k = [k/n, (k+1)/n). Returns k.
    """
    return int(pos * n) % n


def occupancy_vector(speeds, n, t):
    """Compute the occupancy vector at time t."""
    occ = [0] * n
    for v in speeds:
        pos = frac(Fraction(v) * t)
        sec = section_of(pos, n)
        occ[sec] += 1
    return tuple(occ)


def section_tournament(occ, n):
    """Build tournament on n sections from occupancy vector.
    Section i → section j iff c_i > c_j.
    Tie-break: i → j if i has the runner closest to section boundary
    (proxy: i < j for simplicity).
    """
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if occ[i] > occ[j]:
                adj[i][j] = 1
            elif occ[i] == occ[j]:
                adj[i][j] = 1 if i < j else 0
    return tuple(tuple(r) for r in adj)


def canonicalize(adj, n):
    best = adj
    for perm in permutations(range(n)):
        new = tuple(tuple(adj[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if new < best:
            best = new
    return best


def pointed_canon(adj, n, obs_sec=0):
    """Canonicalize with observer's section fixed."""
    # The observer straddles sections 0 and n-1.
    # Fix section 0 (or the observer's position). Permute only sections 1..n-1.
    best = adj
    for perm in permutations(range(1, n)):
        full = (0,) + perm
        new = tuple(tuple(adj[full[i]][full[j]] for j in range(n)) for i in range(n))
        if new < best:
            best = new
    return best


def compute_section_walls(speeds, n):
    """Walls where a runner crosses a section boundary.
    Runner v crosses section boundary k/n at time k/(nv) for integer k.
    """
    walls = set([ZERO])
    for v in speeds:
        for k in range(1, n * v):
            t = Fraction(k, n * v)
            if ZERO < t < ONE:
                walls.add(t)
    return sorted(walls)


# ═══════════════════════════════════════════════════════════════
# PART 1: Section tournament class count
# ═══════════════════════════════════════════════════════════════

def section_class_count(n_values=[4, 5, 6, 7]):
    """Count distinct section tournament iso classes."""
    print("=" * 70)
    print("PART 1: Section tournament iso class counts")
    print("=" * 70)
    print()

    for n in n_values:
        max_speed = {4: 15, 5: 12, 6: 10, 7: 9}[n]

        sec_classes = set()
        sec_pointed = set()
        occ_vectors = set()
        lonely_sec = set()
        lonely_occ = set()
        total_cells = 0
        lonely_cells = 0

        for combo in combinations(range(1, max_speed + 1), n - 1):
            if reduce(gcd, combo) != 1:
                continue
            speeds = combo

            walls = compute_section_walls(speeds, n)
            walls_ext = walls + [ONE]

            for idx in range(len(walls)):
                if walls_ext[idx + 1] <= walls_ext[idx]:
                    continue
                t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
                total_cells += 1

                occ = occupancy_vector(speeds, n, t_mid)
                occ_vectors.add(occ)

                adj = section_tournament(occ, n)
                if n <= 7:
                    canon = canonicalize(adj, n)
                    pointed = pointed_canon(adj, n)
                    sec_classes.add(canon)
                    sec_pointed.add(pointed)

                # Check lonely: c_0 = 0 AND c_{n-1} = 0
                lonely = (occ[0] == 0 and occ[n - 1] == 0)
                if lonely:
                    lonely_cells += 1
                    lonely_sec.add(canon if n <= 7 else adj)
                    lonely_occ.add(occ)

        A000568 = {4: 4, 5: 12, 6: 56, 7: 456}

        print(f"n={n}:")
        print(f"  cells: {total_cells}, lonely: {lonely_cells}")
        print(f"  distinct occupancy vectors: {len(occ_vectors)}")
        print(f"  section tournament classes (unpointed): {len(sec_classes)}")
        print(f"  section tournament classes (pointed): {len(sec_pointed)}")
        print(f"  A000568({n}): {A000568.get(n, '?')}")
        print(f"  RESTRICTION: {len(sec_classes)}/{A000568.get(n,1)} = "
              f"{len(sec_classes)/A000568.get(n,1):.1%}")
        print(f"  lonely section classes: {len(lonely_sec)}")
        print(f"  lonely occupancy vectors: {len(lonely_occ)}")
        print()

        # Show the lonely occupancy vectors
        print(f"  LONELY OCCUPANCY VECTORS (sample):")
        for occ in sorted(lonely_occ)[:8]:
            print(f"    {occ}")
        print()


# ═══════════════════════════════════════════════════════════════
# PART 2: The chip-firing structure
# ═══════════════════════════════════════════════════════════════

def chip_firing_analysis(n_values=[4, 5, 6]):
    """Analyze the section occupancy as a chip-firing process.

    At t=0: all chips in section 0 (all runners at observer).
    Chips move clockwise: runner v crosses section boundaries at rate v.

    The chip-firing dynamics: at each wall, one chip moves from section k
    to section k+1 (mod n). The number of moves per period = Σ v_i.

    Key: the occupancy vector traces a PATH through the partition lattice.
    The lonely states (c_0 = c_{n-1} = 0) are specific vertices of this lattice.
    """
    print("=" * 70)
    print("PART 2: Chip-firing dynamics on sections")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))

        walls = compute_section_walls(speeds, n)
        walls_ext = walls + [ONE]

        print(f"n={n}, initial segment speeds={speeds}")
        print(f"  section boundary crossings (walls): {len(walls)}")
        print(f"  total chip moves per period: Σ v_i = {sum(speeds)}")
        print()

        # Trace the occupancy walk
        occ_walk = []
        for idx in range(len(walls)):
            if walls_ext[idx + 1] <= walls_ext[idx]:
                continue
            t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
            occ = occupancy_vector(speeds, n, t_mid)
            lonely = (occ[0] == 0 and occ[n - 1] == 0)
            occ_walk.append((float(t_mid), occ, lonely))

        # Show the walk
        print(f"  occupancy walk (first 15 states):")
        for t, occ, lonely in occ_walk[:15]:
            mark = " *** LONELY ***" if lonely else ""
            print(f"    t≈{t:.4f}: {occ}{mark}")

        # Find lonely states
        lonely_times = [(t, occ) for t, occ, lonely in occ_walk if lonely]
        print(f"\n  lonely states: {len(lonely_times)}")
        for t, occ in lonely_times[:5]:
            print(f"    t≈{t:.4f}: {occ}")

        # The TRANSITION GRAPH on occupancy vectors
        transitions = set()
        for i in range(len(occ_walk) - 1):
            transitions.add((occ_walk[i][1], occ_walk[i + 1][1]))

        # How many vectors differ by exactly one chip move?
        one_chip_moves = 0
        for a, b in transitions:
            diff = sum(abs(a[k] - b[k]) for k in range(n))
            if diff == 2:  # one chip moved: one section -1, adjacent section +1
                one_chip_moves += 1

        print(f"\n  transitions: {len(transitions)}")
        print(f"  single-chip moves: {one_chip_moves}/{len(transitions)}")
        print()


# ═══════════════════════════════════════════════════════════════
# PART 3: The occupancy vector as the NATURAL LRC state
# ═══════════════════════════════════════════════════════════════

def occupancy_as_state(n_values=[4, 5, 6, 7]):
    """The occupancy vector (c_0,...,c_{n-1}) with Σ c_i = n-1 is the
    NATURAL state for LRC:

    1. It's a composition of n-1 into n parts
    2. The dynamics are chip moves (one chip hops clockwise)
    3. The lonely condition is c_0 = c_{n-1} = 0
    4. The number of states = C(2n-2, n-1) (stars and bars)

    But the REALIZABLE states are fewer — the chip-firing dynamics
    constrain which compositions are accessible from the initial state
    (all chips in section 0).
    """
    print("=" * 70)
    print("PART 3: Occupancy vector as the natural LRC state")
    print("=" * 70)
    print()

    for n in n_values:
        max_speed = {4: 15, 5: 12, 6: 10, 7: 9}[n]

        # Total possible compositions of n-1 into n non-negative parts
        from math import comb
        total_compositions = comb(2 * n - 2, n - 1)

        # Lonely compositions: c_0 = c_{n-1} = 0, Σ c_i = n-1 with i=1..n-2
        # = compositions of n-1 into n-2 parts
        lonely_compositions = comb(2 * n - 4, n - 3) if n >= 3 else 0

        # Realizable compositions (from actual speed sets)
        realized = set()
        realized_lonely = set()

        for combo in combinations(range(1, max_speed + 1), n - 1):
            if reduce(gcd, combo) != 1:
                continue
            speeds = combo
            walls = compute_section_walls(speeds, n)
            walls_ext = walls + [ONE]

            for idx in range(len(walls)):
                if walls_ext[idx + 1] <= walls_ext[idx]:
                    continue
                t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
                occ = occupancy_vector(speeds, n, t_mid)
                realized.add(occ)
                if occ[0] == 0 and occ[n - 1] == 0:
                    realized_lonely.add(occ)

            for t in walls:
                occ = occupancy_vector(speeds, n, t)
                realized.add(occ)
                if occ[0] == 0 and occ[n - 1] == 0:
                    realized_lonely.add(occ)

        print(f"n={n}:")
        print(f"  total compositions: C({2*n-2},{n-1}) = {total_compositions}")
        print(f"  lonely compositions: C({2*n-4},{n-3}) = {lonely_compositions}")
        print(f"  realized compositions: {len(realized)}")
        print(f"  realized lonely: {len(realized_lonely)}")
        print(f"  RESTRICTION: {len(realized)}/{total_compositions} = "
              f"{len(realized)/total_compositions:.1%}")
        print(f"  lonely RESTRICTION: {len(realized_lonely)}/{lonely_compositions} = "
              f"{len(realized_lonely)/lonely_compositions:.1%}" if lonely_compositions > 0 else "")
        print()

        # The LONELY occupancy patterns
        if realized_lonely:
            print(f"  lonely occupancy patterns:")
            for occ in sorted(realized_lonely)[:10]:
                # How many sections are empty?
                empty = sum(1 for c in occ if c == 0)
                max_occ = max(occ)
                print(f"    {occ}  (empty={empty}, max={max_occ})")
            print()


# ═══════════════════════════════════════════════════════════════
# PART 4: The section-tournament view of the cascade
# ═══════════════════════════════════════════════════════════════

def section_cascade(n_values=[5, 6, 7]):
    """Reinterpret the cascade proof in section-tournament language.

    The cascade processes runners from slowest to fastest.
    In section language:
    - Runner v=1 moves at speed 1: crosses one section boundary per 1/n time
    - Runner v=k moves at speed k: crosses k boundaries per 1/n time

    The cascade's image-wrapping condition:
    - After constraining the first k runners, the feasible set I_{1..k}
      has measure μ_k
    - The next runner (speed k+1) moves fast enough to "visit" all sections
      within I_{1..k} → the section occupancy changes enough to reach
      the lonely configuration

    In chip-firing terms:
    - Each fast runner is a fast-moving chip that explores many sections
    - If it moves fast enough, it can be in section j at some time
      within the feasible set → any target section is reachable
    """
    print("=" * 70)
    print("PART 4: Section-tournament view of the cascade")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))
        print(f"n={n}: sections of width 1/{n}")
        print()

        # For each runner: how many section boundaries does it cross per period?
        for v in speeds:
            crossings = v  # crosses v boundaries in [0, 1)
            sections_visited = v  # visits v distinct sections
            print(f"  runner v={v}: {crossings} crossings/period, "
                  f"visits {min(sections_visited, n)} of {n} sections")

        print()

        # The TOTAL section crossings per period = Σ v_i
        total = sum(speeds)
        print(f"  total crossings: {total}")
        print(f"  crossings per section: {total/n:.1f}")
        print(f"  lonely requires: sections 0,{n-1} simultaneously empty")
        print()

        # KEY: the slowest runner (v=1) visits EXACTLY 1 section per 1/n time.
        # It's in section k at time t ∈ [k/n, (k+1)/n).
        # It's in section 0 for t ∈ [0, 1/n) and in section n-1 for t ∈ [(n-1)/n, 1).
        # For lonely: need runner 1 NOT in section 0 or n-1.
        # Runner 1 is safe iff t ∈ [1/n, (n-1)/n) = sections 1 through n-2.
        # This is the APEX condition (c=1): runner 1 in the interior.

        print(f"  Runner v=1 (the apex):")
        print(f"    in section 0 for t ∈ [0, 1/{n})")
        print(f"    in section k for t ∈ [{k}/{n}, {k+1}/{n})")
        print(f"    safe (not in 0 or {n-1}) for t ∈ [1/{n}, {n-1}/{n})")
        print(f"    → the APEX condition: t ∈ [1/{n}, {n-1}/{n})")
        print()

        # The FASTEST runner (v=n-1) visits n-1 sections per period.
        # It's in section (n-1)*k mod n for the k-th section visit.
        # For n-1 coprime to n: it visits ALL n sections.
        g = gcd(n - 1, n)
        visits = n // g
        print(f"  Runner v={n-1} (the fastest):")
        print(f"    visits {visits} of {n} sections per period")
        print(f"    gcd({n-1},{n}) = {g}")
        if visits == n:
            print(f"    EQUIDISTRIBUTED: visits every section once per period")
        else:
            print(f"    visits only {visits} sections (misses {n - visits})")
        print()

    print("THE SECTION-CASCADE PICTURE:")
    print("  The cascade processes runners as chips:")
    print("  - Runner 1 (slowest): a slow chip that visits 1 section per 1/n time.")
    print("    The APEX condition = this chip is in the interior (not section 0 or n-1).")
    print("  - Runner k: visits k sections per 1/n time.")
    print("    If k is coprime to n: equidistributed across all sections.")
    print("    If k shares a factor with n: visits only n/gcd(k,n) sections.")
    print()
    print("  The cascade succeeds when:")
    print("  - The slow chip (apex) is in the interior (t ∈ [1/n, (n-1)/n])")
    print("  - The fast chips collectively AVOID sections 0 and n-1")
    print("  - This is possible because fast chips visit many sections")
    print("    and can be 'steered' into the interior by choosing t within")
    print("    the apex-conditioned feasible set.")
    print()
    print("  The chip-firing on sections IS the cascade in disguise:")
    print("  each cascade step constrains one chip, shrinking the feasible set.")
    print("  The image-wrapping condition = the chip visits enough sections")
    print("  to guarantee a feasible t exists.")
    print()


def main():
    print("Section Tournament: Nodes = Arcs of the Circle — opus-S539")
    print()

    section_class_count()
    chip_firing_analysis()
    occupancy_as_state()
    section_cascade()


if __name__ == "__main__":
    main()
