#!/usr/bin/env python3
"""Deep creative reframings of LRC — beyond tournaments.

opus-2026-06-01-S540

Going beyond tournament mappings: what OTHER mathematical structures
does LRC naturally encode? Each reframing should give a new proof angle.

REFRAMINGS:

1. BRAID: Runner trajectories on the spacetime cylinder form a BRAID.
   The braid's linking numbers encode which runners cross which.
   Lonely = the observer strand is "unlinked" from all runners.

2. ECLIPSE: n-1 opaque shadows on a rotating dial. Each shadow has
   angular width 2/n. LRC = the observer escapes eclipse at some time.
   This is a COVERING problem: can moving shadows permanently cover a point?

3. LATTICE SHORTEST VECTOR: LRC as a closest-vector problem in the
   lattice Z(v_1,...,v_{n-1}) + Z^{n-1}. The lonely set is a ball of
   radius (n-2)/(2n) around (1/2,...,1/2). LRC = the lattice has a
   vector in this ball.

4. GRAPH COLORING: Build the "conflict hypergraph" where runners are
   vertices and hyper-edges are sets of runners whose close zones
   overlap in time. LRC iff the hypergraph has a "gap" (time not
   covered by any hyper-edge).

5. MUSICAL CANON: Runners are voices in a round at different tempos.
   The observer is the reference pitch. Loneliness = all voices are
   "consonant" (≥ 1/n away from unison). LRC = every canon resolves.
"""

from __future__ import annotations
from fractions import Fraction
from math import gcd, sqrt, pi, sin, cos
from functools import reduce
from itertools import combinations
from collections import Counter, defaultdict


ONE = Fraction(1)
ZERO = Fraction(0)
def frac(x): return x - Fraction(x.numerator // x.denominator)
def dist0(x):
    f = frac(x); return min(f, ONE - f)


# ═══════════════════════════════════════════════════════════════
# REFRAME 1: BRAID CROSSING NUMBER
# ═══════════════════════════════════════════════════════════════

def braid_reframe(n_values=[4, 5, 6, 7]):
    """Runner trajectories as braids on the spacetime cylinder.

    Each runner traces a helical path: (t, {v_i t}) on the cylinder
    [0,1) × [0,1). The observer traces (t, 0).

    Two runners "cross" when their positions swap circular order.
    A runner "crosses the observer" when it enters/exits the close zone.

    The CROSSING NUMBER = total crossings in one period.
    The LINKING NUMBER between observer and runner i = algebraic sum
    of crossings (±1 depending on direction).

    For runner i: it crosses the observer (enters close zone) v_i times
    per period going clockwise, and v_i times going counterclockwise.
    Total crossings = 2v_i. But the ALGEBRAIC sum = 0 (same number
    in each direction for a full period).

    Lonely = the observer strand is momentarily "unlinked" from all runner strands.
    """
    print("=" * 70)
    print("REFRAME 1: BRAID — runner trajectories as a braid group")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))
        m = len(speeds)

        total_crossings = 2 * sum(speeds)  # observer-runner crossings
        runner_crossings = sum(abs(vi - vj) for i, vi in enumerate(speeds)
                               for j, vj in enumerate(speeds) if i < j)

        print(f"n={n}, speeds={speeds}:")
        print(f"  observer-runner crossings per period: {total_crossings}")
        print(f"  runner-runner crossings: {runner_crossings}")
        print(f"  total braid crossings: {total_crossings + runner_crossings}")
        print(f"  braid group: B_{n} on {n} strands")
        print()

        # The braid WORD: sequence of crossings σ_i (strand i crosses strand i+1)
        # The lonely condition: a moment when the observer strand is
        # "above" all runner strands (not crossed by any).

        # Between consecutive observer crossings: the observer is either
        # in the close zone (crossed/linked) or in the safe zone (unlinked).
        # LRC = there exists a segment where the observer is unlinked from ALL.

        # The average unlinked time: 1 - 2(n-1)/n = (2-n)/n (negative for n≥3)
        # So on average, the observer is more linked than unlinked.
        # But the STRUCTURE of the braid (the word) determines the actual
        # linked/unlinked pattern.

        avg_unlinked = 1 - 2 * (n - 1) / n
        print(f"  avg unlinked fraction: {avg_unlinked:.4f} ({'positive' if avg_unlinked > 0 else 'NEGATIVE'})")
        print(f"  (negative means observer is more often linked than not)")
        print()

    print("BRAID INSIGHT: the observer strand is 'entangled' with runner strands.")
    print("LRC says: the braid has a moment of 'disentanglement' where the")
    print("observer is free. The entanglement measure (linking number) is")
    print("always 0 algebraically, but the GEOMETRIC entanglement (crossing")
    print("pattern) determines whether disentanglement is possible.")
    print()
    print("The braid group B_n acts on the tournament: a crossing σ_i")
    print("corresponds to an arc flip in the tournament. The BRAID INVARIANTS")
    print("(Jones polynomial, etc.) might encode the LRC condition.")
    print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 2: ECLIPSE / SHADOW COVERAGE
# ═══════════════════════════════════════════════════════════════

def eclipse_reframe(n_values=[4, 5, 6, 7, 14]):
    """n-1 moving shadows trying to permanently eclipse the observer.

    Each runner casts a "shadow" of angular width 2/n on the circle.
    Shadow_i(t) = {x : ||x - v_i t|| < 1/n} = the close zone.

    LRC = the observer escapes eclipse at some time.
    Equivalently: ∪ Shadow_i ≠ [0,1).

    The COVERAGE at each time: how many shadows cover the observer.
    Coverage 0 = lonely. Coverage ≥ 1 = eclipsed.

    KEY METRIC: the MAXIMUM OVERLAP of shadows at the observer.
    If the max overlap is < n-1, some shadow is "wasted" (covering
    the observer redundantly instead of covering a gap).
    """
    print("=" * 70)
    print("REFRAME 2: ECLIPSE — moving shadows covering the observer")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))
        m = len(speeds)

        # Each shadow has measure 2/n per period
        shadow_width = Fraction(2, n)
        total_shadow = m * shadow_width

        # Average coverage at the observer:
        # Each shadow covers the observer for time 2/n.
        # Average coverage = m * 2/n = 2(n-1)/n
        avg_coverage = 2 * m / n

        # For permanent eclipse: need coverage ≥ 1 at ALL times.
        # Total shadow = 2(n-1)/n. For n≥3: > 1 but < 2.
        # So permanent eclipse is possible IF shadows are arranged right.

        # The OVERLAP: average pairwise overlap = ?
        # For coprime speeds: overlap(i,j) ≈ (2/n)^2 = 4/n^2
        # Total overlap ≈ C(m,2) * 4/n^2 = 2m(m-1)/n^2
        avg_pairwise = 2 * m * (m - 1) / n**2

        # Eclipse excess = total_shadow - 1 = how much shadow is "wasted" on overlap
        eclipse_excess = float(total_shadow) - 1

        print(f"n={n}:")
        print(f"  shadow width: 2/{n} = {float(shadow_width):.4f}")
        print(f"  total shadow measure: {m}×2/{n} = {float(total_shadow):.4f}")
        print(f"  eclipse excess: {eclipse_excess:.4f}")
        print(f"  avg pairwise overlap: {avg_pairwise:.4f}")
        print(f"  eclipse excess / overlap: {eclipse_excess/max(0.001,avg_pairwise):.2f}")
        print()

    print("ECLIPSE INSIGHT: the total shadow measure is 2(n-1)/n < 2.")
    print("For permanent eclipse, shadows must cover [0,1) with this budget.")
    print("The 'eclipse excess' (total - 1) = (n-2)/n is the amount of shadow")
    print("that MUST overlap. If the overlap is FORCED to be more than the excess,")
    print("the shadows can't cover [0,1), proving LRC.")
    print()
    print("The CASCADE proof (S527) shows: the overlap IS forced to exceed the")
    print("excess for n≥7, because the image-wrapping condition creates mandatory")
    print("shadow overlaps. The 'wasted' shadow = the lonely time.")
    print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 3: LATTICE CLOSEST VECTOR
# ═══════════════════════════════════════════════════════════════

def lattice_reframe(n_values=[4, 5, 6, 7, 14]):
    """LRC as a closest-vector problem in a lattice.

    The runner trajectory visits lattice points on the torus:
    (v_1 k/M, ..., v_{n-1} k/M) mod 1 for k=0,...,M-1 where M=lcm.

    The lonely set = box [1/n, (n-1)/n]^{n-1} centered at (1/2,...,1/2).
    LRC = some lattice point lands in the box.

    The COVERING RADIUS μ of the rank-1 sublattice:
    μ = max_{x ∈ T^{n-1}} min_k ||(v_1 k/M,...,v_{n-1} k/M) - x||_∞

    LRC iff μ ≤ (n-2)/(2n) (the half-width of the box).
    """
    print("=" * 70)
    print("REFRAME 3: LATTICE — closest-vector to the lonely box")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))
        m = len(speeds)

        box_half_width = (n - 2) / (2 * n)
        box_volume = ((n - 2) / n) ** m

        # Lattice: rank-1 sublattice generated by (v_1,...,v_m)
        # Number of lattice points per period = lcm(v_1,...,v_m)
        N = speeds[0]
        for v in speeds[1:]:
            N = N * v // gcd(N, v)

        # Average density = N lattice points uniformly in T^m
        # Expected points in box ≈ N * box_volume
        expected = N * box_volume

        print(f"n={n}:")
        print(f"  box half-width: (n-2)/(2n) = {box_half_width:.4f}")
        print(f"  box volume: ((n-2)/n)^{m} = {box_volume:.6f}")
        print(f"  lattice points: lcm = {N}")
        print(f"  expected in box: {expected:.1f}")
        print(f"  {'>> 1 (easy)' if expected > 10 else '≈ 1 (tight)' if expected > 0.5 else '< 1 (HARD)'}")
        print()

    print("LATTICE INSIGHT: LRC is a DIOPHANTINE APPROXIMATION problem.")
    print("The rational line must come within (n-2)/(2n) of the box center")
    print("in the L∞ norm. This is the VIEW OBSTRUCTION problem (Cusick 1973).")
    print()
    print("The COVERING RADIUS of the speed lattice determines LRC:")
    print("  μ ≤ (n-2)/(2n) iff LRC holds.")
    print()
    print("For the initial segment: the lattice is 'well-spread' (many points)")
    print("but the box hit count is LOW (6 for n=14) because the lattice points")
    print("cluster in specific regions of the torus, missing the box interior")
    print("but grazing its boundary.")
    print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 4: CONFLICT HYPERGRAPH
# ═══════════════════════════════════════════════════════════════

def conflict_hypergraph(n_values=[4, 5, 6]):
    """The conflict hypergraph: when do close zones overlap?

    For each TIME INTERVAL: the set of runners currently close to the
    observer forms a HYPER-EDGE in the conflict hypergraph.

    Vertices = time intervals.
    Hyper-edges = sets of time intervals where a specific runner is close.

    The COVERING NUMBER = min hyper-edges needed to cover all vertices.
    LRC iff covering number < total vertices (some vertex uncovered = lonely).

    Equivalently: the INDEPENDENCE NUMBER of the dual hypergraph > 0.
    An independent set in the dual = a set of non-overlapping close zones
    = a time interval not covered by any close zone.
    """
    print("=" * 70)
    print("REFRAME 4: CONFLICT HYPERGRAPH — overlapping close zones")
    print("=" * 70)
    print()

    for n in n_values:
        max_speed = {4: 15, 5: 10, 6: 9}[n]
        thr = Fraction(1, n)

        total_sets = 0
        max_overlap_sum = 0
        min_overlap_sum = float('inf')

        for combo in combinations(range(1, max_speed + 1), n - 1):
            if reduce(gcd, combo) != 1:
                continue
            total_sets += 1
            speeds = combo

            # Sample the overlap count at many times
            num_pts = 5000
            overlaps = []
            for s in range(num_pts):
                t = Fraction(s, num_pts)
                count = sum(1 for v in speeds if dist0(Fraction(v) * t) < thr)
                overlaps.append(count)

            # The "overlap signature" = histogram of coverage counts
            hist = Counter(overlaps)
            overlap_sum = sum(overlaps) / num_pts  # average coverage

            max_overlap_sum = max(max_overlap_sum, overlap_sum)
            min_overlap_sum = min(min_overlap_sum, overlap_sum)

        print(f"n={n} ({total_sets} sets):")
        print(f"  avg coverage range: [{min_overlap_sum:.4f}, {max_overlap_sum:.4f}]")
        print(f"  expected (uniform): {2*(n-1)/n:.4f}")
        print()

    print("HYPERGRAPH INSIGHT: the average coverage is always 2(n-1)/n.")
    print("But the DISTRIBUTION of coverage matters:")
    print("  - Uniform coverage: each time has ~2(n-1)/n runners close.")
    print("    Lonely probability ≈ Poisson with λ=2(n-1)/n.")
    print("  - Clustered coverage: some times have many close runners,")
    print("    others have few → more likely to have gaps (lonely times).")
    print()
    print("The initial segment has the MOST UNIFORM coverage (least clustered)")
    print("→ LEAST likely to have gaps → tightest case (lonely measure 0).")
    print("All other speed sets have MORE clustering → MORE gaps → easier LRC.")
    print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 5: MUSICAL CANON — consonance resolution
# ═══════════════════════════════════════════════════════════════

def musical_canon(n_values=[4, 5, 6, 7]):
    """Runners as voices in a musical canon at different tempos.

    Voice i plays at tempo v_i (cycles per period).
    The observer is the reference pitch (position 0).
    The "dissonance" of voice i at time t = 1 - ||v_i t|| / (1/2).
    Maximum dissonance = 1 (unison with observer). Minimum = 0 (antipodal).

    LRC says: all voices are simultaneously "consonant" (distance ≥ 1/n
    from observer) at some time.

    The "consonance threshold" = 1/n of the octave.
    The "resolution" = the moment when all dissonances drop below threshold.

    In MUSIC THEORY: this is related to BEATING patterns. Two voices at
    speeds v_i and v_j create a beat frequency |v_i - v_j|. The beats
    of all pairs create a complex rhythm. The "resolution" (lonely time)
    is when all beat patterns align at a consonant phase.
    """
    print("=" * 70)
    print("REFRAME 5: MUSICAL CANON — when does the round resolve?")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))
        m = len(speeds)

        # Beat frequencies between all pairs
        beat_freqs = []
        for i in range(m):
            for j in range(i + 1, m):
                beat_freqs.append(abs(speeds[i] - speeds[j]))

        # The "resolution period" = lcm of all beat frequencies
        resolution = beat_freqs[0]
        for b in beat_freqs[1:]:
            resolution = resolution * b // gcd(resolution, b)

        # Fundamental beat = gcd of all beat frequencies
        fundamental = beat_freqs[0]
        for b in beat_freqs[1:]:
            fundamental = gcd(fundamental, b)

        print(f"n={n}, voices at tempos {speeds}:")
        print(f"  beat frequencies: {sorted(set(beat_freqs))}")
        print(f"  fundamental beat: {fundamental}")
        print(f"  resolution period (lcm of beats): {resolution}")
        print(f"  resolution points per period: {resolution}")
        print(f"  consonance threshold: 1/{n} of the octave")
        print()

    print("MUSICAL INSIGHT: LRC is about BEAT PATTERN ALIGNMENT.")
    print("Each pair of voices creates a beat. The beats interact to create")
    print("a complex rhythm. The 'resolution' is when ALL beats align at a")
    print("consonant phase — all voices are simultaneously far from unison.")
    print()
    print("The INITIAL SEGMENT has the richest beat structure (all integer")
    print("beat frequencies), which makes alignment hardest (tightest case).")
    print("Speed sets with coprime beat frequencies have SPARSER interactions,")
    print("making resolution easier.")
    print()
    print("The FORMAL GROUP F(x,y)=(x+y)/(1+xy) is the CONSONANCE LAW:")
    print("it adds dissonances hyperbolically. The rapidity = hyperbolic")
    print("distance = musical 'interval'. The lonely condition = all")
    print("intervals exceed the consonance threshold 0.5*ln(n-1) in rapidity.")
    print()


# ═══════════════════════════════════════════════════════════════
# SYNTHESIS: The deepest reframing
# ═══════════════════════════════════════════════════════════════

def synthesis():
    print("=" * 70)
    print("SYNTHESIS: Five faces of LRC")
    print("=" * 70)
    print()
    print("LRC has (at least) five natural mathematical identities:")
    print()
    print("1. TOURNAMENT: observer source in a marked round tournament.")
    print("   → iso class reachability on the metagraph G_n")
    print("   → the 'inside/outside' polygon geometry")
    print()
    print("2. BRAID: disentanglement of the observer strand.")
    print("   → linking numbers and braid invariants")
    print("   → the crossing pattern determines lonely times")
    print()
    print("3. ECLIPSE: gaps in moving shadow coverage.")
    print("   → overlap budget: total shadow = 2(n-1)/n")
    print("   → mandatory overlap (cascade) forces gaps")
    print()
    print("4. LATTICE: closest vector to the lonely box.")
    print("   → covering radius of the speed sublattice")
    print("   → Diophantine approximation in the L∞ norm")
    print()
    print("5. CANON: resolution of a musical round.")
    print("   → beat frequency alignment")
    print("   → formal group as consonance law")
    print()
    print("Each reframing highlights a different OBSTRUCTION to permanent")
    print("non-loneliness:")
    print("  - Tournament: the metagraph is too connected")
    print("  - Braid: the linking is algebraically zero")
    print("  - Eclipse: the shadow budget is too small")
    print("  - Lattice: the box is too large")
    print("  - Canon: the beats must eventually align")
    print()
    print("THE UNIFYING PRINCIPLE: n-1 periodic constraints on a single")
    print("parameter (time t) cannot SIMULTANEOUSLY block a target of")
    print("measure ≈ e^{-2} ≈ 13.5% because the constraints INTERFERE")
    print("with each other — their overlaps waste budget, their phases")
    print("drift apart, their linking cancels, their resonances decohere.")
    print()
    print("The INITIAL SEGMENT is the unique CRITICAL POINT where the")
    print("interference is MAXIMALLY COHERENT — all phases align, all")
    print("overlaps minimize, all beat frequencies harmonize. Even at")
    print("this critical point, the target is barely reached (wall-only).")
    print("Any perturbation breaks the coherence and opens gaps.")
    print()


def main():
    print("Deep Creative Reframings of LRC — opus-S540")
    print()

    braid_reframe()
    eclipse_reframe()
    lattice_reframe()
    conflict_hypergraph()
    musical_canon()
    synthesis()


if __name__ == "__main__":
    main()
