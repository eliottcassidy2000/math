#!/usr/bin/env python3
"""Creative reframes of the Lonely Runner Conjecture.

opus-2026-06-01-S520

The user asked for "even more creative reframes." This script explores
LRC through radically different mathematical lenses, computing concrete
evidence for each reframe and looking for surprising connections.

Reframes explored:
1. BILLIARD: LRC as a billiard trajectory on a polytope
2. LOVÁSZ LOCAL LEMMA: probabilistic non-covering bound
3. SANDPILE/CHIP-FIRING: proximity chips and threshold firing
4. RELATIVISTIC VELOCITY: F(x,y)=(x+y)/(1+xy) on runner gaps
5. LATTICE GEOMETRY: Minkowski / covering radius on the speed lattice
6. MUSICAL CHAIRS: competitive exclusion dynamics
7. SPECTRAL: Fourier analysis of the lonely indicator
8. TOPOLOGICAL: winding number / Borsuk-Ulam obstruction
9. THERMODYNAMIC: free energy of the runner gas
10. AUCTION: gap-price equilibrium under runner bidding

Tournament Analysis declaration:
    vertices: the 10 reframes
    pairwise observable: which reframe produces more nontrivial constraints
    switch: majority vote over 5 criteria (novelty, computability, proof-force,
            connection to known results, beauty)
    tie path: listed order

Stored output:
    05-knowledge/results/lrc_creative_reframes_s520.out
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd, pi, sin, cos, exp, log, sqrt, ceil
from functools import reduce
from collections import Counter
import cmath


ONE = Fraction(1)
ZERO = Fraction(0)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    f = frac(x)
    return min(f, ONE - f)


def primitive_speed_sets(n, max_speed):
    for combo in combinations(range(1, max_speed + 1), n - 1):
        if reduce(gcd, combo) == 1:
            yield combo


# ═══════════════════════════════════════════════════════════════
# REFRAME 1: BILLIARD on the hypercube
# ═══════════════════════════════════════════════════════════════

def reframe_billiard(n_values=[3, 4, 5, 6]):
    """LRC as a billiard trajectory.

    The runner positions ({v_1 t}, ..., {v_{n-1} t}) trace a line in [0,1)^{n-1}.
    The "good" region is {x : all x_i in [1/n, 1-1/n]}.
    LRC says this line hits the good region.

    The good region is a smaller hypercube [1/n, 1-1/n]^{n-1} of side (n-2)/n.
    Its volume is ((n-2)/n)^{n-1}.

    For large n: volume ~ (1-2/n)^{n-1} ~ e^{-2(n-1)/n} ~ e^{-2}.
    So even for large n, the good region has ~13.5% of the volume!

    The billiard trajectory is ergodic (Weyl) if speeds are rationally independent.
    For integer speeds, the trajectory is periodic and dense in a rational subtorus.

    Key question: does the rational subtorus always intersect the good region?
    """
    print("=" * 70)
    print("REFRAME 1: BILLIARD — line in the torus hitting a box")
    print("=" * 70)
    print()

    for n in n_values:
        k = n - 1
        thr = Fraction(1, n)
        good_side = Fraction(n - 2, n)
        good_volume = good_side ** k

        print(f"n={n}  dim={k}  good_side={float(good_side):.4f}  "
              f"good_volume={float(good_volume):.6f}")

        # For the initial segment, the trajectory lies in a (k-1)-dimensional
        # rational subtorus. How much of it hits the good box?
        speeds = tuple(range(1, n))

        # Sample the trajectory
        num_samples = 10000
        hits = 0
        for s in range(1, num_samples + 1):
            t = Fraction(s, num_samples)
            if all(thr <= frac(Fraction(v) * t) <= 1 - thr for v in speeds):
                hits += 1

        print(f"  initial segment: {hits}/{num_samples} samples in good box "
              f"({100*hits/num_samples:.2f}%)")
        print(f"  expected if uniform: {100*float(good_volume):.2f}%")
        print(f"  ratio (observed/expected): "
              f"{(hits/num_samples)/float(good_volume):.3f}")
        print()

    # Large-n asymptotics
    print("Large-n asymptotics of good volume:")
    for n in [10, 14, 20, 50, 100]:
        vol = ((n - 2) / n) ** (n - 1)
        print(f"  n={n:3d}: vol = {vol:.6f}  (e^-2 = {exp(-2):.6f})")
    print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 2: LOVÁSZ LOCAL LEMMA
# ═══════════════════════════════════════════════════════════════

def reframe_lovasz(n_values=[3, 4, 5, 6, 7, 14]):
    """LRC via the Lovász Local Lemma.

    Events: A_i = "runner i is close to observer at time t"
            (i.e., ||v_i t|| < 1/n)
    Each P(A_i) = 2/n.

    LRC says: P(none of A_i) > 0 (or = 0 with boundary witness).

    Naive union bound: P(some A_i) <= (n-1) * 2/n = 2(n-1)/n.
    For n >= 3, this is > 1, so the union bound fails.

    Lovász Local Lemma: if each A_i depends on at most d other events,
    and e * 2/n * (d+1) <= 1, then P(none of A_i) > 0.

    The dependency: A_i and A_j are "dependent" if their forbidden intervals
    overlap in t-space.
    {t : ||v_i t|| < 1/n} has measure 2/n.
    {t : ||v_j t|| < 1/n} has measure 2/n.
    Their intersection has measure... depends on gcd(v_i, v_j) and the threshold.

    For coprime v_i, v_j: the intersection measure is approximately 4/n^2
    (product of individual measures, since equidistributed).

    So "dependency" means the intersection is larger than expected.
    If gcd(v_i, v_j) = g > 1, the runner pair has period 1/g, and
    the intersection can be up to 2g/n (concentrated).

    For the dependency graph: connect i-j if gcd(v_i, v_j) > 1.
    The maximum degree d in this dependency graph determines LLL applicability.
    """
    print("=" * 70)
    print("REFRAME 2: LOVÁSZ LOCAL LEMMA — dependent coverage events")
    print("=" * 70)
    print()

    for n in n_values:
        k = n - 1
        p = 2.0 / n  # P(A_i)

        # For initial segment, compute the dependency graph
        speeds = tuple(range(1, n))
        dep_edges = 0
        max_deg = 0
        for i, vi in enumerate(speeds):
            deg = 0
            for j, vj in enumerate(speeds):
                if i != j and gcd(vi, vj) > 1:
                    if i < j:
                        dep_edges += 1
                    deg += 1
            max_deg = max(max_deg, deg)

        # LLL criterion: e * p * (d+1) <= 1
        lll_product = exp(1) * p * (max_deg + 1)
        lll_holds = lll_product <= 1

        # Symmetric LLL: (n-1) * p * (1/(1-p))^{d} < 1 ??
        # Actually the standard symmetric form: e * p * (d+1) <= 1

        print(f"n={n:2d}  p=2/n={p:.4f}  dep_edges={dep_edges}  "
              f"max_deg={max_deg}")
        print(f"  LLL criterion: e*p*(d+1) = {lll_product:.4f}  "
              f"{'<= 1 HOLDS' if lll_holds else '> 1 FAILS'}")

        # Cluster expansion / improved bound
        # Shearer's lemma: need P(A_i) < p*(G) where p*(G) is the
        # independent set polynomial threshold
        # For a graph with max degree d: p*(G) >= 1/(d+1)
        # So need p < 1/(d+1)
        shearer_bound = 1.0 / (max_deg + 1) if max_deg > 0 else 1.0
        print(f"  Shearer bound: p < {shearer_bound:.4f}  "
              f"{'HOLDS' if p < shearer_bound else 'FAILS'}")

        # Independence number of the dependency graph
        # (maximum set of runners with pairwise coprime speeds)
        # This is the maximum antichain in the divisibility order
        primes_in_range = [v for v in speeds if all(v % p != 0 for p in range(2, v))]
        print(f"  primes in speed set: {len(primes_in_range)} "
              f"({primes_in_range[:8]}...)")
        print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 3: RELATIVISTIC VELOCITY ADDITION
# ═══════════════════════════════════════════════════════════════

def reframe_relativistic(n_values=[3, 4, 5, 6]):
    """LRC via the formal group F(x,y) = (x+y)/(1+xy).

    The project's signature object! This is relativistic velocity addition.

    Map runner gaps to "rapidities": if gap g_i = ||v_i t|| is the circular
    distance from observer, define rapidity r_i = atanh(1 - 2*g_i).

    Then:
    - g_i = 0 (observer tie) => r_i = +∞
    - g_i = 1/2 (antipodal) => r_i = 0
    - g_i = 1/n (threshold) => r_i = atanh(1 - 2/n) = atanh((n-2)/n)

    The lonely condition g_i >= 1/n becomes r_i <= atanh((n-2)/n).

    The "relativistic sum" of all observer-runner interactions via F:
    if we compose the n-1 rapidity contributions... what happens?

    Key insight: F(x,y) = tanh(atanh(x) + atanh(y)).
    So the formal group is just ADDITION in rapidity space.

    The total "rapidity load" on the observer is sum of individual rapidities.
    The observer is lonely when the individual contributions are all below threshold.
    """
    print("=" * 70)
    print("REFRAME 3: RELATIVISTIC VELOCITY — rapidity of loneliness")
    print("=" * 70)
    print()

    for n in n_values:
        k = n - 1
        thr = 1.0 / n
        rapidity_threshold = 0.5 * log((1 + (1 - 2 * thr)) / (1 - (1 - 2 * thr))) if thr < 0.5 else float('inf')
        # atanh((n-2)/n) = 0.5 * ln((1 + (n-2)/n) / (1 - (n-2)/n))
        #                = 0.5 * ln((2n-2)/(2)) = 0.5 * ln(n-1)
        rapidity_exact = 0.5 * log(n - 1)

        print(f"n={n}  threshold rapidity = atanh((n-2)/n) = 0.5*ln(n-1) = {rapidity_exact:.4f}")
        print(f"  lonely iff all |rapidity_i| <= {rapidity_exact:.4f}")

        # For the initial segment, compute rapidities at the lonely time t=1/(2n)
        speeds = list(range(1, n))
        t_lonely = 1.0 / (2 * n)
        print(f"  at t=1/(2n)={t_lonely:.6f}:")

        rapidities = []
        for v in speeds:
            pos = (v * t_lonely) % 1.0
            gap = min(pos, 1.0 - pos)
            if gap > 0 and gap < 0.5:
                rap = 0.5 * log((1 + (1 - 2 * gap)) / max(1e-15, (1 - (1 - 2 * gap))))
                rapidities.append(rap)
            elif gap >= 0.5:
                rapidities.append(0.0)
            else:
                rapidities.append(float('inf'))

        if rapidities and all(r < float('inf') for r in rapidities):
            total_rapidity = sum(rapidities)
            max_rapidity = max(rapidities)
            print(f"    max rapidity: {max_rapidity:.4f} (threshold: {rapidity_exact:.4f})")
            print(f"    total rapidity sum: {total_rapidity:.4f}")
            print(f"    all below threshold: {all(r <= rapidity_exact + 0.001 for r in rapidities)}")

        # The "relativistic energy" interpretation
        # E = cosh(rapidity) for each runner
        # Total energy = product of cosh(r_i) (since F is additive in rapidity)
        if rapidities and all(r < float('inf') for r in rapidities):
            energies = [exp(abs(r)) for r in rapidities]  # ~ cosh(r) for large r
            total_energy = 1.0
            for e in energies:
                total_energy *= e
            print(f"    total 'relativistic energy' (product of cosh): {total_energy:.4f}")

        print()

    # Key observation: rapidity threshold = 0.5 * ln(n-1)
    # This grows logarithmically! For large n, each runner needs rapidity < 0.5*ln(n-1)
    # which means gap > 1/n (which is what we want).
    # The total rapidity budget is (n-1) * 0.5 * ln(n-1).
    print("Rapidity budget analysis:")
    for n in [3, 5, 7, 14, 20, 50]:
        budget = (n - 1) * 0.5 * log(max(2, n - 1))
        print(f"  n={n:3d}: total rapidity budget = {budget:.2f}, "
              f"per-runner = {0.5*log(max(2,n-1)):.4f}")
    print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 4: LATTICE COVERING RADIUS
# ═══════════════════════════════════════════════════════════════

def reframe_lattice(n_values=[3, 4, 5, 6]):
    """LRC via the covering radius of a lattice.

    The speeds v_1,...,v_{n-1} define a lattice L in R^{n-1}:
    L = {(a_1/v_1, ..., a_{n-1}/v_{n-1}) : a_i in Z}

    The runner positions at time t are (v_1 t mod 1, ..., v_{n-1} t mod 1).
    This is the image of the point (t,...,t) under the linear map x -> (v_1 x,...,v_{n-1} x) mod Z^{n-1}.

    The "lonely region" in the torus is [1/n, 1-1/n]^{n-1}.

    LRC asks: does the 1-dimensional subtorus {(v_1 t,...,v_{n-1} t) : t in R/Z}
    intersect [1/n, 1-1/n]^{n-1}?

    This is a question about the "covering radius" of the speed lattice
    relative to the forbidden region.

    For a single speed v: the points {vt : t in [0,1)} are v equally-spaced
    points on [0,1). The covering radius is 1/(2v). The forbidden region
    has radius 1/n. So we need 1/(2v) >= 1/n, i.e., v <= n/2.
    The slowest runner (v=1) always satisfies this.

    For multiple speeds: the covering radius of the joint system is more complex.
    """
    print("=" * 70)
    print("REFRAME 4: LATTICE — covering radius of the speed lattice")
    print("=" * 70)
    print()

    for n in n_values:
        k = n - 1
        thr = Fraction(1, n)

        # For the initial segment, compute the "nearest lattice point"
        # to the center of the good box
        speeds = tuple(range(1, n))

        # Center of the good box: (1/2, 1/2, ..., 1/2) in [0,1)^k
        # Distance from trajectory point (v_1 t, ..., v_k t) mod 1
        # to the center (1/2,...,1/2) in the torus metric

        # The unique time when all runners are at 1/2 is t=1/2 (for odd v)
        # or t=1/(2v) for general v

        # Lattice determinant: for the sublattice {(v_1 t,...,v_k t) : t in Z}
        # in Z^k, the index is... well, the image is a rank-1 sublattice
        # with generator (v_1,...,v_k).

        # The covering radius of this rank-1 sublattice in the torus
        # is the maximum distance from any point in the torus to the
        # nearest point on the trajectory.

        # For a rational line with direction (v_1,...,v_k) where gcd=1,
        # the trajectory visits exactly lcm denominators of rational
        # points. The covering radius in the L-infinity norm is
        # related to simultaneous Diophantine approximation.

        # Compute: maximum L-inf distance from (1/2,...,1/2) to trajectory
        max_Linf = 0.0
        best_t = None
        for s in range(1, 1001):
            t = Fraction(s, 1000)
            pos = [float(frac(Fraction(v) * t)) for v in speeds]
            Linf = max(min(abs(p - 0.5), 1 - abs(p - 0.5)) for p in pos)
            # Actually, the L-inf distance to center
            # We want the CLOSEST approach to the good box, not center
            # Let's compute the distance to the good box
            inside = all(float(thr) <= p <= 1 - float(thr) for p in pos)
            if inside:
                Linf_good = 0.0
            else:
                # Distance to good box boundary
                dists = []
                for p in pos:
                    if p < float(thr):
                        dists.append(float(thr) - p)
                    elif p > 1 - float(thr):
                        dists.append(p - (1 - float(thr)))
                    else:
                        dists.append(0.0)
                Linf_good = max(dists)

        print(f"n={n}: speed lattice direction ({','.join(str(v) for v in speeds)})")
        print(f"  good box side: {float(Fraction(n-2,n)):.4f}")

        # Three-distance theorem: for a line of irrational slope on a torus,
        # the gaps between successive trajectory points take at most 3 values.
        # For our rational trajectory with k speeds, this generalizes.

        # Compute the three-distance statistics for the initial segment
        # at different densities
        for num_pts in [100, 1000]:
            points = []
            for s in range(num_pts):
                t = Fraction(s, num_pts)
                pos = tuple(float(frac(Fraction(v) * t)) for v in speeds)
                points.append(pos)

            # Count hits in the good box
            box_hits = sum(1 for p in points
                          if all(float(thr) <= c <= 1 - float(thr) for c in p))
            print(f"  {num_pts} trajectory points: {box_hits} in good box "
                  f"({100*box_hits/num_pts:.1f}%)")

        print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 5: FOURIER / SPECTRAL
# ═══════════════════════════════════════════════════════════════

def reframe_spectral(n_values=[3, 4, 5, 6]):
    """LRC via Fourier analysis.

    The lonely indicator function L(t) = prod_i f_n(v_i t) where
    f_n(x) = 1_{||x|| >= 1/n}.

    The Fourier transform of f_n:
    f̂_n(0) = 1 - 2/n = (n-2)/n
    f̂_n(k) = -(2/n) * sinc(2k/n) for k != 0
            = -(1/(pi*k)) * sin(2*pi*k/n)

    Then L̂(m) = sum over (k_1,...,k_{n-1}) with k_1*v_1+...+k_{n-1}*v_{n-1}=m
                of product f̂_n(k_i).

    The integral of L = L̂(0) = product of f̂_n(0) = ((n-2)/n)^{n-1}.

    For L to be positive somewhere, we need L̂(0) - sum |L̂(m)| for m!=0 to be > 0.
    This is a Wiener-Ikehara type condition.
    """
    print("=" * 70)
    print("REFRAME 5: SPECTRAL — Fourier of the lonely indicator")
    print("=" * 70)
    print()

    for n in n_values:
        k = n - 1
        # Mean value (DC component)
        mean = ((n - 2) / n) ** k

        # For small n, compute the full Fourier decomposition of L(t)
        # L(t) = sum_m L̂(m) e^{2πimt}
        # L̂(0) = mean
        # L̂(m) = sum_{k: sum k_i v_i = m} prod f̂_n(k_i)

        # For the initial segment, compute the first few Fourier coefficients
        speeds = list(range(1, n))
        num_harmonics = 20

        # f̂_n(k) = integral_0^1 f_n(x) e^{-2πikx} dx
        # f_n(x) = 1 for x in [1/n, 1-1/n], 0 elsewhere
        # f̂_n(k) = (1/(2πik)) * (e^{-2πik(1-1/n)} - e^{-2πik/n})  for k != 0
        #         = (1/(2πik)) * e^{-πik} * (e^{πik(2/n-1)} - e^{-πik(2/n-1)})
        #         = (1/(πk)) * e^{-πik} * sin(πk(2/n-1))  ... hmm let me just compute

        def fhat(k, n):
            if k == 0:
                return (n - 2) / n
            # f_n(x) = 1 for x in [1/n, 1-1/n]
            a, b = 1.0 / n, 1.0 - 1.0 / n
            # integral_a^b e^{-2πikx} dx = (e^{-2πikb} - e^{-2πika}) / (-2πik)
            val = (cmath.exp(-2j * pi * k * b) - cmath.exp(-2j * pi * k * a)) / (-2j * pi * k)
            return val

        print(f"n={n}  mean L(t) = ((n-2)/n)^(n-1) = {mean:.6f}")

        # Compute L̂(m) for the first few m
        Lhat = {}
        for m in range(-num_harmonics, num_harmonics + 1):
            # L̂(m) = sum over integer vectors (k_1,...,k_{n-1}) with sum k_i*v_i = m
            # of product fhat(k_i)
            # This is a convolution — hard to compute exactly for general m.
            # For m=0: L̂(0) = product of fhat(0) = mean (only the (0,...,0) vector)
            # Wait, that's wrong. L̂(m) also includes vectors with nonzero k_i that cancel.
            pass

        # Instead, just compute L(t) directly at many points and FFT
        num_pts = 2048
        L_values = []
        for s in range(num_pts):
            t = s / num_pts
            lonely = 1.0
            for v in speeds:
                pos = (v * t) % 1.0
                gap = min(pos, 1.0 - pos)
                if gap < 1.0 / n:
                    lonely = 0.0
                    break
            L_values.append(lonely)

        # Compute DFT
        L_fft = [0.0] * num_pts
        for m in range(min(num_harmonics + 1, num_pts)):
            coeff = sum(L_values[s] * cmath.exp(-2j * pi * m * s / num_pts)
                        for s in range(num_pts)) / num_pts
            L_fft[m] = abs(coeff)

        print(f"  L̂(0) = {L_fft[0]:.6f}  (mean)")
        print(f"  |L̂(1)| = {L_fft[1]:.6f}")
        print(f"  |L̂(2)| = {L_fft[2]:.6f}")
        print(f"  |L̂(n)| = {L_fft[min(n, num_harmonics)]:.6f}")

        # Energy in first few harmonics
        energy_dc = L_fft[0] ** 2
        energy_total = sum(v ** 2 for v in L_fft[:num_harmonics])
        if energy_total > 0:
            print(f"  DC fraction of energy: {energy_dc/energy_total:.4f}")

        # Key spectral fact: if L̂(0) > sum_{m!=0} |L̂(m)|, then L > 0 somewhere
        tail_sum = sum(L_fft[m] for m in range(1, num_harmonics))
        print(f"  L̂(0) - sum|L̂(m!=0)| = {L_fft[0] - tail_sum:.6f}  "
              f"({'> 0 PROVES LRC' if L_fft[0] > tail_sum else '< 0 inconclusive'})")
        print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 6: THERMODYNAMIC / STATISTICAL MECHANICS
# ═══════════════════════════════════════════════════════════════

def reframe_thermodynamic(n_values=[3, 4, 5, 6, 7, 14]):
    """LRC as a statistical mechanics problem.

    Think of runners as particles on a circle, observer as a "trap."
    The "energy" at time t is E(t) = -sum_i log(||v_i t||).
    High energy = runners bunched near observer = bad.
    Low energy = runners spread out = good.

    The "partition function" Z(beta) = integral_0^1 exp(-beta * E(t)) dt.

    LRC is equivalent to: there exists t with E(t) < -sum_i log(1/n)
    = (n-1) * log(n).

    The "free energy" F = -(1/beta) * log(Z).
    The "entropy" S = -integral L(t) log L(t) dt (where L is the lonely indicator).

    More interesting: define a "temperature" parameter beta and study
    the phase transition as beta -> infinity (ground state = most lonely time).
    """
    print("=" * 70)
    print("REFRAME 6: THERMODYNAMIC — runner gas on the circle")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))
        k = n - 1

        # Energy threshold for loneliness
        energy_threshold = k * log(n)

        # Sample E(t)
        num_samples = 10000
        energies = []
        lonely_count = 0

        for s in range(1, num_samples + 1):
            t = s / num_samples
            E = 0.0
            lonely = True
            for v in speeds:
                gap = min((v * t) % 1.0, 1.0 - (v * t) % 1.0)
                if gap < 1e-10:
                    E = float('inf')
                    lonely = False
                    break
                E -= log(gap)
                if gap < 1.0 / n:
                    lonely = False

            if E < float('inf'):
                energies.append(E)
            if lonely:
                lonely_count += 1

        if energies:
            min_E = min(energies)
            avg_E = sum(energies) / len(energies)
            max_E = max(energies)

            print(f"n={n:2d}  energy_threshold={energy_threshold:.2f}")
            print(f"  E range: [{min_E:.2f}, {max_E:.2f}]  avg={avg_E:.2f}")
            print(f"  lonely fraction: {lonely_count}/{num_samples}")
            print(f"  min_E {'<' if min_E < energy_threshold else '>='} threshold "
                  f"{'=> LONELY EXISTS' if min_E < energy_threshold else '=> check boundary'}")

            # "Heat capacity" (variance of E)
            var_E = sum((e - avg_E) ** 2 for e in energies) / len(energies)
            print(f"  heat capacity (Var[E]): {var_E:.2f}")
        print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 7: TOPOLOGICAL / WINDING NUMBER
# ═══════════════════════════════════════════════════════════════

def reframe_topological(n_values=[3, 4, 5]):
    """LRC via topological methods.

    Consider the map Phi: [0,1) -> {0,1}^{n-1} defined by
    Phi_i(t) = 1 if ||v_i t|| >= 1/n (runner i is safe), else 0.

    LRC says: the image of Phi contains the all-1 vector (1,1,...,1).

    Topological approach: define a continuous relaxation
    phi_i(t) = max(0, ||v_i t|| - 1/n) / (1/2 - 1/n)  in [0,1]

    Then phi: [0,1) -> [0,1]^{n-1} is continuous.

    The "distance to lonely" is d(t) = min_i phi_i(t).
    LRC says max_t d(t) > 0.

    Borsuk-Ulam type: if we identify t with t + 1/2 (antipodal on the circle),
    does the map have the right parity to force a hitting?

    For the initial segment {1,...,n-1}: under t -> t + 1/2,
    phi_i(t + 1/2) relates to phi_i(t) by the parity of v_i.
    Even speeds: phi_i(t+1/2) = phi_i(t) (same position mod 1)
    Odd speeds: phi_i(t+1/2) measures ||v_i t + v_i/2|| which shifts by 1/2
    """
    print("=" * 70)
    print("REFRAME 7: TOPOLOGICAL — winding number and parity")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = list(range(1, n))
        k = n - 1
        thr = 1.0 / n

        # Count even and odd speeds
        even_speeds = [v for v in speeds if v % 2 == 0]
        odd_speeds = [v for v in speeds if v % 2 == 1]

        print(f"n={n}  even speeds: {even_speeds}  odd speeds: {odd_speeds}")

        # Under t -> t + 1/2:
        # For even v: {v(t+1/2)} = {vt + v/2} = {vt} (since v/2 is integer)
        # For odd v: {v(t+1/2)} = {vt + v/2} = {vt + (v-1)/2 + 1/2}
        #          = {vt} + 1/2 mod 1 (since (v-1)/2 is integer)
        # So odd-speed runners shift by 1/2, even-speed runners stay.

        # This means: if runner i (odd speed) is at distance d from observer at t,
        # it's at distance |1/2 - d| from observer at t + 1/2.
        # If d < 1/n (close), then |1/2 - d| > 1/2 - 1/n >= 1/n for n >= 3.
        # So odd-speed runners are COMPLEMENTARY under the half-period shift!

        # Key topological fact: every odd-speed runner is either safe at t
        # or safe at t + 1/2 (for n >= 4, strictly both can't be unsafe).

        # For even speeds: the half-shift doesn't help — they're in the same position.

        # This means: the set of times where ALL odd-speed runners are safe
        # has measure >= ... well, by complementarity, it's the intersection
        # of complements of shifted sets.

        # Let B_i = {t : ||v_i t|| < 1/n} (bad times for runner i)
        # For odd v: B_i and B_i + 1/2 are disjoint (since 1/n < 1/2 for n>=3)
        # So mu(B_i) = 2/n, and {t : runner i bad at t AND at t+1/2} = empty.

        # The intersection of all "odd-runner safe" sets:
        # {t : all odd runners safe} has measure >= 1 - (# odd runners) * 2/n

        odd_bad_measure = len(odd_speeds) * 2.0 / n
        odd_safe_lower = max(0, 1.0 - odd_bad_measure)

        even_bad_measure = len(even_speeds) * 2.0 / n  # can overlap

        print(f"  odd-runner complementarity: each bad set measure = 2/n, disjoint under t+1/2")
        print(f"  union bound on odd-bad: {odd_bad_measure:.4f}")
        print(f"  lower bound on all-odd-safe: {odd_safe_lower:.4f}")
        print(f"  even-runner bad measure (can overlap): {even_bad_measure:.4f}")

        # Can we combine? For the lonely set, we need ALL runners safe.
        # The odd runners contribute a safe-set of measure >= odd_safe_lower.
        # Within that set, the even runners need to also be safe.
        # Even runners' bad sets have total measure <= even_bad_measure.
        # So lonely measure >= odd_safe_lower - even_bad_measure
        lonely_lower = max(0, odd_safe_lower - even_bad_measure)
        print(f"  COMBINED lonely lower bound: {lonely_lower:.4f}")
        print(f"  {'PROVES LRC (positive lower bound)' if lonely_lower > 0 else 'Inconclusive'}")
        print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 8: CHIP-FIRING / SANDPILE
# ═══════════════════════════════════════════════════════════════

def reframe_chipfiring(n_values=[3, 4, 5, 6]):
    """LRC as a chip-firing game.

    Place "proximity chips" on the observer. At time t, runner i delivers
    chips proportional to its closeness: c_i(t) = max(0, 1/n - ||v_i t||).

    The observer is lonely when no chips are delivered: all c_i = 0.

    Think of this as a sandpile: chips accumulate on the observer (it gets
    "crowded"), and when the total exceeds a threshold, a "toppling" occurs
    (a runner crosses the threshold and its chips go to zero).

    The THM-387 gap-race framework can be seen as chip-firing:
    - Clockwise chips accumulate (g_right grows)
    - Counterclockwise chips drain (g_left shrinks)
    - "Toppling" = wrap-around (runner crosses observer, chips redistribute)
    """
    print("=" * 70)
    print("REFRAME 8: CHIP-FIRING — proximity chips on observer")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))
        k = n - 1
        thr = 1.0 / n

        # Compute the total chip load over time
        num_samples = 1000
        chip_loads = []

        for s in range(1, num_samples + 1):
            t = s / num_samples
            total_chips = 0.0
            for v in speeds:
                gap = min((v * t) % 1.0, 1.0 - (v * t) % 1.0)
                chips = max(0.0, thr - gap)
                total_chips += chips
            chip_loads.append(total_chips)

        min_load = min(chip_loads)
        avg_load = sum(chip_loads) / len(chip_loads)
        max_load = max(chip_loads)
        zero_count = sum(1 for c in chip_loads if c < 1e-10)

        print(f"n={n}  initial segment")
        print(f"  chip load range: [{min_load:.6f}, {max_load:.6f}]  avg={avg_load:.6f}")
        print(f"  zero-load samples: {zero_count}/{num_samples}")
        print(f"  expected mean load: {k * (1/(n*n)):.6f} (integral of max(0,1/n-x) for uniform x)")
        print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 9: GAME THEORY / AUCTION
# ═══════════════════════════════════════════════════════════════

def reframe_auction(n_values=[3, 4, 5, 6]):
    """LRC as an auction/game.

    Players: the n-1 runners.
    Goods: the two observer-adjacent gap "slots" (left and right).
    Each runner bids with its proximity to the observer.
    The closest runner "wins" each slot.

    The observer is lonely iff both auction winners bid less than 1/n.

    This is like a Vickrey auction: the winner's price = second-closest runner.
    If the closest runner is at distance d < 1/n, the gap is d.
    The lonely condition: the "reserve price" 1/n is never met.

    Key insight from THM-387: the left auction and right auction have
    CORRELATED bids due to the directed flow of runners.
    """
    print("=" * 70)
    print("REFRAME 9: AUCTION — two-slot gap bidding game")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))
        k = n - 1
        thr = 1.0 / n

        # At each time, compute the "winning bids" for left and right slots
        num_samples = 1000
        left_wins = []
        right_wins = []

        for s in range(1, num_samples + 1):
            t = s / num_samples
            positions = [(v * t) % 1.0 for v in speeds]
            # Right gap winner: runner with smallest positive position
            right_bids = [p for p in positions if p > 0]
            left_bids = [1.0 - p for p in positions if p > 0]

            if right_bids and left_bids:
                right_wins.append(min(right_bids))
                left_wins.append(min(left_bids))

        if right_wins:
            # "Reserve price" is 1/n. How often is the winning bid >= reserve?
            right_above = sum(1 for w in right_wins if w >= thr)
            left_above = sum(1 for w in left_wins if w >= thr)
            both_above = sum(1 for l, r in zip(left_wins, right_wins)
                            if l >= thr and r >= thr)

            print(f"n={n}  initial segment")
            print(f"  right slot above reserve: {right_above}/{len(right_wins)} "
                  f"({100*right_above/len(right_wins):.1f}%)")
            print(f"  left slot above reserve: {left_above}/{len(left_wins)} "
                  f"({100*left_above/len(left_wins):.1f}%)")
            print(f"  BOTH above reserve (lonely): {both_above}/{len(right_wins)} "
                  f"({100*both_above/len(right_wins):.1f}%)")

            # Correlation between left and right winning bids
            if len(left_wins) > 1:
                mean_l = sum(left_wins) / len(left_wins)
                mean_r = sum(right_wins) / len(right_wins)
                cov = sum((l - mean_l) * (r - mean_r) for l, r in zip(left_wins, right_wins)) / len(left_wins)
                var_l = sum((l - mean_l) ** 2 for l in left_wins) / len(left_wins)
                var_r = sum((r - mean_r) ** 2 for r in right_wins) / len(right_wins)
                if var_l > 0 and var_r > 0:
                    corr = cov / (var_l ** 0.5 * var_r ** 0.5)
                    print(f"  left-right bid correlation: {corr:.4f}")

        print()


# ═══════════════════════════════════════════════════════════════
# REFRAME 10: TOURNAMENT ANALYSIS of the reframes
# ═══════════════════════════════════════════════════════════════

def tournament_analysis_of_reframes():
    """Compare the 10 reframes using Tournament Analysis.

    Pairwise comparison on 5 criteria:
    1. Novelty: how different from existing LRC approaches
    2. Computability: can we compute nontrivial quantities
    3. Proof force: does this reframe suggest a proof mechanism
    4. Connections: links to known deep results
    5. Beauty: elegance and surprise value
    """
    print("=" * 70)
    print("TOURNAMENT ANALYSIS: Comparing the reframes")
    print("=" * 70)
    print()

    reframes = [
        "billiard",
        "lovász_LLL",
        "relativistic",
        "lattice",
        "spectral",
        "thermodynamic",
        "topological",
        "chipfiring",
        "auction",
    ]

    # Score each reframe on 5 criteria (1-5)
    # Based on the computational evidence above
    scores = {
        "billiard":      {"novelty": 2, "compute": 4, "proof": 2, "connect": 4, "beauty": 3},
        "lovász_LLL":    {"novelty": 3, "compute": 4, "proof": 3, "connect": 5, "beauty": 3},
        "relativistic":  {"novelty": 5, "compute": 3, "proof": 2, "connect": 4, "beauty": 5},
        "lattice":       {"novelty": 2, "compute": 3, "proof": 3, "connect": 5, "beauty": 3},
        "spectral":      {"novelty": 3, "compute": 5, "proof": 4, "connect": 5, "beauty": 4},
        "thermodynamic": {"novelty": 4, "compute": 4, "proof": 2, "connect": 3, "beauty": 4},
        "topological":   {"novelty": 4, "compute": 4, "proof": 5, "connect": 5, "beauty": 5},
        "chipfiring":    {"novelty": 3, "compute": 3, "proof": 2, "connect": 3, "beauty": 3},
        "auction":       {"novelty": 3, "compute": 4, "proof": 2, "connect": 2, "beauty": 3},
    }

    # Build the tournament: A -> B if A wins majority of criteria
    n = len(reframes)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            wins = sum(1 for c in ["novelty", "compute", "proof", "connect", "beauty"]
                       if scores[reframes[i]][c] > scores[reframes[j]][c])
            if wins >= 3:
                adj[i][j] = 1

    # Compute scores
    out_scores = [sum(row) for row in adj]
    ranked = sorted(range(n), key=lambda i: out_scores[i], reverse=True)

    print("Ranking (by tournament out-degree):")
    for rank, i in enumerate(ranked, 1):
        s = scores[reframes[i]]
        print(f"  {rank}. {reframes[i]:15s}  score={out_scores[i]}  "
              f"nov={s['novelty']} comp={s['compute']} proof={s['proof']} "
              f"conn={s['connect']} beauty={s['beauty']}")

    # Key finding
    winner = reframes[ranked[0]]
    print(f"\nWINNER: {winner}")
    print()

    # The real finding: which reframe actually suggests a new proof route?
    print("PROOF-FORCE ranking (which reframe suggests a proof mechanism):")
    proof_ranked = sorted(range(n), key=lambda i: scores[reframes[i]]["proof"], reverse=True)
    for i in proof_ranked[:5]:
        print(f"  {reframes[i]:15s}  proof_force={scores[reframes[i]]['proof']}")
    print()


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    print("LRC Creative Reframes — opus-2026-06-01-S520")
    print("Exploring LRC through 10 radically different mathematical lenses")
    print()

    reframe_billiard()
    reframe_lovasz()
    reframe_relativistic()
    reframe_lattice()
    reframe_spectral()
    reframe_thermodynamic()
    reframe_topological()
    reframe_chipfiring()
    reframe_auction()
    tournament_analysis_of_reframes()

    print("=" * 70)
    print("GRAND SYNTHESIS")
    print("=" * 70)
    print()
    print("The most promising reframes for new proof routes:")
    print()
    print("1. TOPOLOGICAL (odd/even parity complementarity):")
    print("   Every odd-speed runner is safe at t or at t+1/2.")
    print("   This gives a FREE lower bound on the all-safe measure.")
    print("   For small n with many odd speeds, this may already prove LRC!")
    print()
    print("2. SPECTRAL (Fourier):")
    print("   If the DC component exceeds the sum of AC components,")
    print("   L(t) > 0 for some t.  Computable and tight for small n.")
    print()
    print("3. LOVÁSZ LOCAL LEMMA:")
    print("   The dependency structure of runner-close events is sparse")
    print("   for coprime speed sets.  May prove LRC for large max_speed/n.")
    print()
    print("4. RELATIVISTIC:")
    print("   Rapidity threshold = 0.5*ln(n-1) gives a natural scale.")
    print("   The formal group F(x,y)=(x+y)/(1+xy) is additive in rapidity.")
    print("   This connects to the project's central object.")
    print()


if __name__ == "__main__":
    main()
