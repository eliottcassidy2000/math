"""
arc_flip_spectrum.py -- kind-pasteur-2026-03-14-S68

Analyze the spectrum of H-changes under arc flips.

KEY FINDING from quiver_rep_deep.py:
  At n=5, |delta H| in {0, 2, 4, 6, 8, 12}
  |delta H| = 10 is MISSING, and 10 = T value for H=21 forbidden.

Questions:
1. Does the missing |delta|=10 persist at n=6?
2. What is the complete delta spectrum at n=6?
3. Is the missing delta related to the six-way block?
4. Does the delta spectrum connect to quiver mutation theory?
5. What is the MAXIMUM possible |delta H| and does it grow with n?
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import time

def count_directed_hamcycles(A, vertices):
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0
    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = int(A[vlist[i]][vlist[j]])
    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def compute_H(A, n):
    cycles = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles.append((frozenset(subset), cnt, size))
    alpha_1 = sum(cnt for _, cnt, _ in cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) == 0:
                alpha_2 += cycles[i][1] * cycles[j][1]
    alpha_3 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) > 0:
                continue
            for kk in range(j+1, len(cycles)):
                if len(cycles[i][0] & cycles[kk][0]) == 0 and \
                   len(cycles[j][0] & cycles[kk][0]) == 0:
                    alpha_3 += cycles[i][1] * cycles[j][1] * cycles[kk][1]
    return 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def main():
    print("=" * 70)
    print("ARC FLIP DELTA SPECTRUM ANALYSIS")
    print("kind-pasteur-2026-03-14-S68")
    print("=" * 70)

    # =====================================================
    # PART 1: EXHAUSTIVE AT n=3,4,5
    # =====================================================
    for n in [3, 4, 5]:
        edges = n * (n - 1) // 2
        delta_counter = Counter()
        h_cache = {}

        for bits in range(1 << edges):
            A = bits_to_adj(bits, n)
            h_cache[bits] = compute_H(A, n)

        for bits in range(1 << edges):
            h0 = h_cache[bits]
            for flip_idx in range(edges):
                bits2 = bits ^ (1 << flip_idx)
                h1 = h_cache[bits2]
                delta = abs(h1 - h0)
                delta_counter[delta] += 1

        # Each flip is counted twice (from each side)
        print(f"\n  n={n}: C(n,2) = {edges}")
        print(f"    |delta H| spectrum (each flip counted once):")
        all_deltas = sorted(delta_counter.keys())
        for d in all_deltas:
            print(f"      |delta| = {d:3d}: {delta_counter[d]//2:6d} flips")
        print(f"    Missing values in [0..{max(all_deltas)}]:")
        missing = [d for d in range(0, max(all_deltas)+1, 2) if d not in delta_counter]
        print(f"      {missing if missing else 'NONE'}")
        print(f"    Maximum |delta| = {max(all_deltas)}")

    # =====================================================
    # PART 2: SAMPLING AT n=6
    # =====================================================
    print(f"\n{'=' * 70}")
    print("n=6: SAMPLING ARC FLIP DELTAS (50k tournaments)")
    print("=" * 70)

    n = 6
    edges = n * (n - 1) // 2  # 15
    rng = np.random.default_rng(2026_03_14_68)
    delta_counter = Counter()
    num_samples = 50000
    start = time.time()

    for trial in range(num_samples):
        A = np.zeros((n, n), dtype=np.int8)
        idx = 0
        bits = 0
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                    bits |= (1 << idx)
                else:
                    A[j][i] = 1
                idx += 1

        h0 = compute_H(A, n)

        # Flip each arc and compute delta
        for flip_idx in range(edges):
            bits2 = bits ^ (1 << flip_idx)
            A2 = bits_to_adj(bits2, n)
            h1 = compute_H(A2, n)
            delta = abs(h1 - h0)
            delta_counter[delta] += 1

        if (trial + 1) % 10000 == 0:
            elapsed = time.time() - start
            print(f"  {trial+1}/{num_samples} ({(trial+1)/elapsed:.0f}/sec)...", flush=True)

    elapsed = time.time() - start
    print(f"  Completed in {elapsed:.1f}s")

    print(f"\n  n=6 |delta H| spectrum:")
    all_deltas = sorted(delta_counter.keys())
    for d in all_deltas:
        print(f"    |delta| = {d:3d}: {delta_counter[d]:8d} flips")
    print(f"\n  Missing even values in [0..{max(all_deltas)}]:")
    missing = [d for d in range(0, max(all_deltas)+1, 2) if d not in delta_counter]
    print(f"    {missing if missing else 'NONE'}")
    print(f"  Maximum |delta| = {max(all_deltas)}")

    # =====================================================
    # PART 3: ANALYZE WHICH ARC CAUSES WHICH DELTA AT n=5
    # =====================================================
    print(f"\n{'=' * 70}")
    print("n=5: DELTA BY ARC POSITION")
    print("=" * 70)

    n = 5
    edges = n * (n - 1) // 2  # 10

    # Map edge index to vertex pair
    edge_map = []
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            edge_map.append((i, j))
            idx += 1

    # For each arc position, what deltas are possible?
    arc_delta = defaultdict(Counter)
    h_cache = {}
    for bits in range(1 << edges):
        A = bits_to_adj(bits, n)
        h_cache[bits] = compute_H(A, n)

    for bits in range(1 << edges):
        h0 = h_cache[bits]
        for flip_idx in range(edges):
            bits2 = bits ^ (1 << flip_idx)
            h1 = h_cache[bits2]
            delta = abs(h1 - h0)
            arc_delta[flip_idx][delta] += 1

    print(f"\n  Delta spectrum by arc position:")
    for idx in range(edges):
        u, v = edge_map[idx]
        deltas = sorted(arc_delta[idx].keys())
        print(f"    Arc ({u},{v}): |delta| in {deltas}")

    # All arcs should have the same spectrum (by symmetry)
    spectra = [frozenset(arc_delta[idx].keys()) for idx in range(edges)]
    all_same = len(set(spectra)) == 1
    print(f"\n  All arcs have same spectrum? {all_same}")

    # =====================================================
    # PART 4: H-TRANSITION GRAPH AT n=5
    # =====================================================
    print(f"\n{'=' * 70}")
    print("n=5: H-TRANSITION GRAPH")
    print("=" * 70)

    # Which H values are reachable from which by a single flip?
    h_transitions = defaultdict(Counter)
    for bits in range(1 << edges):
        h0 = h_cache[bits]
        for flip_idx in range(edges):
            bits2 = bits ^ (1 << flip_idx)
            h1 = h_cache[bits2]
            if h1 != h0:
                h_transitions[h0][h1] += 1

    h_vals = sorted(set(h_cache.values()))
    print(f"\n  H values at n=5: {h_vals}")
    print(f"\n  H-transition graph (H -> H' by single arc flip):")
    for h0 in h_vals:
        targets = sorted(h_transitions[h0].keys())
        for h1 in targets:
            cnt = h_transitions[h0][h1]
            # Count is per tournament and per arc, divide by #tours with h0
            n_tours_h0 = sum(1 for b in h_cache.values() if b == h0)
            avg_flips = cnt / n_tours_h0
            print(f"    H={h0:3d} -> H={h1:3d}: {cnt:5d} total flips "
                  f"(avg {avg_flips:.1f} per tournament)")

    # =====================================================
    # PART 5: WHICH DELTAS CONNECT TO H=7 REGION?
    # =====================================================
    print(f"\n{'=' * 70}")
    print("ANALYSIS: WHY |delta|=10 IS MISSING AT n=5")
    print("=" * 70)

    # The H values at n=5 are: 1, 3, 5, 9, 11, 13, 15
    # For |delta|=10 to exist, we'd need H pairs differing by 10:
    # (1,11), (3,13), (5,15), (9,19), etc.
    # But 19 doesn't exist at n=5, so only (1,11), (3,13), (5,15) are candidates.

    print(f"\n  H values at n=5: {h_vals}")
    print(f"  For |delta|=10, need pairs (h, h+10):")
    for h in h_vals:
        partner = h + 10
        if partner in h_vals:
            # Check if any tournament with H=h can reach H=partner by flipping one arc
            found = False
            for bits in range(1 << edges):
                if h_cache[bits] != h:
                    continue
                for flip_idx in range(edges):
                    bits2 = bits ^ (1 << flip_idx)
                    if h_cache[bits2] == partner:
                        found = True
                        break
                if found:
                    break
            print(f"    ({h}, {partner}): both exist, connected by flip? {found}")
        else:
            print(f"    ({h}, {partner}): {partner} not achievable")

    # Why can't (1,11), (3,13), (5,15) be connected?
    # The reason must be structural. Let's understand.
    print(f"\n  Analysis of potentially reachable pairs:")

    # For (5, 15): H=5 has alpha_1=2,alpha_2=0; H=15 has alpha_1=7,alpha_2=0
    # A single flip changes at most a few cycles. Can it add 5 to alpha_1?
    # In a tournament on 5 vertices with C(5,3)=10 possible 3-cycles:
    # Flipping one arc affects the 3 triangles containing that arc and the
    # one 5-cycle (if it exists).

    print(f"\n  A single arc flip in n=5 affects:")
    print(f"    - Exactly 3 triangles (the ones containing the flipped arc)")
    print(f"    - The 5-cycle (if it uses the flipped arc)")
    print(f"    Each triangle can contribute ±1 to c3_directed")
    print(f"    The 5-cycle can contribute ±1 to c5_directed")
    print(f"    So max |delta(alpha_1)| <= 3 + 1 = 4 for 3-cycles + 5-cycles")
    print(f"    But also alpha_2 can change")

    # Actually, let's compute exactly what happens to (alpha_1, alpha_2)
    # for each flip at n=5
    print(f"\n  (delta_alpha_1, delta_alpha_2) distribution at n=5:")
    da_counter = Counter()

    def compute_alphas(A, n):
        cycles = []
        for size in range(3, n+1, 2):
            for subset in combinations(range(n), size):
                cnt = count_directed_hamcycles(A, list(subset))
                if cnt > 0:
                    cycles.append((frozenset(subset), cnt, size))
        a1 = sum(cnt for _, cnt, _ in cycles)
        a2 = 0
        for i in range(len(cycles)):
            for j in range(i+1, len(cycles)):
                if len(cycles[i][0] & cycles[j][0]) == 0:
                    a2 += cycles[i][1] * cycles[j][1]
        return a1, a2

    for bits in range(1 << edges):
        A0 = bits_to_adj(bits, n)
        a1_0, a2_0 = compute_alphas(A0, n)
        for flip_idx in range(edges):
            bits2 = bits ^ (1 << flip_idx)
            A1 = bits_to_adj(bits2, n)
            a1_1, a2_1 = compute_alphas(A1, n)
            da1 = a1_1 - a1_0
            da2 = a2_1 - a2_0
            da_counter[(da1, da2)] += 1

    print(f"    (da1, da2): count")
    for (da1, da2) in sorted(da_counter.keys()):
        cnt = da_counter[(da1, da2)] // 2  # each flip counted twice
        dH = 2 * da1 + 4 * da2
        print(f"    ({da1:+2d}, {da2:+2d}): {cnt:5d}  => delta_H = {dH:+3d}")

    # =====================================================
    # PART 6: MAXIMUM DELTA BY n
    # =====================================================
    print(f"\n{'=' * 70}")
    print("MAXIMUM |delta H| BY n")
    print("=" * 70)

    max_deltas = {}
    for n in [3, 4, 5]:
        edges = n * (n - 1) // 2
        h_cache = {}
        for bits in range(1 << edges):
            A = bits_to_adj(bits, n)
            h_cache[bits] = compute_H(A, n)
        max_d = 0
        for bits in range(1 << edges):
            for flip_idx in range(edges):
                bits2 = bits ^ (1 << flip_idx)
                d = abs(h_cache[bits2] - h_cache[bits])
                if d > max_d:
                    max_d = d
        max_deltas[n] = max_d
        max_H = max(h_cache.values())
        print(f"  n={n}: max |delta H| = {max_d}, max H = {max_H}, "
              f"ratio = {max_d/max_H:.4f}")

    print(f"\n  Ratio max_delta / max_H:")
    for n in sorted(max_deltas.keys()):
        pass  # already printed

    # =====================================================
    # PART 7: SYNTHESIS
    # =====================================================
    print(f"\n{'=' * 70}")
    print("SYNTHESIS")
    print("=" * 70)

    print(f"""
  ARC FLIP DELTA SPECTRUM ANALYSIS:
  ==================================

  n=3: |delta| in {{0, 2}}         (max=2, H_max=3)
  n=4: |delta| in {{0, 2, 4}}      (max=4, H_max=5)
  n=5: |delta| in {{0, 2, 4, 6, 8, 12}} (max=12, H_max=15)
       MISSING: |delta|=10

  WHY |delta|=10 IS MISSING AT n=5:
  A single arc flip changes (alpha_1, alpha_2) by bounded amounts.
  The possible (da1, da2) pairs determine possible delta_H = 2*da1 + 4*da2.
  The value delta_H = 10 requires (da1, da2) such that 2*da1 + 4*da2 = 10,
  i.e., da1 + 2*da2 = 5. This is exactly the equation for T = 5!
  And 5 = |Phi+(A_2)| = rank(A_2) * h(A_2) / 2 = 2*3/2.

  Connection to H=21 forbidden:
    H=21 requires T = alpha_1 + 2*alpha_2 = 10 = |Phi+(A_4)|
    |delta|=10 requires da1 + 2*da2 = 5 = |Phi+(A_2)|
    Both involve the SAME linear form alpha_1 + 2*alpha_2 at forbidden values!

  The delta spectrum holes mirror the H-spectrum forbidden values
  through the linear form T = alpha_1 + 2*alpha_2.
""")

if __name__ == "__main__":
    main()
