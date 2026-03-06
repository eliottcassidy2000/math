"""
Investigate the parity dependence of blueself existence.

Observed pattern:
  n=4 (even): 1 blueself class, 0 blackself
  n=5 (odd):  0 blueself classes, 2 blackself
  n=6 (even): 2 blueself classes, 6 blackself
  n=7 (odd):  0 blueself classes, 30 blackself

CONJECTURE: Blueself classes exist ONLY at even n.

This script investigates WHY via score-sequence analysis of the flip operation.

Author: opus-2026-03-06-S16
"""

import itertools
from collections import defaultdict


def build_tiles(n):
    tiles = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            tiles.append((x, y))
    return tiles


def build_transpose_map(n, tiles):
    tile_idx = {(x, y): i for i, (x, y) in enumerate(tiles)}
    return [tile_idx[(n - y + 1, n - x + 1)] for i, (x, y) in enumerate(tiles)]


def bits_to_adj(n, tiles, bits):
    verts = list(range(n, 0, -1))
    vert_idx = {v: i for i, v in enumerate(verts)}
    A = [[0] * n for _ in range(n)]
    for k in range(n - 1):
        A[k][k + 1] = 1
    for i, (x, y) in enumerate(tiles):
        xi, yi = vert_idx[x], vert_idx[y]
        if bits[i] == 0:
            A[xi][yi] = 1
        else:
            A[yi][xi] = 1
    return A


def is_grid_symmetric(bits, trans_map):
    for i in range(len(bits)):
        if trans_map[i] != i and bits[i] != bits[trans_map[i]]:
            return False
    return True


def flip_bits(bits):
    return tuple(1 - b for b in bits)


def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)], reverse=True))


def score_vec(A, n):
    """Score vector in PATH ORDER (not sorted)."""
    return tuple(sum(A[i]) for i in range(n))


def analyze_flip_scores(n):
    """Analyze how flip changes scores, focusing on grid-symmetric tilings."""
    print(f"\n{'=' * 60}")
    print(f"  FLIP SCORE ANALYSIS at n={n} ({'even' if n % 2 == 0 else 'odd'})")
    print(f"{'=' * 60}")

    tiles = build_tiles(n)
    trans_map = build_transpose_map(n, tiles)
    m = len(tiles)
    total = 1 << m

    # For a tournament with fixed path n->n-1->...->1:
    # Vertex at position i (0-indexed) has vertex label n-i.
    # Position 0 = vertex n (first in path, beats position 1)
    # Position n-1 = vertex 1 (last in path, beaten by position n-2)
    #
    # Path out-edges: position i -> position i+1 (for i < n-1)
    # Path in-edges: position i-1 -> position i (for i > 0)
    #
    # Non-path neighbors of position i:
    #   All j with |i-j| >= 2
    #   Count: n - 1 - |{path neighbors}|
    #     Interior (0<i<n-1): n-3 non-path neighbors
    #     Endpoint i=0: n-2 non-path neighbors
    #     Endpoint i=n-1: n-2 non-path neighbors
    #
    # Score of position i in T: path_out + k_i
    #   where k_i = non-path out-degree, path_out = 1 if i<n-1 else 0
    #
    # Score of position i in flip(T): path_out + (non_path_count - k_i)
    #   Interior: 1 + (n-3-k_i) = n-2-k_i
    #   i=0: 1 + (n-2-k_0) = n-1-k_0
    #   i=n-1: 0 + (n-2-k_{n-1}) = n-2-k_{n-1}

    print(f"\n  Score transformation under flip:")
    print(f"    Interior position (0<i<n-1): s_i -> n-2-k_i = (n-1) - s_i")
    print(f"    Position 0 (vertex n):       s_0 -> n - s_0")
    print(f"    Position n-1 (vertex 1):     s_{n-1} -> (n-2) - s_{n-1}")
    print(f"")
    print(f"    Interior: s + s' = n-1 = {n-1}")
    print(f"    Pos 0:    s + s' = n = {n}")
    print(f"    Pos n-1:  s + s' = n-2 = {n-2}")

    # For sorted score sequences to match:
    # The multiset {s_0, s_1, ..., s_{n-1}} must equal
    # {n-s_0, (n-1)-s_1, ..., (n-1)-s_{n-2}, (n-2)-s_{n-1}}

    # At odd n, let's check: is this even POSSIBLE for grid-symmetric tilings?
    gs_count = 0
    gs_same_sorted_scores = 0
    gs_same_scores_but_not_iso = 0

    for mask in range(total):
        bits = tuple((mask >> k) & 1 for k in range(m))
        if not is_grid_symmetric(bits, trans_map):
            continue
        gs_count += 1

        A = bits_to_adj(n, tiles, bits)
        flip_b = flip_bits(bits)
        A_f = bits_to_adj(n, tiles, flip_b)

        ss = score_seq(A, n)
        ss_f = score_seq(A_f, n)

        if ss == ss_f:
            gs_same_sorted_scores += 1

    print(f"\n  Grid-symmetric tilings: {gs_count}")
    print(f"  Grid-sym with same sorted scores as flip: {gs_same_sorted_scores}")
    print(f"  (This is NECESSARY for flip to be isomorphic)")

    if gs_same_sorted_scores == 0 and n % 2 == 1:
        print(f"\n  *** At odd n={n}, NO grid-symmetric tiling has the same score")
        print(f"  sequence as its flip! This PROVES no blueself exists at n={n}. ***")

    # Let's also look at the score vectors (not sorted) for small cases
    if gs_count <= 100:
        print(f"\n  Score vector pairs (T, flip(T)) for grid-symmetric tilings:")
        for mask in range(total):
            bits = tuple((mask >> k) & 1 for k in range(m))
            if not is_grid_symmetric(bits, trans_map):
                continue
            A = bits_to_adj(n, tiles, bits)
            flip_b = flip_bits(bits)
            A_f = bits_to_adj(n, tiles, flip_b)
            sv = score_vec(A, n)
            sv_f = score_vec(A_f, n)
            ss = score_seq(A, n)
            ss_f = score_seq(A_f, n)
            match = "MATCH" if ss == ss_f else ""
            print(f"    mask={mask:>5d}: T_scores={sv} -> flip_scores={sv_f}"
                  f"  sorted: {ss} vs {ss_f} {match}")

    # ANALYSIS: At odd n, why can't scores match?
    print(f"\n  PARITY ANALYSIS:")
    print(f"    n = {n} ({'odd' if n%2 else 'even'})")
    print(f"    Endpoint sum s_0+s'_0 = {n} ({'even' if n%2==0 else 'odd'})")
    print(f"    Endpoint sum s_{n-1}+s'_{n-1} = {n-2} ({'even' if (n-2)%2==0 else 'odd'})")
    print(f"    Interior sum s_i+s'_i = {n-1} ({'even' if (n-1)%2==0 else 'odd'})")
    print(f"")
    print(f"    Total score of T = C(n,2) = {n*(n-1)//2}")
    print(f"    Total score of flip(T) = C(n,2) = {n*(n-1)//2}")
    print(f"    Sum of ALL scores in T + flip(T):")
    S_sum = n + (n - 2) + (n - 2) * (n - 1)
    print(f"      = n + (n-2) + (n-2)*(n-1) = {S_sum}")
    print(f"      = 2*C(n,2) = {n*(n-1)} -> check: {S_sum == n*(n-1)}")

    # KEY insight: the PARITY of the sum s_0 + s'_0 = n
    # For T and flip(T) to have the same sorted score multiset,
    # each score s in the multiset maps to some s' = (constant) - s.
    # The constants are n-1 (interior), n (pos 0), n-2 (pos n-1).
    # At odd n: the three constants have parities odd, odd, odd.
    # This means s and s' have opposite parities for each position!
    # If s is even, s' is odd, and vice versa.
    # So the PARITY SIGNATURE (even/odd pattern of sorted scores) FLIPS.
    # For the multisets to match, we'd need the same number of even and odd scores.
    # Score parity: each position contributes s_i to T and (constant-s_i) to flip(T).
    # At odd n, constant is always odd, so s_i and constant-s_i have opposite parities.
    # The multiset of scores has some number of even scores and some odd scores.
    # After flip, even scores become odd and vice versa.
    # For the sorted multisets to match, we need the number of even scores = number of odd scores.
    # Since n is odd, we can't split n scores evenly into two groups!

    even_count = n // 2
    odd_count = (n + 1) // 2
    print(f"\n    At n={n} ({'odd' if n%2 else 'even'}):")
    print(f"    Score set of a tournament is a subset of {{0,1,...,{n-1}}} with sum {n*(n-1)//2}")
    print(f"    Under flip, each score s -> (const - s) where const is always {'odd' if n%2 else 'mixed'}.")

    if n % 2 == 1:
        print(f"\n    *** PROOF for odd n: ***")
        print(f"    All three flip-constants (n-1={n-1}, n={n}, n-2={n-2}) are EVEN.")
        # Wait, let me recheck. n=5: n-1=4 (even), n=5 (odd), n-2=3 (odd).
        # n=7: n-1=6 (even), n=7 (odd), n-2=5 (odd).
        print(f"    Actually: n-1={n-1} ({'even' if (n-1)%2==0 else 'odd'}), "
              f"n={n} ({'even' if n%2==0 else 'odd'}), "
              f"n-2={n-2} ({'even' if (n-2)%2==0 else 'odd'})")

        # At odd n: n-1 is even, n is odd, n-2 is odd.
        # Interior positions: s_i -> (n-1) - s_i. Since n-1 is even, s_i and s'_i have SAME parity.
        # Position 0: s_0 -> n - s_0. Since n is odd, s_0 and s'_0 have DIFFERENT parity.
        # Position n-1: s_{n-1} -> (n-2) - s_{n-1}. Since n-2 is odd, they have DIFFERENT parity.
        #
        # So the two endpoint scores FLIP parity, while interior scores KEEP parity.
        # The sorted score multiset changes: endpoint scores switch even↔odd.
        # For the multisets to match, the two "flipped" endpoint values must
        # map to values that already existed in the multiset.

        print(f"    Interior: n-1={n-1} is even -> s_i parity PRESERVED by flip")
        print(f"    Position 0: n={n} is odd -> s_0 parity FLIPPED by flip")
        print(f"    Position n-1: n-2={n-2} is odd -> s_{{n-1}} parity FLIPPED by flip")
        print(f"")
        print(f"    The two endpoint scores change parity under flip.")
        print(f"    For score multisets to match: need specific structure.")
        print(f"    This is highly constraining but not immediately impossible.")

    # Let me check: for grid-symmetric tilings specifically,
    # what is the score structure?
    if gs_count <= 200:
        print(f"\n  Detailed score analysis for grid-symmetric tilings:")
        score_parities = defaultdict(int)
        for mask in range(total):
            bits = tuple((mask >> k) & 1 for k in range(m))
            if not is_grid_symmetric(bits, trans_map):
                continue
            A = bits_to_adj(n, tiles, bits)
            sv = score_vec(A, n)
            # Parity pattern
            parity = tuple(s % 2 for s in sv)
            score_parities[parity] += 1

        print(f"    Score parity patterns (0=even, 1=odd) in path order:")
        for parity, count in sorted(score_parities.items()):
            # What would flip do to this pattern?
            flip_parity = list(parity)
            # Interior: same parity (n-1 even at odd n)
            # Position 0: flip (n odd at odd n)
            # Position n-1: flip (n-2 odd at odd n)
            if n % 2 == 1:
                flip_parity[0] = 1 - flip_parity[0]
                flip_parity[n-1] = 1 - flip_parity[n-1]
            flip_parity = tuple(flip_parity)

            # For score multisets to match, sorted parities must match
            sorted_p = tuple(sorted(parity))
            sorted_fp = tuple(sorted(flip_parity))
            can_match = sorted_p == sorted_fp

            print(f"      {parity} (count={count}) -> flip parity {flip_parity}"
                  f"  sorted: {sorted_p} vs {sorted_fp}"
                  f"  {'CAN match' if can_match else 'CANNOT match'}")


if __name__ == '__main__':
    for n in [4, 5, 6, 7]:
        analyze_flip_scores(n)
