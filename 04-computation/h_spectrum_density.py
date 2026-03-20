#!/usr/bin/env python3
"""h_spectrum_density.py -- Toward proving 7 and 21 are the only permanent gaps.

Session: kind-pasteur-2026-03-20-S3

APPROACH: Show that the achievable H-spectrum becomes arbitrarily dense
as n grows, with ONLY H=7 and H=21 as permanent gaps.

KEY IDEAS:
1. MONOTONICITY: If H=k is achievable at n, it's achievable at n+1
   (add an isolated source/sink vertex: H stays the same)
   Wait — can't add isolated vertex to tournament. Instead: add vertex
   that dominates all (source) or is dominated by all (sink).
   Source v: H(T+source) = H(T) (the source must be first in any HP)
   Actually that's wrong too. Adding a source: score(v)=n-1, and
   H(T+source) = H(T) + Delta. If v is a source, then v can only be
   at the START of any HP. So H(T+source) = (# of HP starting from v).

2. FILLING THEOREM: If H=k and H=k+2 are both achievable at n,
   then H=k+2j is achievable at n for all j in some range.
   This doesn't quite work. Instead:

3. ARC-FLIP INCREMENTS: Flipping one arc changes H by an even amount.
   If we can show the arc-flip increment can be +2 or -2 from any
   tournament, then we can reach ALL odd values in [H_min, H_max].

4. FROM OCF: H = 1 + 2*S where S = alpha_1 + 2*alpha_2 + ...
   Flipping an arc changes S by some amount Delta_S. If Delta_S can
   be +1 or -1 (reaching adjacent S values), the spectrum is dense.

Let me test approach 3: what is the distribution of arc-flip increments?
"""

from itertools import combinations
from collections import Counter, defaultdict
from math import comb

def adj_from_bits(bits, n):
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def held_karp_H(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def analyze_arc_flips():
    """For each tournament, compute H after flipping each arc.
    The arc-flip increment Delta_H = H(T') - H(T) where T' differs
    by one arc reversal.
    """
    print("=" * 70)
    print("ARC-FLIP INCREMENT ANALYSIS")
    print("=" * 70)

    for n in [5, 6]:
        m = n * (n - 1) // 2
        N = 1 << m

        print(f"\n  --- n={n} ---")

        # For each tournament, find all arc-flip neighbors
        all_deltas = Counter()
        can_reach_plus2 = 0
        can_reach_minus2 = 0
        total = 0

        for bits in range(N):
            adj = adj_from_bits(bits, n)
            H = held_karp_H(adj, n)

            deltas_this = set()
            for edge in range(m):
                # Flip edge 'edge'
                flipped_bits = bits ^ (1 << edge)
                adj_flip = adj_from_bits(flipped_bits, n)
                H_flip = held_karp_H(adj_flip, n)
                delta = H_flip - H
                deltas_this.add(delta)
                all_deltas[delta] += 1

            if 2 in deltas_this:
                can_reach_plus2 += 1
            if -2 in deltas_this:
                can_reach_minus2 += 1
            total += 1

        print(f"  Arc-flip delta distribution:")
        for d in sorted(all_deltas.keys()):
            print(f"    Delta={d:+4d}: {all_deltas[d]} occurrences")

        print(f"\n  Tournaments that can reach H+2: {can_reach_plus2}/{total} "
              f"({can_reach_plus2/total:.1%})")
        print(f"  Tournaments that can reach H-2: {can_reach_minus2}/{total} "
              f"({can_reach_minus2/total:.1%})")

        # KEY: Is every achievable H value reachable from an adjacent one?
        # If every tournament with H=k can reach some tournament with H=k+2 or H=k-2,
        # then the achievable set is an interval (connected).

        H_vals = defaultdict(int)
        for bits in range(N):
            adj = adj_from_bits(bits, n)
            H = held_karp_H(adj, n)
            H_vals[H] += 1

        achievable = sorted(H_vals.keys())
        print(f"\n  Achievable H values: {achievable}")
        print(f"  Odd gaps: {sorted(set(range(1, max(achievable)+1, 2)) - set(achievable))}")

        # For each achievable H, can we reach H-2 and H+2?
        print(f"\n  Reachability via arc flips:")
        for h in achievable:
            can_up = False
            can_down = False
            for bits in range(N):
                adj = adj_from_bits(bits, n)
                if held_karp_H(adj, n) != h:
                    continue
                for edge in range(m):
                    H_flip = held_karp_H(adj_from_bits(bits ^ (1 << edge), n), n)
                    if H_flip == h + 2:
                        can_up = True
                    if H_flip == h - 2:
                        can_down = True
                if can_up and can_down:
                    break
            print(f"    H={h}: can_reach(H+2)={can_up}, can_reach(H-2)={can_down}")


def analyze_vertex_addition_spectrum():
    """Track how vertex addition fills in the spectrum."""
    print("\n" + "=" * 70)
    print("VERTEX ADDITION SPECTRUM FILLING")
    print("=" * 70)

    # For each n, track which NEW odd values become achievable
    cumulative = set()
    for n in range(3, 8):
        m = n * (n - 1) // 2
        new_this_n = set()
        for bits in range(1 << m):
            adj = adj_from_bits(bits, n)
            H = held_karp_H(adj, n)
            if H % 2 == 1 and H not in cumulative:
                new_this_n.add(H)
                cumulative.add(H)

        max_H = max(cumulative)
        missing = sorted(set(range(1, max_H + 1, 2)) - cumulative)
        density = len(cumulative) / ((max_H + 1) // 2)

        print(f"\n  n={n}: max H = {max_H}")
        print(f"    New values: {sorted(new_this_n)}")
        print(f"    Cumulative: {len(cumulative)} values, density = {density:.1%}")
        print(f"    Still missing: {missing[:20]}{'...' if len(missing) > 20 else ''}")


if __name__ == "__main__":
    analyze_vertex_addition_spectrum()
    analyze_arc_flips()
