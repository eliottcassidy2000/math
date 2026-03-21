#!/usr/bin/env python3
"""layers_extended.py — Push the layer decomposition further.

Session: kind-pasteur-2026-03-20-S7

EXTENDING THM-260:
1. Compute P(n) for n=9..12 using the Burnside formula
2. Study the excess sequence E(n) = D(n)/2 - C(2k,k) — does IT have layers?
3. The "layer generating function" — can we write D as an EGF/OGF?
4. Cheeger constant of the tournament arc-flip graph
5. Second-order layers: the excess E(n) as a deficit of a deficit
"""

from math import factorial, gcd, comb, log, sqrt, pi
from fractions import Fraction
from collections import defaultdict

def partitions_odd(n):
    """Generate all partitions of n into odd parts."""
    def _gen(n, max_part):
        if n == 0:
            yield ()
            return
        for first in range(min(n, max_part), 0, -1):
            if first % 2 == 0:
                continue
            for rest in _gen(n - first, first):
                yield (first,) + rest
    yield from _gen(n, n)

def cycle_type_to_a(partition, n):
    a = [0] * (n + 1)
    for part in partition:
        a[part] += 1
    return a

def num_perms(a, n):
    result = factorial(n)
    for k in range(1, n + 1):
        if a[k] > 0:
            result //= (k ** a[k]) * factorial(a[k])
    return result

def count_free_edge_orbits(partition, n):
    """Count free edge-orbit-pairs for a permutation with given odd partition.
    Returns number of free binary choices, or -1 if impossible (self-reverse orbit).
    """
    a = cycle_type_to_a(partition, n)

    # Build concrete permutation
    perm = list(range(n))
    pos = 0
    for k in range(1, n + 1):
        for _ in range(a[k]):
            for i in range(k - 1):
                perm[pos + i] = pos + i + 1
            perm[pos + k - 1] = pos
            pos += k

    visited = [[False]*n for _ in range(n)]
    free = 0
    for i in range(n):
        for j in range(n):
            if i == j or visited[i][j]:
                continue
            orbit = []
            ci, cj = i, j
            while not visited[ci][cj]:
                visited[ci][cj] = True
                orbit.append((ci, cj))
                ci, cj = perm[ci], perm[cj]
            reverse_in = any(aa == j and bb == i for aa, bb in orbit)
            if not reverse_in:
                ci, cj = j, i
                while not visited[ci][cj]:
                    visited[ci][cj] = True
                    ci, cj = perm[ci], perm[cj]
                free += 1
            else:
                return -1
    return free


def compute_P_and_layers(max_n):
    """Compute P(n), D(n), and full layer decomposition."""
    results = {}

    for n in range(1, max_n + 1):
        total_P = 0
        total_T = 0
        layers = {}

        for partition in partitions_odd(n):
            a = cycle_type_to_a(partition, n)
            N = num_perms(a, n)
            fix = a[1]
            f = count_free_edge_orbits(partition, n)

            if f < 0:
                continue
            c_T = 2 ** f

            total_T += N * c_T
            total_P += N * c_T * fix

            is_identity = (partition == (1,) * n)
            if not is_identity:
                non_one = tuple(p for p in partition if p > 1)
                layers[non_one] = layers.get(non_one, Fraction(0)) + \
                    Fraction(N * c_T * (n - fix), factorial(n))

        T_n = total_T // factorial(n)
        P_n = total_P // factorial(n)
        D_n = n * T_n - P_n

        results[n] = {
            'T': T_n, 'P': P_n, 'D': D_n,
            'layers': dict(sorted(layers.items())),
        }

        print(f"  n={n}: T={T_n}, P={P_n}, D={D_n}")

    return results


def main():
    print("=" * 70)
    print("EXTENDED LAYER ANALYSIS — P(n) for n up to 10")
    print("=" * 70)

    results = compute_P_and_layers(10)

    # Summary table
    print(f"\n  {'='*80}")
    print(f"  COMPLETE P(n) TABLE")
    print(f"  {'='*80}")
    print(f"\n  {'n':>3} {'T(n)':>10} {'P(n)':>12} {'D(n)':>10} {'D/2':>10} {'C(2k,k)':>10} {'Excess':>10}")

    known_T = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536, 10:9733056}

    for n in sorted(results.keys()):
        r = results[n]
        k = n - 3
        cb = comb(2*k, k) if k >= 0 else 0
        D_half = r['D'] // 2 if r['D'] % 2 == 0 else r['D'] / 2
        excess = D_half - cb if k >= 0 else 0
        T_check = known_T.get(n, '?')
        t_match = ' ok' if r['T'] == T_check else f' MISMATCH({T_check})'
        print(f"  {n:3d} {r['T']:10d} {r['P']:12d} {r['D']:10d} {int(D_half):10d} "
              f"{cb:10d} {int(excess):10d}{t_match}")

    # Layer contributions table
    print(f"\n  {'='*80}")
    print(f"  LAYER CONTRIBUTIONS")
    print(f"  {'='*80}")

    all_layer_keys = set()
    for n in results:
        all_layer_keys.update(results[n]['layers'].keys())

    for key in sorted(all_layer_keys, key=lambda k: (sum(k), len(k), k)):
        print(f"\n  Layer {key} (activates at n={sum(key)}):")
        for n in sorted(results.keys()):
            if key in results[n]['layers']:
                val = results[n]['layers'][key]
                print(f"    n={n}: {float(val):14.4f}")

    # EXCESS SEQUENCE ANALYSIS
    print(f"\n  {'='*80}")
    print(f"  THE EXCESS SEQUENCE E(n) = D(n)/2 - C(2(n-3), n-3)")
    print(f"  {'='*80}")

    excess_seq = []
    for n in range(3, min(11, max(results.keys()) + 1)):
        if n in results:
            k = n - 3
            cb = comb(2*k, k)
            E = results[n]['D'] // 2 - cb
            excess_seq.append(E)
            print(f"  E({n}) = {E}")

    print(f"\n  Excess sequence: {excess_seq}")

    # Does the excess have its OWN layer structure?
    # E(n) = D(n)/2 - C(2k,k)
    # If we define E as a "second-order deficit", does it decompose further?
    print(f"\n  Excess ratios:")
    for i in range(1, len(excess_seq)):
        if excess_seq[i-1] > 0:
            print(f"  E({i+3})/E({i+2}) = {excess_seq[i]}/{excess_seq[i-1]} = "
                  f"{excess_seq[i]/excess_seq[i-1]:.4f}")

    # CHEEGER CONSTANT CONNECTION
    print(f"\n  {'='*80}")
    print(f"  CHEEGER CONSTANT OF TOURNAMENT SPACE")
    print(f"  {'='*80}")

    # The Cheeger constant h(G) of a graph G measures the worst-case
    # bottleneck for expansion: h = min_{S} |boundary(S)| / |S|
    # where the min is over sets S with |S| <= |V|/2.
    #
    # For the tournament arc-flip graph:
    # - Vertices = tournaments on n vertices (2^{C(n,2)} total)
    # - Edges = arc flips (each tournament has C(n,2) neighbors)
    # - The Cheeger constant measures how well-connected this graph is
    #
    # Cheeger's inequality: lambda_1/2 <= h <= sqrt(2 * lambda_1)
    # where lambda_1 is the spectral gap.
    #
    # For our purposes: the arc-flip graph is the HYPERCUBE {0,1}^m
    # where m = C(n,2). The Cheeger constant of the m-dimensional
    # hypercube is EXACTLY 1 (achieved by any coordinate hyperplane cut).
    #
    # But the QUOTIENT graph (tournament iso classes connected by arc flips
    # on representatives) has a DIFFERENT Cheeger constant.
    # The "bottleneck" tournaments are those that are hardest to reach
    # from the rest — potentially the Paley tournaments or transitive ones.

    print(f"""
  The tournament arc-flip graph is the m-hypercube {{0,1}}^m, m = C(n,2).
  Cheeger constant h(Q_m) = 1 (each hyperplane cut removes 2^{{m-1}} edges
  from 2^m vertices, giving boundary/volume = 1).

  The QUOTIENT by S_n symmetry changes the Cheeger constant.
  The quotient graph has T(n) vertices (iso classes) with edges
  between classes connected by arc flips.

  KEY INSIGHT: P(n) = sum of vertex orbits relates to the Cheeger
  constant via the WEIGHTED expansion:

  h_weighted = min_S sum_{{v in boundary(S)}} (1/|Aut(v)|) / sum_{{v in S}} (1/|Aut(v)|)

  The automorphism-weighted expansion measures how well the
  "perspective space" is connected. P(n) is the total weight
  of all vertices in this space.

  For the ROOTED tournament graph (P(n) vertices):
  each edge in the quotient graph lifts to O(n) edges in the rooted graph,
  so the Cheeger constant of the rooted graph is approximately
  h_rooted ~ h_quotient * n.
""")

    # SELF-SIMILARITY TEST
    print(f"  {'='*80}")
    print(f"  SELF-SIMILARITY: Does P(n) / n! approach a limit?")
    print(f"  {'='*80}")

    for n in sorted(results.keys()):
        r = results[n]
        ratio = r['P'] / factorial(n)
        ratio_T = r['T'] / factorial(n) if r['T'] > 0 else 0
        log_P = log(r['P']) if r['P'] > 0 else 0
        log_nfac = log(factorial(n))
        print(f"  n={n}: P/n! = {ratio:.8f}, T/n! = {ratio_T:.10f}, "
              f"log(P)/log(n!) = {log_P/log_nfac:.6f}" if log_nfac > 0 else "")

    # P(n)/n! should approach the "probability that a random labeled
    # tournament has a distinguished vertex up to automorphism"
    # For large n, most tournaments have trivial Aut, so P(n) ~ n*T(n),
    # and P(n)/n! ~ n*T(n)/n! = T(n)/(n-1)!

    print(f"\n  P(n) / (n * T(n)):")
    for n in sorted(results.keys()):
        r = results[n]
        if r['T'] > 0:
            ratio = r['P'] / (n * r['T'])
            print(f"  n={n}: P/(nT) = {ratio:.8f} = 1 - D/(nT) = 1 - {r['D']/(n*r['T']):.8f}")

    # THE BIG PICTURE
    print(f"\n  {'='*80}")
    print(f"  THE BIG PICTURE: P(n) AS A MULTI-SCALE EXPANSION")
    print(f"  {'='*80}")
    print(f"""
  P(n) = n*T(n) * [1 - delta(n)]

  where delta(n) = D(n) / (n*T(n)) is the "automorphism fraction."

  delta(n) measures what fraction of (vertex, class) pairs are
  redundant due to automorphisms. As n grows:

  delta(n) -> 0  (most tournaments have trivial Aut)

  But the RATE of approach has the layer structure:
  delta(n) ~ D_3(n)/(n*T(n)) + D_5(n)/(n*T(n)) + ...

  Each layer decays at a different rate, creating the multi-scale
  structure we observe.

  The dominant layer D_3 ~ 2^{{C(n-2,2)}} / (n-3)! grows as
  2^{{n^2/2}} / n!, while n*T(n) grows as n * 2^{{C(n,2)}} / n!
  (Burnside's lemma, dominated by identity permutation).

  So delta(n) ~ 2^{{C(n-2,2)}} / (n * 2^{{C(n,2)}} / n!)
  = n! * 2^{{C(n-2,2) - C(n,2)}} / n
  = (n-1)! * 2^{{-C(n,2) + C(n-2,2)}}
  = (n-1)! * 2^{{-(2n-3)}}

  Since C(n,2) - C(n-2,2) = (n-1) + (n-2) = 2n-3.

  So delta(n) ~ (n-1)! / 2^{{2n-3}} -> 0 super-exponentially.

  This SUPER-EXPONENTIAL decay of the automorphism correction
  is why the near-formula works so well: each layer decays
  faster than the previous, and the total correction D(n) is
  dominated by finitely many layers at any given n.
""")


if __name__ == "__main__":
    main()
