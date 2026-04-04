"""
disjoint_pair_formula.py — opus-2026-04-04-S1

Find the exact formula for c₂ of vertex-disjoint tile pairs.

From the n=7 data, the ratio c₂/2^(s₁+s₂-2) depends on "overlap"
(min(x₁,x₂)-max(y₁,y₂)), but overlap alone isn't sufficient — pairs
with the same (skips, overlap) can have different c₂.

Key observations:
- When ranges don't overlap (overlap ≤ 0): c₂ = 2^(s₁+s₂-2) (ratio=1)
  WAIT: actually overlap ≤ 0 pairs at n=7 have c₂=4, and 2^(2+2-2)=4, so ratio=1.
  But these are s=(2,2) with overlap=-1 or -2.

- When ranges overlap but tiles are "far": ratio = 1/4 mostly
- Intermediate: ratio = 5/16, 3/8, 1/2, etc.

Hypothesis: c₂ depends on the number of Hamiltonian paths through BOTH arcs
simultaneously. Let me compute this directly.
"""
import numpy as np
from collections import defaultdict
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles, tiling_to_tournament, count_hamiltonian_paths
from fractions import Fraction

def count_paths_using_arc(n, tiles, tiling_bits, arc):
    """Count Hamiltonian paths in tournament T(tiling_bits) that use the arc (u,v)."""
    T = tiling_to_tournament(n, tiles, tiling_bits)
    u, v = arc
    if not T[u][v]:
        return 0  # arc not present

    # Count HP using u→v: DP over (mask, endpoint) tracking if u→v used
    # Actually, simpler: count all HP containing edge u→v
    # = count paths where u immediately precedes v

    # DP: dp[mask][last] = count of paths visiting mask ending at last
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for s in range(n):
        dp[1 << s][s] = 1

    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)) or dp[mask][last] == 0:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if T[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]

    # Count paths using arc u→v: those where v immediately follows u
    # = sum over all full paths where at position k we have u, at k+1 we have v
    # In DP terms: dp[mask_without_v][u] * (paths from v to end using remaining)
    # But this is complicated. Instead, use: remove edge u→v, count HP with and without.

    # Total HP
    total = sum(dp[full])

    # HP NOT using u→v: in a tournament where u→v is removed (replaced by v→u)
    T2 = T.copy()
    T2[u][v] = 0
    T2[v][u] = 1
    total_without = count_hamiltonian_paths(T2, n)

    # HP using u→v
    using = total - total_without
    # Hmm, this is the number of HP where u→v appears as a consecutive pair.
    # Actually, removing u→v and adding v→u changes the tournament.
    # Some paths that used u→v are lost, some new paths using v→u are gained.
    # total_without counts paths in the MODIFIED tournament, not in the original
    # minus those using u→v.

    # Better approach: count directly from DP
    # Count paths where u is at position k and v at position k+1
    count_uv = 0
    for k in range(n-1):
        # Paths where position k = u, position k+1 = v
        # = (paths visiting some mask ending at u at position k) *
        #   (v is reachable from u, i.e., T[u][v]) *
        #   (paths from v through remaining vertices)
        pass

    # Actually the simplest: modify DP to track whether u→v was used
    # Use a separate DP
    # dp2[mask][last][used_uv] where used_uv ∈ {0,1}
    dp2 = np.zeros((1 << n, n, 2), dtype=np.int64)
    for s in range(n):
        dp2[1 << s][s][0] = 1

    for mask in range(1, 1 << n):
        for last in range(n):
            for used in range(2):
                if dp2[mask][last][used] == 0:
                    continue
                for nxt in range(n):
                    if mask & (1 << nxt):
                        continue
                    if T[last][nxt]:
                        new_used = used
                        if last == u and nxt == v:
                            new_used = 1
                        dp2[mask | (1 << nxt)][nxt][new_used] += dp2[mask][last][used]

    count_uv = sum(dp2[full][e][1] for e in range(n))
    return count_uv

def analyze_disjoint_formula(n=6):
    """For each disjoint pair, compute c₂ and the number of paths using both arcs."""
    H, tiles, m = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    print(f"n={n}, m={m}")
    print(f"\nDisjoint pair analysis:")
    print(f"{'Tiles':<22} {'s1,s2':<8} {'ov':>4} {'c₂':>6} {'2^(Σs-2)':>10} "
          f"{'ratio':>10} {'gap':>4}")

    results = []
    for S, c in sorted(coeffs.items()):
        if bin(S).count('1') != 2:
            continue
        bits = [i for i in range(m) if S & (1 << i)]
        t1, t2 = tiles[bits[0]], tiles[bits[1]]
        x1, y1 = t1
        x2, y2 = t2

        # vertex disjoint?
        if {x1,y1} & {x2,y2}:
            continue

        s1, s2 = x1-y1, x2-y2
        ov = min(x1,x2) - max(y1,y2)
        base = 2**(s1+s2-2)
        ratio = Fraction(c, base)

        # "Gap" = number of vertices strictly between the two tile ranges
        # that are not endpoints of either tile
        range1 = set(range(y1, x1+1))
        range2 = set(range(y2, x2+1))
        # Number of shared internal vertices
        shared_internal = len(range1 & range2) - len({x1,y1} & {x2,y2})

        gap = max(0, max(y1,y2) - min(x1,x2))  # vertices between ranges

        results.append({
            'tiles': (t1,t2), 's': (s1,s2), 'ov': ov, 'c2': c,
            'base': base, 'ratio': ratio, 'gap': gap,
            'shared_internal': shared_internal,
        })

        print(f"({x1},{y1})-({x2},{y2})  ({s1},{s2})  {ov:>4} {c:>6} {base:>10} "
              f"{str(ratio):>10} {gap:>4}")

    # Test hypothesis: c₂ = 2·N_{both} where N_{both} = expected number of
    # Hamiltonian paths using both arcs in the transitive tournament
    # Actually, c₂ = H(both flipped) - H(first flipped) - H(second flipped) + H(none flipped)
    # This IS the inclusion-exclusion for "interaction"

    # Key insight: maybe c₂ = 2 * (number of independent sets in Omega containing
    # cycles that use both flipped arcs)

    # Let me try a different formula. For completely non-overlapping tiles:
    # c₂ = 2^(s₁+s₂-2) because flipping both is like two independent flips
    # For overlapping tiles, some paths are shared, reducing c₂

    # Hypothesis: c₂ = 2^(s₁+s₂-2) - 2^(something involving overlap)
    print(f"\nDifference from 2^(Σs-2):")
    for r in sorted(results, key=lambda x: x['ov']):
        diff = r['base'] - r['c2']
        if diff == 0:
            print(f"  {r['tiles']}: c₂ = base (no correction)")
        else:
            diff_frac = Fraction(diff, r['base'])
            print(f"  {r['tiles']} s={r['s']} ov={r['ov']}: "
                  f"c₂={r['c2']}, base={r['base']}, "
                  f"deficit={diff}={diff_frac}*base")

def test_formula_at_n5_n6():
    """Look for closed-form formula across multiple n."""
    for n in [5, 6, 7]:
        H, tiles, m = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        print(f"\n{'='*60}")
        print(f"DISJOINT PAIR FORMULA CHECK n={n}")
        print(f"{'='*60}")

        for S, c in sorted(coeffs.items()):
            if bin(S).count('1') != 2:
                continue
            bits = [i for i in range(m) if S & (1 << i)]
            t1, t2 = tiles[bits[0]], tiles[bits[1]]
            x1, y1 = t1
            x2, y2 = t2
            if {x1,y1} & {x2,y2}:
                continue

            s1, s2 = x1-y1, x2-y2
            ov = min(x1,x2) - max(y1,y2)

            # Hypothesis: c₂ = 2 * (2^(s₁-1)-1) * (2^(s₂-1)-1) + 2
            # when ov >= 1, and = 2^(s₁+s₂-2) when ov <= 0
            h1 = 2 * (2**(s1-1) - 1) * (2**(s2-1) - 1) + 2 if ov >= 1 else 2**(s1+s2-2)
            match1 = "✓" if h1 == c else f"✗ (pred={h1})"

            # Hypothesis 2: c₂ depends on the number of "shared base-path vertices"
            # i.e., vertices that both tiles "span" on the base path
            shared_count = max(0, min(x1,x2) - max(y1,y2) + 1)  # = ov + 1 if ov >= 0
            h2 = 2**(s1+s2-2) - (shared_count-1) * 2**(s1+s2-3) if shared_count > 1 else 2**(s1+s2-2)
            match2 = "✓" if h2 == c else f"✗ (pred={h2})"

            if match1 != "✓" and c != 2**(s1+s2-2):
                print(f"  ({x1},{y1})-({x2},{y2}) s=({s1},{s2}) ov={ov}: "
                      f"c₂={c}, h1={match1}, h2={match2}")

if __name__ == "__main__":
    test_formula_at_n5_n6()
    analyze_disjoint_formula(n=7)
