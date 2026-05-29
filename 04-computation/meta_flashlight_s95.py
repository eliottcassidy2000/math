#!/usr/bin/env python3
"""
meta_flashlight_s95.py

Small exhaustive probes for cross-thread hypotheses:

1. Bare Cartan decomposition of a 0/1 tournament adjacency is inert.
   The symmetric part is always (J-I)/2, and the skew norm is fixed.
   Any useful "dark mode" must therefore live in a weighted or lifted object.

2. The OCF lift T -> Omega(T) is a concrete variable symmetric sector.
   Omega is undirected, varies strongly with T, and exactly recovers H by
   I(Omega, 2) for n <= 6.

3. Look for small-n residue/fiber anomalies: what H values and Omega features
   slip through after conditioning on t3?
"""

from collections import defaultdict, Counter
from itertools import combinations, permutations
from math import comb


def bits_to_adj(bits, n):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def hamiltonian_paths(A):
    n = len(A)
    N = 1 << n
    dp = [[0] * n for _ in range(N)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(N):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if A[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[N - 1])


def directed_cycle_count(A, k):
    n = len(A)
    count = 0
    for cyc in permutations(range(n), k):
        if min(cyc) != cyc[0]:
            continue
        ok = True
        for i in range(k):
            if not A[cyc[i]][cyc[(i + 1) % k]]:
                ok = False
                break
        if ok:
            count += 1
    return count


def canonical_cycle(cyc):
    cyc = tuple(cyc)
    rots = [cyc[i:] + cyc[:i] for i in range(len(cyc))]
    return min(rots)


def odd_cycles(A):
    n = len(A)
    cycles = set()
    for k in range(3, n + 1, 2):
        for subset in combinations(range(n), k):
            for cyc in permutations(subset):
                if min(cyc) != cyc[0]:
                    continue
                if all(A[cyc[i]][cyc[(i + 1) % k]] for i in range(k)):
                    cycles.add(canonical_cycle(cyc))
    return sorted(cycles)


def omega_stats(cycles):
    m = len(cycles)
    masks = []
    for cyc in cycles:
        mask = 0
        for v in cyc:
            mask |= 1 << v
        masks.append(mask)

    conflict = [0] * m
    edges = 0
    indep_pairs = 0
    for i in range(m):
        for j in range(i + 1, m):
            if masks[i] & masks[j]:
                conflict[i] |= 1 << j
                conflict[j] |= 1 << i
                edges += 1
            else:
                indep_pairs += 1

    indep_counts = Counter()

    def backtrack(start, forbidden, size):
        indep_counts[size] += 1
        for v in range(start, m):
            if forbidden & (1 << v):
                continue
            backtrack(v + 1, forbidden | (1 << v) | conflict[v], size + 1)

    backtrack(0, 0, 0)
    I2 = sum((2 ** k) * c for k, c in indep_counts.items())
    max_deg = max((conflict[i].bit_count() for i in range(m)), default=0)
    min_deg = min((conflict[i].bit_count() for i in range(m)), default=0)
    density = edges / comb(m, 2) if m >= 2 else 0.0
    return {
        "omega_v": m,
        "omega_e": edges,
        "omega_density": density,
        "omega_min_deg": min_deg,
        "omega_max_deg": max_deg,
        "alpha2": indep_pairs,
        "I2": I2,
        "indep_counts": dict(sorted(indep_counts.items())),
    }


def odd_part(x):
    while x % 2 == 0 and x:
        x //= 2
    return x


def main():
    print("=" * 72)
    print("META FLASHLIGHT S95: WHERE THE CURRENT ANALOGIES HAVE BLIND SPOTS")
    print("=" * 72)

    print("\n1. Bare Cartan decomposition of 0/1 adjacency is inert")
    print("   For every tournament A, S=(A+A^T)/2=(J-I)/2.")
    print("   Therefore ||S||_F^2 = ||K||_F^2 = C(n,2)/2 for all T.")
    for n in range(3, 8):
        invariant = comb(n, 2) / 2
        print(f"   n={n}: symmetric_norm_sq=skew_norm_sq={invariant:.1f}")

    print("\n2. OCF lift T -> Omega(T) is a variable symmetric dark sector")
    for n in range(3, 7):
        N = 1 << comb(n, 2)
        mismatches = 0
        by_h = Counter()
        by_t3 = defaultdict(list)
        omega_v_by_t3 = defaultdict(set)
        alpha2_by_t3 = defaultdict(set)
        residue_by_h = defaultdict(set)
        D = odd_part(comb(n, 2))

        for bits in range(N):
            A = bits_to_adj(bits, n)
            H = hamiltonian_paths(A)
            t3 = directed_cycle_count(A, 3)
            cycles = odd_cycles(A)
            stats = omega_stats(cycles)
            if stats["I2"] != H:
                mismatches += 1
            by_h[H] += 1
            by_t3[t3].append(H)
            omega_v_by_t3[t3].add(stats["omega_v"])
            alpha2_by_t3[t3].add(stats["alpha2"])
            residue_by_h[H].add(H % D if D > 1 else 0)

        missing_odds = [h for h in range(1, max(by_h) + 1, 2) if h not in by_h]
        ambiguous_t3 = {
            t3: (min(vals), max(vals), len(set(vals)))
            for t3, vals in by_t3.items()
            if len(set(vals)) > 1
        }

        print(f"\n   n={n}: tournaments={N}, OCF mismatches={mismatches}")
        print(f"     H spectrum size={len(by_h)}, range=[{min(by_h)}, {max(by_h)}]")
        print(f"     odd gaps through max={missing_odds[:20]}")
        print(f"     conductor D=odd_part(C(n,2))={D}; H residues seen={sorted({h % D if D > 1 else 0 for h in by_h})}")
        print(f"     t3 fibers with multiple H values={len(ambiguous_t3)}/{len(by_t3)}")
        for t3 in sorted(ambiguous_t3)[:6]:
            lo, hi, distinct = ambiguous_t3[t3]
            print(
                f"       t3={t3:2d}: H range [{lo:3d},{hi:3d}], "
                f"distinct={distinct:2d}, |Omega|={sorted(omega_v_by_t3[t3])}, "
                f"alpha2={sorted(alpha2_by_t3[t3])[:8]}"
            )

    print("\n3. New meta-hypothesis suggested by the data")
    print("   The useful symmetric sector is not the raw Cartan S of A; it is")
    print("   the functorial lift Omega(T), transfer M(T), or weighted attention.")
    print("   Conditioning on t3 removes the first centered skew-chaos term,")
    print("   but H variation remains exactly where Omega gains independent sets.")
    print("   Search target: residual H within fixed t3 fibers, ordered by alpha2.")


if __name__ == "__main__":
    main()
