#!/usr/bin/env python3
"""
eisenstein_forbidden.py -- Eisenstein integers and forbidden residues

Two deep threads from two_three_tower_deep.py:
1. The Eisenstein integer I(w) in Z[w] and its norm -- what does N(I(w)) encode?
2. The forbidden H values mod 216 -- why only 61% of odd residues achieved?
3. The Eisenstein factorization of the resultant I(1)*N(I(w))
4. Connection between Eisenstein primes and tournament structure
5. G(1/2) = H + 1 + 2*a2 -- can we find G at other algebraic points?
6. The full achievable set of (a1, a2) pairs at n=7

Author: kind-pasteur-2026-03-14-S63
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
from math import gcd, comb
import cmath


def random_tournament(n, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_directed_cycles(A, n):
    cycles = []
    for k in range(3, n+1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            nc = count_ham_cycles(A, verts)
            for _ in range(nc):
                cycles.append(frozenset(subset))
    return cycles


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp or dp[(mask, v)] == 0:
                continue
            cnt = dp[(mask, v)]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and dp[(full, v)] > 0:
            if A[verts[v]][verts[0]]:
                total += dp[(full, v)]
    return total


def compute_alpha_1_2(cycles):
    alpha_1 = len(cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                alpha_2 += 1
    return alpha_1, alpha_2


def main():
    n = 7
    rng = np.random.default_rng(42)

    # Collect larger sample
    num_samples = 1000
    tournaments = []
    print(f"Collecting {num_samples} tournaments at n={n}...")
    for i in range(num_samples):
        A = random_tournament(n, rng)
        cycles = get_directed_cycles(A, n)
        a1, a2 = compute_alpha_1_2(cycles)
        H = 1 + 2*a1 + 4*a2
        tournaments.append({'a1': a1, 'a2': a2, 'H': H})
        if (i+1) % 200 == 0:
            print(f"  {i+1}/{num_samples} done")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 1: The achievable (a1, a2) lattice at n=7")
    print("=" * 70)

    a1_vals = sorted(set(t['a1'] for t in tournaments))
    a2_vals = sorted(set(t['a2'] for t in tournaments))
    pairs = set((t['a1'], t['a2']) for t in tournaments)

    print(f"\n  a1 range: [{min(a1_vals)}, {max(a1_vals)}]")
    print(f"  a2 range: [{min(a2_vals)}, {max(a2_vals)}]")
    print(f"  Distinct (a1,a2) pairs: {len(pairs)}")
    print(f"  Distinct H values: {len(set(t['H'] for t in tournaments))}")

    # Check parity constraints on a1
    a1_parity = defaultdict(int)
    for t in tournaments:
        a1_parity[t['a1'] % 2] += 1
    print(f"\n  a1 parity: {dict(a1_parity)}")

    # Check if a1 mod 7 has structure (since n=7)
    a1_mod7 = defaultdict(int)
    for t in tournaments:
        a1_mod7[t['a1'] % 7] += 1
    print(f"  a1 mod 7: {dict(sorted(a1_mod7.items()))}")

    # Check if a2 has parity constraint
    a2_parity = defaultdict(int)
    for t in tournaments:
        a2_parity[t['a2'] % 2] += 1
    print(f"  a2 parity: {dict(a2_parity)}")

    # a2 mod 3
    a2_mod3 = defaultdict(int)
    for t in tournaments:
        a2_mod3[t['a2'] % 3] += 1
    print(f"  a2 mod 3: {dict(sorted(a2_mod3.items()))}")

    # a1 + a2 parity (this determines whether b0 is odd or even)
    a1_a2_parity = defaultdict(int)
    for t in tournaments:
        a1_a2_parity[(t['a1'] + t['a2']) % 2] += 1
    print(f"  (a1+a2) parity: {dict(a1_a2_parity)}")

    # The constraint: a2 <= C(a1, 2) / something... what are the bounds?
    # At n=7: max a2 seen?
    print(f"\n  a2 vs a1 bounds:")
    for a1_val in [7, 14, 21, 28, 35, 42, 49, 54]:
        a2_for_a1 = [t['a2'] for t in tournaments if t['a1'] == a1_val]
        if a2_for_a1:
            print(f"    a1={a1_val:2d}: a2 in [{min(a2_for_a1)}, {max(a2_for_a1)}], "
                  f"a2/a1={max(a2_for_a1)/a1_val:.3f}, count={len(a2_for_a1)}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 2: Forbidden H values mod 6^k")
    print("=" * 70)

    # H = 1 + 2*a1 + 4*a2
    # What are the ACHIEVABLE H values?
    all_H = sorted(set(t['H'] for t in tournaments))
    print(f"\n  Achievable H values: {len(all_H)}")
    print(f"  Range: [{min(all_H)}, {max(all_H)}]")

    # H mod 6
    print(f"\n  H mod 6: {sorted(set(h % 6 for h in all_H))} (all odd residues)")

    # H mod 12
    h_mod12 = sorted(set(h % 12 for h in all_H))
    print(f"  H mod 12: {h_mod12}")
    for r in h_mod12:
        cnt = sum(1 for h in all_H if h % 12 == r)
        print(f"    {r:3d}: {cnt} values")

    # H mod 36
    h_mod36 = sorted(set(h % 36 for h in all_H))
    missing_mod36 = sorted(set(r for r in range(1, 36, 2)) - set(h_mod36))
    print(f"\n  H mod 36: {len(h_mod36)}/18 achievable")
    print(f"  Missing: {missing_mod36}")

    # H mod 216
    h_mod216 = sorted(set(h % 216 for h in all_H))
    all_odd_216 = set(r for r in range(1, 216, 2))
    missing_mod216 = sorted(all_odd_216 - set(h_mod216))
    print(f"\n  H mod 216: {len(h_mod216)}/{len(all_odd_216)} achievable ({100*len(h_mod216)/len(all_odd_216):.1f}%)")
    print(f"  Missing residues ({len(missing_mod216)}):")
    # Group missing by mod 36
    for m36 in range(36):
        missing_in_m36 = [r for r in missing_mod216 if r % 36 == m36]
        if missing_in_m36:
            print(f"    mod 36={m36:2d}: {missing_in_m36}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 3: H gaps -- impossible H values")
    print("=" * 70)

    # Which odd values are NOT achievable as H?
    h_set = set(all_H)
    max_H = max(all_H)

    # Check small gaps
    gaps = []
    for h in range(1, max_H + 1, 2):
        if h not in h_set:
            gaps.append(h)

    print(f"\n  Impossible odd H values in [1, {max_H}]: {len(gaps)}")
    if len(gaps) <= 30:
        print(f"  Values: {gaps}")
    else:
        print(f"  First 30: {gaps[:30]}")

    # Known gaps from previous work
    print(f"\n  H=7 known impossible (THM-029). In our gaps: {7 in gaps}")
    print(f"  H=21 known strong evidence for gap. In our gaps: {21 in gaps}")

    # Gap density by range
    for lo, hi in [(1, 50), (51, 100), (101, 157)]:
        n_odd = len(range(lo, hi+1, 2))
        n_gaps = len([g for g in gaps if lo <= g <= hi])
        print(f"  [{lo:3d},{hi:3d}]: {n_gaps}/{n_odd} gaps ({100*n_gaps/n_odd:.1f}%)")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 4: Eisenstein norm N(I(w)) distribution")
    print("=" * 70)

    norms = []
    for t in tournaments:
        a1, a2 = t['a1'], t['a2']
        rp = 1 - a2
        wp = a1 - a2
        N = rp**2 - rp*wp + wp**2
        norms.append(N)

    print(f"\n  N(I(w)) range: [{min(norms)}, {max(norms)}]")
    print(f"  Mean: {np.mean(norms):.1f}")
    print(f"  Std: {np.std(norms):.1f}")

    # N mod 3 -- when is 3|N(I(w))?
    n3 = [N % 3 for N in norms]
    print(f"\n  N(I(w)) mod 3: {{0: {n3.count(0)}, 1: {n3.count(1)}, 2: {n3.count(2)}}}")
    print(f"  3|N(I(w)) iff (1-w)|I(w) in Z[w] iff I(1)=0 mod 3")

    # Verify: 3|N iff I(1)=0 mod 3
    for t, N in zip(tournaments[:20], norms[:20]):
        I1 = 1 + t['a1'] + t['a2']
        print(f"    a1={t['a1']:2d}, a2={t['a2']:2d}: I(1)={I1}, I(1)%3={I1%3}, N={N}, N%3={N%3}, "
              f"match={(I1%3==0)==(N%3==0)}")

    # N mod 7 -- since n=7 is special
    n7_dist = defaultdict(int)
    for N in norms:
        n7_dist[N % 7] += 1
    print(f"\n  N(I(w)) mod 7: {dict(sorted(n7_dist.items()))}")

    # Is N always a sum of form a^2-ab+b^2? (yes by definition, but what values appear?)
    # These are exactly the norms of Eisenstein integers
    # An integer m is a norm iff every prime p=2 mod 3 divides m to an even power
    print("\n  Is N(I(w)) always an Eisenstein norm? (by construction, yes)")
    print("  But which primes divide N?")

    # Factor a few norms
    def small_factor(n):
        factors = {}
        for p in range(2, min(n+1, 1000)):
            while n > 1 and n % p == 0:
                factors[p] = factors.get(p, 0) + 1
                n //= p
        if n > 1:
            factors[n] = 1
        return factors

    for N in sorted(set(norms))[:20]:
        f = small_factor(N)
        f_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))
        # Check if 2-mod-3 primes appear to even power
        eis_ok = all(e % 2 == 0 for p, e in f.items() if p % 3 == 2)
        print(f"    N={N:5d} = {f_str}, Eis-valid={eis_ok}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 5: I(w) in Eisenstein integers -- factorization")
    print("=" * 70)

    # I(w) = (1-a2) + (a1-a2)*w
    # When is I(w) an Eisenstein prime?
    # Eisenstein primes: p with p=2 mod 3 (inert), or pi with N(pi) = p for p=1 mod 3

    print("\n  I(w) Eisenstein factorization (via norm):")
    for t in tournaments[:15]:
        a1, a2 = t['a1'], t['a2']
        rp = 1 - a2
        wp = a1 - a2
        N = rp**2 - rp*wp + wp**2
        f = small_factor(N)
        n_prime_factors = sum(f.values())
        f_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))
        print(f"    a1={a1:2d}, a2={a2:2d}: N={N:5d} = {f_str}, "
              f"H={t['H']}, #factors={n_prime_factors}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 6: Resultant and its factorization")
    print("=" * 70)

    print("\n  Res(I(x), x^3-1) = I(1)*I(w)*I(w^2) = I(1)*N(I(w))")

    for t in tournaments[:10]:
        a1, a2 = t['a1'], t['a2']
        I1 = 1 + a1 + a2
        rp = 1 - a2
        wp = a1 - a2
        N = rp**2 - rp*wp + wp**2
        Res = I1 * N
        f = small_factor(Res)
        f_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))

        # The resultant tells us: I(x) and x^3-1 share a root mod p
        # iff p | Res
        print(f"    a1={a1:2d}, a2={a2:2d}: Res={Res:8d} = {f_str}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 7: G(z) at algebraic points")
    print("=" * 70)

    print("\n  G(z) = 1/(1-z) + a1*z/(1-z)^2 + a2*z*(1+z)/(1-z)^3")

    print("\n  Known: G(1/2) = H + 1 + 2*a2")
    print("         G(1/3) = 3/2 + 3a1/4 + 3a2/2")
    print()

    # At z = -1 (alternating series, divergent but formally):
    # G(-1) = 1/2 - a1/4 + 0 = 1/2 - a1/4 (if you use (-1)^n Abel sum)
    # Actually: z=-1 is outside radius, but the algebraic formula gives:
    # 1/(1-(-1)) = 1/2
    # a1*(-1)/(1-(-1))^2 = -a1/4
    # a2*(-1)*(1+(-1))/(1-(-1))^3 = 0
    print("  G(-1) [formal]: 1/2 - a1/4 + 0 = (2-a1)/4")
    print("  This is the alternating-sum generating function!")

    for t in tournaments[:10]:
        a1, a2 = t['a1'], t['a2']
        G_neg1 = (2 - a1) / 4
        # sum_{b>=0} I(b) * (-1)^b = I(0) - I(1) + I(2) - I(3) + ...
        # Partial sum
        ps = sum((1 + a1*b + a2*b**2) * (-1)**b for b in range(200))
        # This should NOT converge, but Cesaro/Abel sum might work
        abel = sum((1 + a1*b + a2*b**2) * (-0.99)**b for b in range(500))
        print(f"    a1={a1:2d}, a2={a2:2d}: G(-1)_formal={(2-a1)/4:.3f}, "
              f"Abel(z=-0.99)={abel:.3f}")

    # At z = 1/6 (the natural base):
    print("\n  G(1/6):")
    for t in tournaments[:5]:
        a1, a2 = t['a1'], t['a2']
        z = 1/6
        G = 1/(1-z) + a1*z/(1-z)**2 + a2*z*(1+z)/(1-z)**3
        # Closed form: 6/5 + a1*6/25 + a2*(1/6)*(7/6)/(125/216) = a2*7*216/(6*6*125)
        print(f"    a1={a1:2d}, a2={a2:2d}: G(1/6)={G:.6f}")

    # At z = (sqrt(5)-1)/2 = 1/phi (golden ratio reciprocal):
    phi = (1 + 5**0.5) / 2
    z_gold = 1 / phi
    print(f"\n  G(1/phi) where phi={phi:.6f}, 1/phi={z_gold:.6f}:")
    for t in tournaments[:5]:
        a1, a2 = t['a1'], t['a2']
        z = z_gold
        G = 1/(1-z) + a1*z/(1-z)**2 + a2*z*(1+z)/(1-z)**3
        # 1/(1-1/phi) = phi, z/(1-z)^2 = (1/phi)/(1-1/phi)^2 = phi^2/phi = phi
        # Wait: 1-1/phi = 1/phi^2... no: 1-1/phi = (phi-1)/phi = 1/phi^2? No.
        # phi = (1+sqrt5)/2, 1/phi = (sqrt5-1)/2, 1-1/phi = (3-sqrt5)/2
        print(f"    a1={a1:2d}, a2={a2:2d}: G(1/phi)={G:.4f}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 8: The (a1, a2) -> H map -- fibers and multiplicities")
    print("=" * 70)

    # H = 1 + 2*a1 + 4*a2. Given H, what (a1,a2) are possible?
    # a2 = (H - 1 - 2*a1) / 4, so a1 determines a2 (given H)
    # We need a2 >= 0, a1 >= 0, and (H-1-2*a1) divisible by 4 and >= 0
    # So a1 <= (H-1)/2, and H-1-2*a1 = 0 mod 4, i.e., a1 = (H-1)/2 mod 2

    print("\n  For each achievable H, how many (a1,a2) pairs realize it?")
    h_to_pairs = defaultdict(set)
    for t in tournaments:
        h_to_pairs[t['H']].add((t['a1'], t['a2']))

    fiber_sizes = [(h, len(ps)) for h, ps in h_to_pairs.items()]
    fiber_sizes.sort(key=lambda x: -x[1])

    print(f"\n  Top 20 H values by fiber size:")
    for h, sz in fiber_sizes[:20]:
        pairs_list = sorted(h_to_pairs[h])
        print(f"    H={h:3d}: {sz} pairs", end="")
        if sz <= 5:
            print(f" {pairs_list}")
        else:
            print(f" [{pairs_list[0]},...,{pairs_list[-1]}]")

    # H values with unique (a1,a2)
    unique_fiber = [h for h, sz in fiber_sizes if sz == 1]
    print(f"\n  H values with unique (a1,a2): {len(unique_fiber)}")
    if len(unique_fiber) <= 20:
        for h in sorted(unique_fiber):
            p = list(h_to_pairs[h])[0]
            print(f"    H={h}: (a1,a2)={p}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 9: H mod 2, 3, 4, 6, 8, 12, 24 -- the full cascade")
    print("=" * 70)

    all_H_set = set(t['H'] for t in tournaments)
    for mod in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 24, 36, 72]:
        residues = sorted(set(h % mod for h in all_H_set))
        possible = sorted(r for r in range(mod) if r % 2 == 1) if mod % 2 == 0 else list(range(mod))
        coverage = len(residues) / len(possible) if possible else 0
        print(f"  mod {mod:3d}: {len(residues):3d}/{len(possible):3d} "
              f"({100*coverage:.0f}%) residues={residues}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 10: The constraint b0 + b1 + b2 = 1 in mod arithmetic")
    print("=" * 70)

    print("\n  b0 + b1 + b2 = I(0) = 1 (always)")
    print("  Mod 2: b0+b1+b2 = 1, so an odd number of {b0,b1,b2} are odd")
    print("  Mod 3: b0+b1+b2 = 1 mod 3")

    # Distribution of (b0 mod 2, b1 mod 2, b2 mod 2)
    parity_dist = defaultdict(int)
    for t in tournaments:
        b0 = 1 - t['a1'] + t['a2']
        b1 = t['a1'] - 2*t['a2']
        b2 = t['a2']
        parity_dist[(b0 % 2, b1 % 2, b2 % 2)] += 1

    print(f"\n  (b0%2, b1%2, b2%2) distribution:")
    for key in sorted(parity_dist):
        total_odd = sum(key)
        print(f"    {key}: {parity_dist[key]:3d} (sum mod 2 = {total_odd % 2})")

    # Mod 3 distribution
    mod3_dist = defaultdict(int)
    for t in tournaments:
        b0 = 1 - t['a1'] + t['a2']
        b1 = t['a1'] - 2*t['a2']
        b2 = t['a2']
        mod3_dist[(b0 % 3, b1 % 3, b2 % 3)] += 1

    print(f"\n  (b0%3, b1%3, b2%3) distribution:")
    for key in sorted(mod3_dist):
        s = sum(key) % 3
        if mod3_dist[key] >= 5:
            print(f"    {key}: {mod3_dist[key]:3d} (sum%3={s})")

    # All sums should be 1 mod 3
    all_sum1 = all(sum(k) % 3 == 1 for k in mod3_dist.keys())
    print(f"\n  All (b0+b1+b2) = 1 mod 3: {all_sum1}")

    print("\nDone.")


if __name__ == '__main__':
    main()
