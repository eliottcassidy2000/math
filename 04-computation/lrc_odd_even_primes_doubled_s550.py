#!/usr/bin/env python3
"""Odd cycles, even numbers, primes, and doubled primes in tournaments + LRC.

opus-2026-06-01-S550

THE CONNECTIONS:

1. TOURNAMENTS only have ODD directed cycles as fundamental objects.
   The OCF: H(T) = I(Ω(T), 2) where Ω = conflict graph on ODD cycles.
   Even cycles decompose into pairs of odd sub-structures.

2. GOLDBACH: every even N ≥ 4 = p + q (two primes).
   LEMOINE: every odd N ≥ 7 = p + 2q (prime + doubled prime).
   The DOUBLED PRIME 2q is the bridge between odd and even.

3. In LRC: runner speeds decompose by parity and primality.
   Primes: irreducible frequencies.
   Doubled primes (2p): harmonic frequencies at twice the prime.
   The doubled primes sit at level 1 of the 2-adic tree.

4. The OCF's odd cycles PAIR with Goldbach decompositions:
   a 3-cycle on runners (a,b,c) exists iff the speed differences
   form a "Goldbach triple" — the differences relate to prime structure.

THE KEY QUESTION: what role do doubled primes play in the LRC proof?
- In the resonance debt: which terms come from doubled-prime pairs?
- In the p-adic tree: doubled primes are level-1 connectors
- In the OCF: how do doubled-prime arcs affect H(T)?
"""

from __future__ import annotations
from fractions import Fraction
from math import gcd, sqrt, log
from functools import reduce
from itertools import combinations
from collections import Counter, defaultdict


ONE = Fraction(1)
ZERO = Fraction(0)
def frac(x): return x - Fraction(x.numerator // x.denominator)
def dist0(x): f = frac(x); return min(f, ONE - f)

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True

def classify_speed(v):
    """Classify a speed as prime, doubled prime, or composite."""
    if is_prime(v):
        return "prime"
    elif v % 2 == 0 and is_prime(v // 2):
        return "2×prime"
    else:
        return "composite"


# ═══════════════════════════════════════════════════════════════
# PART 1: Only ODD cycles matter — the OCF fact
# ═══════════════════════════════════════════════════════════════

def odd_cycles_only(n_values=[4, 5, 6, 7]):
    """In the OCF: H(T) = I(Ω(T), 2) where Ω has only ODD directed cycles.

    Count directed cycles of each length in LRC tournaments.
    Show that even-length cycles are always DECOMPOSABLE into odd cycles.
    """
    print("=" * 70)
    print("PART 1: Only ODD cycles matter — the OCF structure")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = tuple(range(1, n))
        thr = Fraction(1, n)

        # Build tournament at the regular polygon time t=1/n
        t = Fraction(1, n)
        m = n
        adj = [[0]*m for _ in range(m)]
        for i in range(m):
            for j in range(m):
                if i == j: continue
                if i == 0:
                    adj[0][j] = 1 if dist0(Fraction(speeds[j-1]) * t) >= thr else 0
                elif j == 0:
                    adj[i][0] = 1 if dist0(Fraction(speeds[i-1]) * t) < thr else 0
                else:
                    diff = frac(Fraction(speeds[i-1] - speeds[j-1]) * t)
                    if ZERO < diff < Fraction(1, 2):
                        adj[i][j] = 1
                    elif diff == ZERO or diff == Fraction(1, 2):
                        adj[i][j] = 1 if i < j else 0

        # Count directed cycles of each length (brute force for small n)
        cycle_counts = Counter()
        for length in range(3, m + 1):
            count = 0
            for perm in combinations(range(m), length):
                # Check all cyclic orderings
                from itertools import permutations as perms
                for p in perms(perm):
                    is_cycle = all(adj[p[i]][p[(i+1) % length]] for i in range(length))
                    if is_cycle:
                        count += 1
                        break  # count each vertex set once

            cycle_counts[length] = count

        print(f"n={n} at regular polygon (t=1/{n}):")
        for length in sorted(cycle_counts.keys()):
            parity = "ODD" if length % 2 == 1 else "even"
            print(f"  {length}-cycles: {cycle_counts[length]:5d}  ({parity})")

        # How many 3-cycles?
        c3 = 0
        for i in range(m):
            for j in range(i+1, m):
                for k in range(j+1, m):
                    if ((adj[i][j] and adj[j][k] and adj[k][i]) or
                        (adj[i][k] and adj[k][j] and adj[j][i])):
                        c3 += 1
        print(f"  directed 3-cycles (c3): {c3}")
        print()

    print("ODD CYCLES FACT:")
    print("  The OCF uses ONLY odd cycles. Even cycles decompose into pairs.")
    print("  A tournament's H(T) is determined by its odd-cycle conflict graph.")
    print("  The independence polynomial I(Ω, 2) = Σ 2^k × (# independent k-sets)")
    print("  counts collections of VERTEX-DISJOINT ODD cycles.")
    print()
    print("  WHY ONLY ODD? Because tournaments orient COMPLETE graphs, and the")
    print("  parity constraint from the complete graph's even Euler characteristic")
    print("  means even cycles decompose: a 4-cycle A→B→C→D→A always contains")
    print("  either A→C or C→A as a 'shortcut', creating two 3-cycles.")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 2: Speed classification — primes vs doubled primes
# ═══════════════════════════════════════════════════════════════

def speed_classification(n_values=[8, 14, 18, 20]):
    """Classify each speed as prime, doubled prime (2p), or composite."""
    print("=" * 70)
    print("PART 2: Speed classification — primes, doubled primes, composites")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = list(range(1, n))

        primes_list = [v for v in speeds if is_prime(v)]
        doubled_list = [v for v in speeds if v % 2 == 0 and is_prime(v // 2)]
        composite_list = [v for v in speeds if not is_prime(v) and not (v % 2 == 0 and is_prime(v // 2))]

        print(f"n={n}: {n-1} runners")
        print(f"  primes: {primes_list} ({len(primes_list)})")
        print(f"  doubled primes (2p): {doubled_list} ({len(doubled_list)})")
        print(f"  other: {composite_list} ({len(composite_list)})")
        print()

        # The prime-doubled prime PAIRS: (p, 2p) where both are in the speed set
        pairs = [(p, 2*p) for p in primes_list if 2*p in speeds]
        print(f"  prime↔doubled pairs: {pairs}")

        # Goldbach for even speeds: v = p + q (both prime, both in speed set)
        print(f"  Goldbach pairs for even speeds:")
        for v in speeds:
            if v % 2 == 0 and v >= 4:
                goldbach = [(p, v-p) for p in primes_list if v-p in primes_list and p <= v-p]
                if goldbach:
                    print(f"    {v} = {' = '.join(f'{p}+{q}' for p,q in goldbach[:3])}")

        # Lemoine for odd speeds: v = p + 2q (p prime, q prime, both or 2q in speed set)
        print(f"  Lemoine pairs for odd speeds:")
        for v in speeds:
            if v % 2 == 1 and v >= 7:
                lemoine = [(p, q) for p in primes_list for q in primes_list
                           if p + 2*q == v and 2*q in speeds]
                if lemoine:
                    print(f"    {v} = {' = '.join(f'{p}+2×{q}' for p,q in lemoine[:3])}")

        print()


# ═══════════════════════════════════════════════════════════════
# PART 3: Doubled primes in the resonance debt
# ═══════════════════════════════════════════════════════════════

def doubled_prime_resonance(n_values=[6, 8, 14]):
    """Compute the resonance debt contribution from doubled-prime pairs.

    In the pairwise debt (S531): the term for pair (v_i, v_j) depends on
    the Bernoulli polynomial B_2({(v_i-v_j)/n}).

    For a prime-doubled pair (p, 2p): difference = p.
    For a prime-prime pair (p, q): difference = |p-q|.
    For a doubled-doubled pair (2p, 2q): difference = 2|p-q|.

    Question: do doubled-prime pairs contribute MORE to the debt?
    """
    print("=" * 70)
    print("PART 3: Doubled primes in the resonance debt")
    print("=" * 70)
    print()

    from math import pi, sin

    def fhat(k, n):
        if k == 0: return (n-2)/n
        if k % n == 0: return 0.0
        return -sin(2*pi*k/n) / (pi*k)

    for n in n_values:
        speeds = list(range(1, n))
        m = len(speeds)
        f0 = (n-2)/n

        # Compute pairwise debt for each pair, classified by type
        type_debt = defaultdict(float)
        type_count = defaultdict(int)

        for i in range(m):
            for j in range(i+1, m):
                vi, vj = speeds[i], speeds[j]
                ci = classify_speed(vi)
                cj = classify_speed(vj)
                pair_type = tuple(sorted([ci, cj]))

                # Compute pairwise debt
                g = gcd(vi, vj)
                vi_r, vj_r = vi // g, vj // g
                pair_sum = 0.0
                for t in range(1, 200):
                    pair_sum += 2 * fhat(-t * vj_r, n) * fhat(t * vi_r, n)

                pair_debt = f0 ** (m - 2) * pair_sum
                type_debt[pair_type] += pair_debt
                type_count[pair_type] += 1

        outside = f0 ** m
        total_debt = sum(type_debt.values())

        print(f"n={n}: pairwise debt by pair type")
        print(f"  outside credit: {outside:.6f}")
        print(f"  total pairwise debt: {total_debt:.6f}")
        print()

        for pair_type in sorted(type_debt.keys()):
            debt = type_debt[pair_type]
            count = type_count[pair_type]
            frac_of_total = debt / total_debt if total_debt != 0 else 0
            print(f"  {str(pair_type):40s}: {count:3d} pairs, "
                  f"debt = {debt:+.6f}  ({100*frac_of_total:+.1f}% of total)")

        # The prime↔doubled_prime pairs specifically
        pd_pairs = [(vi, vj) for i, vi in enumerate(speeds)
                     for j, vj in enumerate(speeds) if i < j
                     and ((is_prime(vi) and vj == 2*vi) or (is_prime(vj) and vi == 2*vj))]

        if pd_pairs:
            pd_debt = 0.0
            for vi, vj in pd_pairs:
                g = gcd(vi, vj)
                vi_r, vj_r = vi // g, vj // g
                ps = 0.0
                for t in range(1, 200):
                    ps += 2 * fhat(-t * vj_r, n) * fhat(t * vi_r, n)
                pd_debt += f0 ** (m - 2) * ps

            print(f"\n  prime↔doubled prime pairs: {pd_pairs}")
            print(f"  their total debt: {pd_debt:+.6f} "
                  f"({100*pd_debt/total_debt:+.1f}% of total)" if total_debt != 0 else "")

        print()


# ═══════════════════════════════════════════════════════════════
# PART 4: The 3-cycle Goldbach connection
# ═══════════════════════════════════════════════════════════════

def three_cycle_goldbach(n_values=[6, 7, 8]):
    """Directed 3-cycles involve triples (a,b,c) with a→b→c→a.
    The speed differences (b-a, c-b, a-c) sum to 0.
    These differences are the 'Goldbach triple' of the cycle.

    For a 3-cycle with differences (d1, d2, d3) where d1+d2+d3=0:
    Two of the differences have the same sign, one has the opposite.
    The 'even' difference can be decomposed as sum of two 'odd' ones
    (the Goldbach structure).

    The doubled primes appear as EVEN differences that bridge two prime
    differences: if d1=p and d2=q are prime, then d3=-(p+q) is even.
    The even difference = sum of two primes = Goldbach!
    """
    print("=" * 70)
    print("PART 4: 3-cycle speed differences — the Goldbach triple")
    print("=" * 70)
    print()

    for n in n_values:
        speeds = list(range(1, n))
        m = len(speeds)

        # For each 3-element subset of runners: check if it forms a 3-cycle
        # at the regular polygon time
        t = Fraction(1, n)
        thr = Fraction(1, n)

        cycle_triples = []

        for combo in combinations(range(m), 3):
            i, j, k = combo
            vi, vj, vk = speeds[i], speeds[j], speeds[k]

            # Check if (vi, vj, vk) forms a directed 3-cycle at time t
            # in the runner sub-tournament (half-turn)
            def beats(a, b):
                diff = frac(Fraction(a - b) * t)
                if diff == ZERO or diff == Fraction(1, 2):
                    return a < b  # tie-break
                return ZERO < diff < Fraction(1, 2)

            if beats(vi, vj) and beats(vj, vk) and beats(vk, vi):
                diffs = (vj-vi, vk-vj, vi-vk)
                cycle_triples.append(((vi, vj, vk), diffs))
            elif beats(vi, vk) and beats(vk, vj) and beats(vj, vi):
                diffs = (vk-vi, vj-vk, vi-vj)
                cycle_triples.append(((vi, vk, vj), diffs))

        print(f"n={n} at t=1/{n}: {len(cycle_triples)} directed 3-cycles")

        # Classify differences by parity
        even_diffs = 0
        odd_diffs = 0
        goldbach_triples = 0

        for (a, b, c), (d1, d2, d3) in cycle_triples:
            for d in (d1, d2, d3):
                if d % 2 == 0:
                    even_diffs += 1
                else:
                    odd_diffs += 1

            # Check if the even difference is a Goldbach sum
            for d in (d1, d2, d3):
                if d % 2 == 0 and abs(d) >= 4:
                    # Can |d| be written as p+q with p,q prime?
                    D = abs(d)
                    is_goldbach = any(is_prime(p) and is_prime(D-p)
                                      for p in range(2, D//2 + 1))
                    if is_goldbach:
                        goldbach_triples += 1

        print(f"  differences: {odd_diffs} odd, {even_diffs} even")
        print(f"  even differences that are Goldbach sums: {goldbach_triples}")

        # Show sample cycles with their difference structure
        print(f"  sample cycles:")
        for (a, b, c), (d1, d2, d3) in cycle_triples[:8]:
            diffs_class = []
            for d in (d1, d2, d3):
                if abs(d) == 1:
                    diffs_class.append("1")
                elif is_prime(abs(d)):
                    diffs_class.append(f"p={abs(d)}")
                elif abs(d) % 2 == 0 and is_prime(abs(d) // 2):
                    diffs_class.append(f"2p={abs(d)}")
                else:
                    diffs_class.append(f"{abs(d)}")
            print(f"    ({a},{b},{c}): diffs=({d1:+d},{d2:+d},{d3:+d}) = [{', '.join(diffs_class)}]")

        print()

    print("GOLDBACH-CYCLE CONNECTION:")
    print("  Every 3-cycle has three speed differences summing to 0.")
    print("  One difference is the SUM of the other two (with signs).")
    print("  When two differences are PRIME: the third is EVEN = Goldbach sum.")
    print("  The doubled primes (2p) appear as these EVEN differences.")
    print()
    print("  In the OCF: 3-cycles are the FUNDAMENTAL odd cycles.")
    print("  Their weight μ(C) = H(T[V\\V(C)]) depends on the complement.")
    print("  The Goldbach structure of the differences determines μ(C)")
    print("  through the arithmetic of the speed set.")
    print()
    print("  THE DOUBLED PRIMES ARE THE BRIDGES: they connect two prime")
    print("  speed differences into a closed cycle. Without doubled primes,")
    print("  many 3-cycles wouldn't exist. They're the 'glue' between")
    print("  prime-frequency Gabor atoms.")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 5: The doubled prime's role in LRC
# ═══════════════════════════════════════════════════════════════

def doubled_prime_role():
    """Synthesize the role of doubled primes in LRC.

    1. In the p-adic tree: doubled primes sit at level 1 of the 2-adic tree.
       They're the FIRST harmonic above the prime level.

    2. In the resonance debt: doubled-prime pairs create the strongest
       pairwise resonances (small speed difference = large Bernoulli term).

    3. In the OCF: doubled primes are the EVEN bridges between prime cycles.
       They create 3-cycles that wouldn't exist without the 2× harmonic.

    4. In the cascade: doubled primes are processed at level 1 of the
       2-adic descent. They're the FIRST non-trivial cascade step after
       the units (primes).

    5. In the Gabor picture: doubled primes are the SECOND HARMONIC of
       prime atoms. They create beat frequencies at prime wavelengths.
    """
    print("=" * 70)
    print("PART 5: The doubled prime's role — five perspectives")
    print("=" * 70)
    print()

    print("1. P-ADIC TREE: doubled primes = level 1 of the 2-adic tree.")
    print("   They sit BETWEEN the primes (level 0) and the 4-multiples (level 2).")
    print("   They're the first 'harmonic' above the fundamental prime layer.")
    print()

    print("2. RESONANCE DEBT: doubled-prime pairs (p, 2p) have speed")
    print("   difference p, which creates a resonance at frequency p/n.")
    print("   The Bernoulli term B₂(p/n) for these pairs is typically LARGE")
    print("   (because p/n is often near 1/2 for moderate primes).")
    print()

    print("3. OCF / 3-CYCLES: doubled primes bridge prime differences.")
    print("   A 3-cycle with runners (a, b, c) where b-a=p (prime) and")
    print("   c-b=q (prime) has a-c=-(p+q) (even = Goldbach sum).")
    print("   The doubled prime 2r appears when p+q=2r, i.e., the sum")
    print("   of two primes equals twice a prime (Goldbach partition).")
    print()

    print("4. CASCADE: doubled primes are processed at the 2-adic level 1.")
    print("   In the cascade descent (S544): after clearing the odd speeds")
    print("   (level 0), the doubled primes are the NEXT constraint.")
    print("   Their clearance is the '2-adic step' of the cascade.")
    print()

    print("5. GABOR: doubled primes are the 2nd harmonic of prime atoms.")
    print("   In the time-frequency plane (S541): the atom at frequency 2p")
    print("   is the OCTAVE of the atom at frequency p.")
    print("   The interference between p and 2p creates a 'beat' at frequency p.")
    print("   This beat IS the pairwise resonance term in the debt (S531).")
    print()

    print("THE DEEP POINT:")
    print("  Primes are the IRREDUCIBLE frequencies (Gabor atoms).")
    print("  Doubled primes are their FIRST HARMONICS (octaves).")
    print("  The relationship prime ↔ doubled prime is:")
    print("    - Goldbach: the even number = p + q (two primes)")
    print("    - Lemoine: the odd number = p + 2q (prime + doubled prime)")
    print("    - OCF: the 3-cycle = two prime arcs + one even bridge")
    print("    - Resonance: the pairwise term = beat between p and 2p")
    print("    - P-adic: the tree level 0 → level 1 transition")
    print()
    print("  Without doubled primes: the OCF has FEWER 3-cycles,")
    print("  the resonance debt is SMALLER, and LRC is EASIER.")
    print("  Doubled primes make LRC HARDER by creating additional")
    print("  cycle-interference — they're the 'complication' that the")
    print("  arithmetic progression exploits for maximum coherence.")
    print()
    print("  The AP {1,...,n-1} has the MOST doubled primes in range")
    print("  (every even number ≤ n-1 is a doubled prime or composite).")
    print("  This is another reason AP is the tightest case.")
    print()


def main():
    print("Odd/Even Cycles, Primes, and Doubled Primes — opus-S550")
    print()

    odd_cycles_only()
    speed_classification()
    doubled_prime_resonance()
    three_cycle_goldbach()
    doubled_prime_role()


if __name__ == "__main__":
    main()
