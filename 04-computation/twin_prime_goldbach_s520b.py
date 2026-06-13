#!/usr/bin/env python3
"""Twin-prime Goldbach: which even numbers can't be written as p+q
where both p,q are members of twin prime pairs?

opus-2026-06-01-S520b

Questions:
1. Find all even N that CANNOT be expressed as sum of two twin primes.
2. Are there exactly 35 such values?
3. What structure do they reveal?
4. Connection to complement necklaces and tournament theory?

A "twin prime" is any prime p such that p-2 or p+2 is also prime.
Twin prime pairs: (3,5), (5,7), (11,13), (17,19), (29,31), ...
Twin prime set: {3, 5, 7, 11, 13, 17, 19, 29, 31, 41, 43, ...}
"""

from math import isqrt, gcd
from collections import Counter, defaultdict


def sieve_primes(limit):
    """Sieve of Eratosthenes."""
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            for j in range(i * i, limit + 1, i):
                is_prime[j] = False
    return is_prime


def twin_primes(limit):
    """Return the set of primes that belong to a twin pair (p, p+2)."""
    is_prime = sieve_primes(limit + 2)
    twins = set()
    for p in range(2, limit + 1):
        if is_prime[p]:
            if (p >= 2 and p + 2 <= limit + 2 and is_prime[p + 2]):
                twins.add(p)
                twins.add(p + 2)
            if (p - 2 >= 2 and is_prime[p - 2]):
                twins.add(p)
                twins.add(p - 2)
    return sorted(twins)


def find_exceptions(search_limit=10000, prime_limit=100000):
    """Find even numbers that can't be written as sum of two twin primes."""
    tp = set(twin_primes(prime_limit))
    tp_list = sorted(tp)

    print(f"Twin primes up to {prime_limit}: {len(tp_list)}")
    print(f"First 30: {tp_list[:30]}")
    print()

    exceptions = []
    representable = []

    for N in range(2, search_limit + 1, 2):
        found = False
        for p in tp_list:
            if p > N:
                break
            q = N - p
            if q in tp:
                found = True
                break
        if not found:
            exceptions.append(N)

    return exceptions, tp_list


def analyze_exceptions(exceptions, tp_list):
    """Deep structural analysis of the exceptions."""
    print("=" * 70)
    print("PART A: The exceptions")
    print("=" * 70)
    print()
    print(f"Total exceptions found: {len(exceptions)}")
    print(f"Exceptions: {exceptions}")
    print()

    # Factorizations
    print("Factorizations:")
    for N in exceptions:
        factors = []
        n = N
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
            while n % p == 0:
                factors.append(p)
                n //= p
            if n == 1:
                break
        if n > 1:
            factors.append(n)
        print(f"  {N:4d} = {' × '.join(str(f) for f in factors)}"
              f"  mod 6={N%6}  mod 12={N%12}  mod 30={N%30}")
    print()

    # Mod distributions
    print("Mod distributions of exceptions:")
    for m in [3, 4, 6, 8, 10, 12, 30]:
        dist = Counter(N % m for N in exceptions)
        print(f"  mod {m:2d}: {dict(sorted(dist.items()))}")
    print()

    # Gap analysis
    print("Gaps between consecutive exceptions:")
    gaps = [exceptions[i+1] - exceptions[i] for i in range(len(exceptions)-1)]
    print(f"  Gaps: {gaps}")
    print(f"  Gap set: {sorted(set(gaps))}")
    print()

    # Which twin primes are "near" the exceptions?
    tp_set = set(tp_list)
    print("Nearest twin primes to each exception:")
    for N in exceptions[:20]:
        below = [p for p in tp_list if p < N]
        above = [p for p in tp_list if p > N]
        near_below = below[-1] if below else None
        near_above = above[0] if above else None
        print(f"  {N:4d}: below={near_below}, above={near_above}, "
              f"gaps=({N - near_below if near_below else '-'}, "
              f"{near_above - N if near_above else '-'})")
    print()


def complement_necklace_analysis(exceptions):
    """Look for complement/necklace structure in the exceptions.

    "Complement" here means: if N is an exception, is max_exception - N also?
    "Necklace" means: do the exceptions form orbits under some group action?
    """
    print("=" * 70)
    print("PART B: Complement and necklace structure")
    print("=" * 70)
    print()

    exc_set = set(exceptions)
    max_exc = max(exceptions) if exceptions else 0

    # Complement pairing: N <-> max_exc - N (if max_exc even)
    # Or N <-> max_exc + 2 - N
    print(f"Max exception: {max_exc}")
    print()

    # Check complement: N <-> (max + 2) - N
    target = max_exc + 2  # Try max+2 as the complement target
    comp_pairs = []
    unpaired = []
    for N in exceptions:
        comp = target - N
        if comp in exc_set and comp != N:
            if (N, comp) not in comp_pairs and (comp, N) not in comp_pairs:
                comp_pairs.append((N, comp))
        elif comp == N:
            comp_pairs.append((N,))
        else:
            unpaired.append(N)

    print(f"Complement target = {target}:")
    print(f"  Paired: {comp_pairs[:15]}...")
    print(f"  Unpaired: {unpaired}")
    print()

    # Try other complement targets
    for target in [max_exc, max_exc + 2, max_exc + 4, 2 * max_exc]:
        paired_count = sum(1 for N in exceptions if target - N in exc_set)
        print(f"  Target {target}: {paired_count}/{len(exceptions)} paired")
    print()

    # Necklace: orbits under N -> N mod various moduli
    # Are the exceptions a union of arithmetic progressions?
    print("Arithmetic progression analysis:")
    for m in [6, 10, 12, 30]:
        residues = set(N % m for N in exceptions)
        # Check if exceptions are exactly the even numbers with these residues
        # up to some bound
        predicted = sorted(N for N in range(2, max_exc + 1, 2) if N % m in residues)
        actual = sorted(exceptions)
        match = sum(1 for a in actual if a in predicted) / len(actual) if actual else 0
        print(f"  mod {m:2d}: residues={sorted(residues)}, "
              f"coverage={match:.2%}")
    print()

    # Connection to primorial structure
    # Primorial(k) = 2·3·5·...·p_k
    # Check if exceptions relate to primorial gaps
    primorials = [2, 6, 30, 210, 2310]
    for P in primorials:
        covered = [N for N in exceptions if N <= P]
        residues = set(N % P for N in exceptions)
        print(f"  Primorial {P}: {len(covered)} exceptions below it, "
              f"residues mod {P}: {sorted(residues)[:10]}...")
    print()


def representation_count_analysis(search_limit, tp_list):
    """How many ways can each even N be written as sum of two twin primes?"""
    print("=" * 70)
    print("PART C: Representation counts")
    print("=" * 70)
    print()

    tp_set = set(tp_list)
    counts = {}

    for N in range(2, min(search_limit, 500) + 1, 2):
        count = 0
        reps = []
        for p in tp_list:
            if p > N // 2:
                break
            q = N - p
            if q in tp_set:
                count += 1
                if len(reps) < 5:
                    reps.append((p, q))
        counts[N] = count

    # Show first values
    print("First even numbers and their representation counts:")
    for N in range(2, 102, 2):
        c = counts.get(N, 0)
        marker = " ***" if c == 0 else ""
        print(f"  {N:3d}: r(N) = {c:3d}{marker}")
    print()

    # Distribution of counts
    count_dist = Counter(counts.values())
    print(f"Distribution of r(N) for N <= {min(search_limit, 500)}:")
    for k in sorted(count_dist.keys())[:15]:
        print(f"  r(N)={k}: {count_dist[k]} values")
    print()

    # Growth of r(N)
    print("Growth of r(N):")
    for N in [50, 100, 200, 300, 400, 500]:
        if N in counts:
            print(f"  r({N}) = {counts[N]}")
    print()


def tournament_connection(exceptions):
    """Explore connections to tournament theory.

    Key question: do the exceptions relate to tournament sizes n where
    something special happens (e.g., A000568, H-spectrum, self-complementary)?

    Also: complement necklace = orbits of binary necklaces under
    complementation. For tournaments: T -> T^op. The number of
    complement-necklace classes is (A000568 + SC) / 2 = |G_n/Z_2|.
    """
    print("=" * 70)
    print("PART D: Tournament / complement necklace connections")
    print("=" * 70)
    print()

    # A000568 values (tournament isomorphism classes)
    A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}

    # Self-complementary counts
    SC = {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8, 7: 8, 8: 48}

    print("Complement necklaces = orbits under T -> T^op:")
    print("  |G_n/Z_2| = (A000568(n) + SC(n)) / 2")
    for n in sorted(A000568.keys()):
        merged = (A000568[n] + SC[n]) // 2
        print(f"  n={n}: A000568={A000568[n]}, SC={SC[n]}, "
              f"|G_n/Z_2|={merged}")
    print()

    # Check if any exception = 2 * tournament_size or similar
    exc_set = set(exceptions)
    print("Exceptions as 2*(n-1) for tournament sizes:")
    for n in range(2, 100):
        if 2 * (n - 1) in exc_set:
            print(f"  n={n}: 2(n-1)={2*(n-1)} IS an exception")
    print()

    # Check if exceptions are related to products of twin prime pairs
    tp_pairs = [(3, 5), (5, 7), (11, 13), (17, 19), (29, 31),
                (41, 43), (59, 61), (71, 73)]
    print("Twin pair products vs exceptions:")
    for p, q in tp_pairs:
        prod = p * q
        in_exc = prod in exc_set
        print(f"  {p}×{q}={prod}: {'EXCEPTION' if in_exc else 'representable'}")
    print()

    # Deeper: complement necklace of the binary representation
    # The exception values in binary — do they have complement structure?
    print("Binary representations of exceptions:")
    for N in exceptions[:20]:
        bits = bin(N)[2:]
        comp_bits = ''.join('1' if b == '0' else '0' for b in bits)
        comp_val = int(comp_bits, 2)
        in_exc = comp_val in exc_set
        print(f"  {N:4d} = {bits:>12s}  complement = {comp_bits} = {comp_val:4d}  "
              f"{'[also exception]' if in_exc else ''}")
    print()


def deeper_structure(exceptions):
    """Look for deeper number-theoretic structure."""
    print("=" * 70)
    print("PART E: Deeper structure")
    print("=" * 70)
    print()

    # Goldbach partition: can exceptions be written as p+q with regular primes?
    limit = max(exceptions) + 100
    is_prime = sieve_primes(limit)

    print("Goldbach partitions of exceptions (standard primes):")
    for N in exceptions:
        # Find the Goldbach partition
        for p in range(2, N):
            if is_prime[p] and is_prime[N - p]:
                q = N - p
                p_twin = is_prime.get(p - 2, False) or is_prime.get(p + 2, False) if p + 2 <= limit else False
                q_twin = is_prime.get(q - 2, False) or is_prime.get(q + 2, False) if q + 2 <= limit else False

                # Check properly
                p_is_twin = (p >= 4 and is_prime[p - 2]) or (p + 2 <= limit and is_prime[p + 2])
                q_is_twin = (q >= 4 and is_prime[q - 2]) or (q + 2 <= limit and is_prime[q + 2])

                print(f"  {N:4d} = {p} + {q}  "
                      f"({'twin' if p_is_twin else 'ISOLATED'}) + "
                      f"({'twin' if q_is_twin else 'ISOLATED'})")
                break
    print()

    # The key: every exception requires at least one ISOLATED prime
    # (a prime not in a twin pair)
    # Isolated primes: 2, 23, 37, 47, 53, 67, 79, 83, 89, 97, ...
    isolated = []
    for p in range(2, 200):
        if is_prime[p]:
            is_twin = (p >= 4 and is_prime[p - 2]) or (p + 2 <= limit and is_prime[p + 2])
            if not is_twin:
                isolated.append(p)

    print(f"Isolated primes (not in twin pair) up to 200: {isolated}")
    print(f"Count: {len(isolated)}")
    print()

    # The 2-structure: 2 is the only even prime and NOT a twin prime
    # (since 2+2=4 is not prime, and 2-2=0 is not prime)
    # Actually: is (2,3) a gap of 1, not 2? No, twin primes are (p, p+2).
    # 2 is only in a "cousin pair" (2,5) with gap 3, not a twin pair.
    # Wait: 3 and 5 are a twin pair (3,5). 2 is NOT part of any twin pair.
    # So 2 is isolated. It's the smallest isolated prime.
    print("Is 2 a twin prime? No — (2,4) fails since 4 is not prime.")
    print("Smallest twin primes: 3, 5 (the pair (3,5))")
    print()

    # Sum of two twin primes is always >= 3+3 = 6
    # So 2 and 4 are automatically exceptions.
    # What about N=6? 6 = 3+3 ✓ (both twin)
    print("Small cases:")
    for N in range(2, 20, 2):
        print(f"  {N}: {N in set(exceptions)}")
    print()

    # OEIS connection: this sequence should be in OEIS
    print("OEIS: search for this sequence")
    print(f"  Values: {exceptions[:15]}")
    print()


def main():
    print("Twin-Prime Goldbach Analysis — opus-2026-06-01-S520b")
    print()

    exceptions, tp_list = find_exceptions(search_limit=5000, prime_limit=100000)

    analyze_exceptions(exceptions, tp_list)
    complement_necklace_analysis(exceptions)
    representation_count_analysis(5000, tp_list)
    tournament_connection(exceptions)
    deeper_structure(exceptions)

    print("=" * 70)
    print("SYNTHESIS")
    print("=" * 70)
    print()
    print(f"Total exceptions: {len(exceptions)}")
    print(f"Largest exception: {max(exceptions) if exceptions else 'none'}")
    print()


if __name__ == "__main__":
    main()
