"""
tiling_count_deep.py — opus-2026-04-04-S9

Deep investigation of tiling counts per iso class.

KEY FINDINGS SO FAR:
- All tiling counts are ODD (at n=4,5,6)
- 7 and 21 are forbidden tiling counts at n≤6
- tiling_count = H(T) / |Aut(T)|

QUESTIONS:
1. Why are ALL tiling counts odd?
2. Is 7 permanently forbidden as a tiling count?
3. What is the structure of the forbidden tiling count set?
4. Connection to the cuboid model {2, 3, 7}?
5. The tiling count sequence: is it in OEIS?
"""
import numpy as np
from itertools import permutations
from collections import defaultdict, Counter
from math import gcd
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, make_tiles, tiling_to_tournament, count_hamiltonian_paths

# ==============================================================
# 1. WHY ARE ALL TILING COUNTS ODD?
# ==============================================================
def prove_tiling_counts_odd():
    """
    THEOREM: Every isomorphism class has an ODD number of tilings.

    PROOF IDEA: H(T) is always odd (Redei). |Aut(T)| is always odd for
    tournaments (because any automorphism of a tournament has odd order —
    a tournament on n vertices cannot have an automorphism of even order
    because an even-order automorphism would fix a vertex or create
    2-cycles which tournaments don't have).

    Wait — is |Aut(T)| always odd? For tournaments, the automorphism group
    has no elements of order 2 (since an involution would swap two vertices
    and their arc, but tournaments are antisymmetric so no 2-cycle is preserved).

    Actually, Moon's theorem: every tournament has an automorphism of odd order,
    but the GROUP ITSELF can have even order... let me check.

    A transposition (2-cycle) applied to a tournament reverses the arc between
    the two vertices. So it's an ANTI-automorphism, not an automorphism.
    Therefore Aut(T) has no elements of order 2.

    More generally: if σ has even order, it contains a transposition in its cycle
    decomposition... no, that's not right. σ can have even order via a 4-cycle etc.

    Let me just verify computationally.
    """
    print(f"{'='*60}")
    print(f"WHY ARE ALL TILING COUNTS ODD?")
    print(f"{'='*60}")

    for n in [4, 5, 6]:
        tiles = make_tiles(n)
        m = len(tiles)
        N = 1 << m

        # Compute |Aut| for each tournament
        aut_sizes = []
        for t in range(N):
            T = tiling_to_tournament(n, tiles, t)
            aut = 0
            for perm in permutations(range(n)):
                if all(T[perm[i]][perm[j]] == T[i][j] for i in range(n) for j in range(n)):
                    aut += 1
            aut_sizes.append(aut)

        # Check: all |Aut| odd?
        all_odd_aut = all(a % 2 == 1 for a in aut_sizes)
        aut_dist = Counter(aut_sizes)
        print(f"\n  n={n}: |Aut| distribution: {dict(sorted(aut_dist.items()))}")
        print(f"  All |Aut| odd? {all_odd_aut}")

        # Since H is odd (Redei) and |Aut| is odd, H/|Aut| is an odd/odd ratio.
        # For H/|Aut| to be an integer AND odd, we need |Aut| | H.

        # Verify: |Aut| divides H for all tournaments
        H, _, _ = compute_all_H(n)
        divides = all(int(H[t]) % aut_sizes[t] == 0 for t in range(N))
        print(f"  |Aut| divides H always? {divides}")

        # H/|Aut| always odd?
        quotients = [int(H[t]) // aut_sizes[t] for t in range(N)]
        all_odd_q = all(q % 2 == 1 for q in quotients)
        print(f"  H/|Aut| always odd? {all_odd_q}")

    print(f"""
  PROOF THAT TILING COUNT IS ALWAYS ODD:

  Step 1: H(T) is always odd (Redei's theorem).
  Step 2: |Aut(T)| is always odd for tournaments.
     Proof: An automorphism σ of tournament T satisfies T(σ(u), σ(v)) = T(u,v).
     If σ has a 2-cycle (u,v): then T(u,v) = T(σ(u), σ(v)) = T(v,u) = 1-T(u,v),
     which gives T(u,v) = 1/2, impossible for binary values.
     So σ has no 2-cycles in its cycle decomposition.
     More generally, any cycle of even length in σ would create a contradiction
     (the arc orientations around the cycle must be consistent but can't be
     for even-length cycles in a tournament).
     Therefore all cycles in σ have ODD length, so |σ| is odd.
     By Cauchy's theorem, if |Aut(T)| were even, there would be an element
     of order 2, which we just showed is impossible.
     Therefore |Aut(T)| is odd.

  Step 3: tiling_count = H(T)/|Aut(T)| = odd/odd.
     Since |Aut(T)| divides H(T) (orbit-stabilizer applied to the
     action of Aut(T) on Hamiltonian paths), the quotient is an integer.
     An odd number divided by an odd divisor is odd.
     Therefore the tiling count is always odd. QED.
    """)

# ==============================================================
# 2. IS 7 PERMANENTLY FORBIDDEN?
# ==============================================================
def analyze_forbidden_7():
    """
    tiling count = 7 requires H = 7 × |Aut|.
    Since |Aut| is always odd: need H ∈ {7, 21, 35, 49, 63, 77, 91, 105, ...}
    with the corresponding |Aut| ∈ {1, 3, 5, 7, 9, 11, 13, 15, ...}.

    H=7: PERMANENTLY FORBIDDEN (THM-029)
    H=21: PERMANENTLY FORBIDDEN (THM-079)
    H=35: achievable at n≥7. Need tournament with |Aut|=5.
    H=49: achievable at n≥7. Need tournament with |Aut|=7.

    Key question: does any tournament with H=35 have |Aut|=5?
    Or H=49 with |Aut|=7? Etc.
    """
    print(f"\n{'='*60}")
    print(f"IS TILING COUNT = 7 PERMANENTLY FORBIDDEN?")
    print(f"{'='*60}")

    # At n=7: compute H and |Aut| for ALL tilings
    n = 7
    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m
    H, _, _ = compute_all_H(n)

    # Target H values for tiling count = 7
    targets = [(7*k, k) for k in range(1, 30) if k % 2 == 1]

    print(f"  Need (H, |Aut|) pairs for tiling count = 7:")
    for h_target, aut_target in targets[:10]:
        count = sum(1 for t in range(N) if int(H[t]) == h_target)
        if count > 0:
            print(f"    H={h_target}, need |Aut|={aut_target}: "
                  f"{count} tilings with this H")
        elif h_target < 200:
            print(f"    H={h_target}: {'FORBIDDEN' if h_target in [7,21] else 'not achieved at n=7'}")

    # For H=35 at n=7: check which tournaments achieve it and their |Aut|
    print(f"\n  Detailed analysis for H=35 (need |Aut|=5 for tc=7):")
    h35_tilings = [t for t in range(N) if int(H[t]) == 35]
    print(f"    {len(h35_tilings)} tilings with H=35")

    if h35_tilings:
        # Compute |Aut| for a sample
        aut_counts_35 = []
        for t in h35_tilings[:20]:  # sample
            T = tiling_to_tournament(n, tiles, t)
            aut = sum(1 for perm in permutations(range(n))
                     if all(T[perm[i]][perm[j]] == T[i][j]
                           for i in range(n) for j in range(n)))
            aut_counts_35.append(aut)

        aut_dist_35 = Counter(aut_counts_35)
        print(f"    |Aut| distribution (sample of {len(aut_counts_35)}): "
              f"{dict(sorted(aut_dist_35.items()))}")
        has_aut5 = 5 in aut_dist_35
        print(f"    Any with |Aut|=5? {has_aut5}")

    # For H=49 at n=7: need |Aut|=7
    print(f"\n  Detailed analysis for H=49 (need |Aut|=7 for tc=7):")
    h49_tilings = [t for t in range(N) if int(H[t]) == 49]
    print(f"    {len(h49_tilings)} tilings with H=49")

    if h49_tilings:
        aut_counts_49 = []
        for t in h49_tilings[:20]:
            T = tiling_to_tournament(n, tiles, t)
            aut = sum(1 for perm in permutations(range(n))
                     if all(T[perm[i]][perm[j]] == T[i][j]
                           for i in range(n) for j in range(n)))
            aut_counts_49.append(aut)

        aut_dist_49 = Counter(aut_counts_49)
        print(f"    |Aut| distribution (sample of {len(aut_counts_49)}): "
              f"{dict(sorted(aut_dist_49.items()))}")
        has_aut7 = 7 in aut_dist_49
        print(f"    Any with |Aut|=7? {has_aut7}")

    # For H=105 at n=7: need |Aut|=15 for tc=7
    print(f"\n  Detailed analysis for H=105 (need |Aut|=15 for tc=7):")
    h105_tilings = [t for t in range(N) if int(H[t]) == 105]
    print(f"    {len(h105_tilings)} tilings with H=105")

    if h105_tilings:
        aut_counts_105 = []
        for t in h105_tilings[:10]:
            T = tiling_to_tournament(n, tiles, t)
            aut = sum(1 for perm in permutations(range(n))
                     if all(T[perm[i]][perm[j]] == T[i][j]
                           for i in range(n) for j in range(n)))
            aut_counts_105.append(aut)

        aut_dist_105 = Counter(aut_counts_105)
        print(f"    |Aut| distribution (sample): {dict(sorted(aut_dist_105.items()))}")

    # ALSO: What about tiling count = 21?
    print(f"\n  IS TILING COUNT = 21 PERMANENTLY FORBIDDEN?")
    targets_21 = [(21*k, k) for k in range(1, 20) if k % 2 == 1]
    for h_target, aut_target in targets_21[:8]:
        count = sum(1 for t in range(N) if int(H[t]) == h_target)
        if count > 0:
            print(f"    H={h_target}, need |Aut|={aut_target}: {count} tilings")
        elif h_target < 200:
            print(f"    H={h_target}: {'FORBIDDEN' if h_target in [7,21] else 'not achieved at n=7'}")

# ==============================================================
# 3. STRUCTURE OF THE TILING COUNT SET
# ==============================================================
def tiling_count_structure():
    """Analyze the set of achieved tiling counts and look for patterns."""
    print(f"\n{'='*60}")
    print(f"TILING COUNT SEQUENCE STRUCTURE")
    print(f"{'='*60}")

    for n in [5, 6]:
        tiles = make_tiles(n)
        m = len(tiles)
        N = 1 << m

        from forbidden_tiling_counts import tournament_canonical
        class_tilings = defaultdict(list)
        for t in range(N):
            T = tiling_to_tournament(n, tiles, t)
            canon = tournament_canonical(T, n)
            class_tilings[canon].append(t)

        tc_set = sorted(set(len(v) for v in class_tilings.values()))
        print(f"\n  n={n}: achieved tiling counts = {tc_set}")

        # Factor each tiling count
        for tc in tc_set:
            factors = []
            temp = tc
            for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
                while temp % p == 0:
                    factors.append(p)
                    temp //= p
            print(f"    tc={tc}: factorization = {factors if factors else [tc]}")

        # The achieved set modulo small primes
        mod3 = Counter(tc % 3 for tc in tc_set)
        mod7 = Counter(tc % 7 for tc in tc_set)
        print(f"\n    mod 3 distribution: {dict(sorted(mod3.items()))}")
        print(f"    mod 7 distribution: {dict(sorted(mod7.items()))}")

        # Gaps between consecutive achieved values
        gaps = [tc_set[i+1] - tc_set[i] for i in range(len(tc_set)-1)]
        print(f"    Gaps: {gaps}")

# ==============================================================
# 4. THE CUBOID CONNECTION
# ==============================================================
def cuboid_connection():
    """
    The cuboid model: Z/42Z = Z/2Z × Z/3Z × Z/7Z.
    Where do the tiling counts live in this cuboid?
    """
    print(f"\n{'='*60}")
    print(f"CUBOID CONNECTION: TILING COUNTS mod {2}, {3}, {7}")
    print(f"{'='*60}")

    for n in [5, 6]:
        tiles = make_tiles(n)
        m = len(tiles)
        N = 1 << m

        from forbidden_tiling_counts import tournament_canonical
        class_tilings = defaultdict(list)
        class_H = {}
        for t in range(N):
            T = tiling_to_tournament(n, tiles, t)
            canon = tournament_canonical(T, n)
            class_tilings[canon].append(t)
            if canon not in class_H:
                class_H[canon] = count_hamiltonian_paths(T, n)

        print(f"\n  n={n}:")
        print(f"  {'tc':>4s} {'H':>4s} {'tc%2':>4s} {'tc%3':>4s} {'tc%7':>4s} "
              f"{'H%3':>4s} {'H%7':>4s} {'tc%42':>5s}")

        for canon in sorted(class_tilings.keys(), key=lambda c: len(class_tilings[c])):
            tc = len(class_tilings[canon])
            h = class_H[canon]
            print(f"  {tc:4d} {h:4d} {tc%2:4d} {tc%3:4d} {tc%7:4d} "
                  f"{h%3:4d} {h%7:4d} {tc%42:5d}")

    # The KEY observation: tc%7 is never 0 at n≤6.
    # tc%7 = 0 would mean tc is divisible by 7.
    # Since tc = H/|Aut| and H is never divisible by 7 (because H%7 ≠ 0 at these n),
    # and |Aut| is not divisible by 7 (at n≤6, since Aut is a subgroup of S_n
    # and 7 > n for n≤6), we have tc = H/|Aut| with gcd(H, 7) = ? and gcd(|Aut|, 7) = 1.

    # So tc ≡ 0 mod 7 iff H ≡ 0 mod 7.
    # At n≤6: is H ever divisible by 7?
    for n in [5, 6]:
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)
        div7 = sum(1 for h in H if int(h) % 7 == 0)
        print(f"\n  n={n}: tilings with H divisible by 7: {div7}/{1<<m}")

if __name__ == "__main__":
    prove_tiling_counts_odd()
    analyze_forbidden_7()
    tiling_count_structure()
    cuboid_connection()
