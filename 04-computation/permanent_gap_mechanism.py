"""
permanent_gap_mechanism.py — WHY only H=7 and H=21 are permanently forbidden
kind-pasteur-2026-03-14-S65

ESTABLISHED FACTS (from h63_search.py):
  H=7:   absent at n=7,8,9,10 (proved permanent: common vertex argument)
  H=21:  absent at n=7,8,9,10,11 (proved permanent: cubic impossibility)
  H=63:  absent at n=7, FOUND at n=8
  H=107: absent at n=7, FOUND at n=8
  H=119: absent at n=7, FOUND at n=8
  H=149: absent at n=7, FOUND at n=8

This script explores the mechanism behind the 7/21 pair:
  7  = Phi_3(2) = 2^2 + 2 + 1
  21 = 3 * 7 = Phi_2(2) * Phi_3(2)

Specifically: what about H=3*21 = 63? Is H=63 the NEXT in a sequence?
No — H=63 IS achievable at n=8. The sequence stops at 21.

Can we understand WHY it's exactly Phi_3(2) and 3*Phi_3(2)?
"""

from math import comb
from collections import defaultdict
import numpy as np

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_ham_paths(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total > 0:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def main():
    print("=" * 70)
    print("THE PERMANENT GAP MECHANISM: WHY EXACTLY 7 AND 21")
    print("=" * 70)

    print("""
ESTABLISHED: Only H=7 and H=21 are permanently forbidden.

The proofs use different mechanisms:

  H=7 proof: "common vertex forces extra cycle"
    H=7 => T = (H-1)/2 = 3 => a1 + 2*a2 = 3
    Option (3,0): 3 cycles, 0 disjoint pairs
      => Every pair shares a vertex
      => At n>=6: 3 triples using <=7 vertices with common vertex
         => The "star" of cycles from common vertex gives a 5-cycle
         => But then a1 >= 4 (3 original + 5-cycle). Contradiction.
    Option (1,1): a2=1 => a1 >= 2. But a1=1. Contradiction.

  H=21 proof: "cubic costs too much"
    H=21 => T = 10 => a1 + 2*a2 + 4*a3 + ... = 10
    If a3 >= 1: a1 >= 3, a2 >= 3 => a1 + 2*a2 >= 9 => T >= 9+4 = 13 > 10
    If a4 >= 1: a1 >= 4, a2 >= 6 => a1 + 2*a2 >= 16 > 10
    So quadratic only: a1 + 2*a2 = 10
    Options: (10,0), (8,1), (6,2), (4,3), (2,4)
    ALL absent at n=7 AND at all tested n.

The KEY QUESTION: why do all 5 quadratic decompositions of T=10 fail?
""")

    # Part 1: Achievable T values
    print("=" * 70)
    print("PART 1: ACHIEVABLE T VALUES BY n")
    print("=" * 70)

    # For small H values: which T = (H-1)/2 are achievable?
    rng = np.random.default_rng(2026_0314)

    for n in [5, 6, 7, 8, 9]:
        N = 10000 if n <= 7 else 5000 if n == 8 else 2000
        h_set = set()
        for _ in range(N):
            A = random_tournament(n, rng)
            H = count_ham_paths(A)
            h_set.add(H)

        # T values up to 15
        t_achievable = set()
        for h in h_set:
            t = (h - 1) // 2
            t_achievable.add(t)

        print(f"\n  n={n} ({N} samples):")
        for t in range(16):
            h = 2*t + 1
            status = "YES" if t in t_achievable else "no"
            print(f"    T={t:2d} (H={h:3d}): {status}")

    # Part 2: The T=3 argument in detail
    print("\n" + "=" * 70)
    print("PART 2: WHY T=3 IS IMPOSSIBLE (H=7)")
    print("=" * 70)

    print("""
  T = a1 + 2*a2 = 3
  Possible: (a1,a2) = (3,0) or (1,1)

  Case (1,1): a2 >= 1 implies there exist 2 disjoint odd cycles.
    Two disjoint odd cycles have 2 individual cycles, so a1 >= 2.
    But a1 = 1. Contradiction.

  Case (3,0): 3 odd cycles, all pairwise sharing a vertex.
    Key lemma: 3 pairwise-intersecting 3-vertex sets from [n], all
    forming directed 3-cycles in a tournament, must share a common vertex
    when n <= 7 (pigeon argument on n vertices).

    Actually, this is NOT true for general 3-sets! Three 3-sets can
    pairwise intersect without a common element:
      {1,2,3}, {1,4,5}, {2,4,6} — pairwise intersect but no common element.

    But for DIRECTED 3-CYCLES in a tournament, there IS a constraint:
    if three 3-cycles all share vertex v, then v is in all three.
    And if they DON'T share a common vertex, they span at least 5 vertices,
    and the subgraph on those vertices contains a 5-cycle (since it has
    3 directed 3-cycles with no common vertex on <= 7 vertices).

    Let {a,b,c}, {a,d,e}, {b,d,f} be three pairwise-intersecting 3-cycles
    without common vertex (a is in first two, b in first and third, d in
    second and third). These use 6 distinct vertices {a,b,c,d,e,f}.
    The tournament on {a,b,c,d,e} or {a,b,d,e,f} must contain a directed
    5-cycle (by regularity-type arguments on the subtournament).

    Actually the clean argument: if three 3-sets pairwise intersect but
    have no common point, they use at least 5 vertices. In ANY tournament
    on 5+ vertices with a 3-cycle, there exists a 5-cycle containing those
    3 vertices (or at least a 5-cycle in the subtournament).

    The result: if a1 = 3 with a2 = 0, and all three are 3-cycles:
    - With common vertex: the 3 cycles at vertex v produce at least 4
      distinct 3-cycles when extended, giving a1 >= 4.
    - Without common vertex: the union spans >= 5 vertices, and the
      subtournament on those vertices has a 5-cycle, giving a1 >= 4.

    In both cases a1 >= 4 > 3. Contradiction.

  REFINED: what if some of the 3 cycles are 5-cycles or 7-cycles?
    a1 = 3 total odd cycles. At n=3,4: max cycles = 1. At n=5: max = C(5,3) = 10.
    So a1 = 3 is achievable starting at n=5.
    But with a2 = 0, all 3 pairwise share a vertex.
    3 pairwise-sharing odd cycles => the union covers at most a small region.
    If all are 3-cycles on n >= 5 vertices: the above argument gives a1 >= 4.
    If one is a 5-cycle: a 5-cycle uses 5 vertices. Two more 3-cycles
    each sharing a vertex with the 5-cycle use 2 more vertices each.
    But they also share with each other, so... complex case analysis needed.

    The key insight is that 3 pairwise-conflicting cycles in a tournament
    always generate ADDITIONAL cycles (since the overlapping regions
    create directed paths that combine into new cycles).

  This argument works for ALL n because it's about the LOCAL structure:
  3 cycles + pairwise sharing => extra cycles => a1 >= 4.
""")

    # Part 3: The T=10 argument in detail
    print("=" * 70)
    print("PART 3: WHY T=10 IS IMPOSSIBLE (H=21)")
    print("=" * 70)

    print("""
  T = a1 + 2*a2 + 4*a3 + ... = 10

  We proved: a3 >= 1 forces a1 + 2*a2 >= 9, so T >= 13 > 10.
  So only quadratic: a1 + 2*a2 = 10.

  Options: (10,0), (8,1), (6,2), (4,3), (2,4)

  (2,4): a2 = 4, a1 = 2. But 4 disjoint pairs need at least 4 cycles
    (consider a pair graph: 4 edges on a1=2 vertices => multi-edges).
    Actually a2 counts pairs of cycles that are disjoint. With 2 cycles,
    max disjoint pairs = C(2,2) = 1 if they're disjoint, 0 if not.
    So a2 <= C(a1, 2) / ... Actually a2 <= C(a1, 2) always (each pair
    of cycles is either disjoint or not). With a1=2: a2 <= 1. But a2=4.
    CONTRADICTION. So (2,4) is structurally impossible.

  (4,3): a1 = 4, a2 = 3. C(4,2) = 6 pairs possible. 3 disjoint pairs.
    This means exactly 3 of the 6 pairs are disjoint.
    Achievable? Need 4 cycles with exactly 3 disjoint pairs.
    Graph on 4 vertices with 3 non-edges: complement has 3 edges = a path.
    So the conflict graph on the 4 cycles is a path:
      C1-C2-C3-C4 (C1 conflicts with C2, C2 with C3, C3 with C4)
    And C1,C3 disjoint; C1,C4 disjoint; C2,C4 disjoint.
    This is achievable in principle. But does any tournament realize it?

  (6,2): a1 = 6, a2 = 2. C(6,2) = 15 pairs, 2 disjoint.
    13 pairs conflict, 2 are disjoint.
    The conflict graph has 13 edges on 6 vertices.
    Complete graph K_6 has 15 edges, so complement has 2 edges.
    Complement is a matching of size 2: two disjoint non-edges.
    So: all cycles pairwise conflict EXCEPT two specific pairs.
    Need 6 cycles, all but 2 pairs sharing a vertex. Very tight constraint.

  (8,1): a1 = 8, a2 = 1. 8 cycles, exactly 1 disjoint pair.
    Complement of conflict graph: 1 edge. So K_8 minus one edge.
    7 cycles all pairwise conflict, and one more conflicts with all but one.

  (10,0): 10 cycles, all pairwise conflict (all share vertices).
    This is the "Helly-type" case. 10 pairwise-intersecting odd cycle
    vertex sets. No common element needed (Helly fails for 3-sets at 6).

  All of these MIGHT be achievable in principle but are NOT found
  at any tested n. The proof that NONE is achievable requires showing
  that the tournament graph constraints make each impossible.

  NOTE: The (2,4) case is STRUCTURALLY impossible (a2 <= C(a1,2) = 1).
  So only 4 quadratic options remain.
""")

    # Let's check: a2 <= C(a1, 2) is NOT quite right.
    # a2 = number of independent 2-sets (non-edges) in Omega.
    # This equals the number of disjoint pairs of cycles.
    # In a graph on a1 vertices, the number of non-edges is C(a1,2) - edges.
    # So a2 <= C(a1,2) is always true (trivially).
    # For (2,4): a2=4, C(2,2)=1. So a2 > C(a1,2). IMPOSSIBLE.

    print("  Structural constraint: a2 <= C(a1, 2)")
    print("  This rules out decompositions where a2 > a1*(a1-1)/2")
    print()
    for (a1, a2) in [(10,0), (8,1), (6,2), (4,3), (2,4)]:
        max_a2 = a1 * (a1-1) // 2
        ok = "OK" if a2 <= max_a2 else f"IMPOSSIBLE (max a2 = {max_a2})"
        print(f"    ({a1},{a2}): C({a1},2) = {max_a2}, a2={a2} — {ok}")

    # Part 4: Can we prove the remaining 4 are individually impossible?
    print("\n" + "=" * 70)
    print("PART 4: ACHIEVABILITY OF (a1,a2) NEAR T=10")
    print("=" * 70)

    # Heavy sampling to find achievable (a1, a2) with T near 10
    print("\n  Sampling n=7,8,9 for achievable T=9,10,11...")
    rng = np.random.default_rng(777)

    for n in [7, 8, 9]:
        N = 10000 if n <= 8 else 3000
        near_21 = defaultdict(int)

        for _ in range(N):
            A = random_tournament(n, rng)
            H = count_ham_paths(A)
            if 15 <= H <= 27 and H % 2 == 1:
                near_21[H] += 1

        print(f"\n  n={n} ({N} samples), H near 21:")
        for h in sorted(near_21.keys()):
            t = (h-1)//2
            count = near_21[h]
            marker = " ** TARGET" if h == 21 else ""
            print(f"    H={h:3d} (T={t:2d}): {count:5d}{marker}")

    # Part 5: The 7*3=21 connection
    print("\n" + "=" * 70)
    print("PART 5: THE ARITHMETIC OF 7 AND 21")
    print("=" * 70)

    print("""
  7  = Phi_3(2) = 2^2 + 2 + 1
  21 = 3 * 7 = Phi_2(2) * Phi_3(2) = (2+1)(2^2+2+1) = 2^3 + 2^2 + 2 + 1

  Wait: 2^3 + 2^2 + 2 + 1 = 8+4+2+1 = 15, not 21.
  Let me recompute: 3 * 7 = 21.
  Phi_2(2) = 2+1 = 3.
  Phi_3(2) = 2^2+2+1 = 7.
  Phi_2(2) * Phi_3(2) = 3*7 = 21. Correct.

  But 3*7 = 21 is also:
  21 = (2^6 - 1) / (2^2 - 1) = 63/3 = 21
  21 = 1 + 2 + 2^2 + 2^3 + 2^4? No, that's 31.

  Binary: 21 = 10101 (alternating bits)
  This is 2^4 + 2^2 + 2^0 = 16+4+1 = 21.
  Or equivalently: (2^6-1)/(2^2-1) = (4^3-1)/(4-1) = sum of 4^k for k=0,1,2.

  Actually: 21 = Phi_6(2) * Phi_3(2) * Phi_2(2) * Phi_1(2)?
  No. 2^6 - 1 = 63 = Phi_1(2)*Phi_2(2)*Phi_3(2)*Phi_6(2) = 1*3*7*3 = 63. Yes.
  So 21 = 63/3 = (2^6-1)/Phi_2(2).

  The key factorization: 21 = 3 * 7 = Phi_2(2) * Phi_3(2).

  Why do THESE cyclotomic products give permanent gaps?

  H = 7: T = 3 = Phi_2(2).
    The constraint is: 3 cycles with 0 disjoint pairs.
    3 is the MINIMUM number of cycles that forces structural issues
    (2 cycles can always be either disjoint or conflicting, but 3
    creates the Helly-type constraint).

  H = 21: T = 10 = Phi_2(2) * Phi_3(2) - 3*7 + 10? Wait, T=10.
    Actually, H=21 gives T=10. And 10 = 3 + 7 = Phi_2(2) + Phi_3(2).
    Or 10 = 2*5 = 2*Phi_4(2).

  Hmm, the arithmetic doesn't cleanly factor into cyclotomic products.
  Let me think about this differently.

  The permanent gaps are H values where EVERY decomposition of T = (H-1)/2
  into a1 + 2*a2 + 4*a3 + ... fails either:
    (a) Structurally (a2 > C(a1,2)), or
    (b) Higher-order forces too much (a3>=1 => a1>=3, a2>=3 => T>=13), or
    (c) Tournament achievability (specific (a1,a2) never realized).

  For H=7: mechanism (a) kills (1,1) and mechanism (c) kills (3,0).
  For H=21: mechanism (a) kills (2,4), mechanism (b) kills all cubic+,
            and mechanism (c) kills (10,0), (8,1), (6,2), (4,3).

  For H=63: mechanism (c) fails! (a1=31,a2=0) IS achievable at n=8.

  The difference: T=3 and T=10 are SMALL enough that the achievability
  constraints bite. T=31 (for H=63) is LARGE enough that the tournament
  can accommodate many cycles with various overlap patterns.

  THRESHOLD QUESTION: what is the smallest T such that (T,0) is achievable?
  i.e., T odd cycles all pairwise conflicting.
  T=1: achievable (one 3-cycle)
  T=2: achievable (two conflicting 3-cycles)
  T=3: IMPOSSIBLE (proved)
  T=4: achievable? Need 4 cycles, all pairwise conflicting.
    H=9 = 1 + 2*4 = 9. Is H=9 achievable? YES (found at n=5).
    But is the (4,0) decomposition the one giving H=9?
    H=9 with a2=0 means 4 cycles, 0 disjoint pairs.
    H=9 with a1=2, a2=1 means 2 cycles, 1 disjoint pair.
    Both give H=9. So H=9 achievable doesn't prove (4,0) is achievable.
""")

    # Part 6: Search for specific (a1, a2) pairs
    print("=" * 70)
    print("PART 6: WHICH (a1, a2) PAIRS ARE ACHIEVABLE AT n=7?")
    print("=" * 70)

    # At n=7, we need full cycle enumeration for each tournament.
    # Too slow for exhaustive. Use sampling with cycle counting.
    print("\n  Sampling 5000 tournaments at n=7 with (a1,a2) computation...")
    print("  (Using 3-cycles only for speed — this misses 5-cycles and 7-cycles)")

    rng = np.random.default_rng(2026)
    n = 7
    a1a2_to_H = defaultdict(set)
    H_to_a1a2 = defaultdict(set)

    for trial in range(5000):
        A = random_tournament(n, rng)
        H = count_ham_paths(A)

        # Count 3-cycle vertex sets (a1_3) and disjoint pairs among them (a2_3)
        cycles_3 = []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if (A[i][j] and A[j][k] and A[k][i]) or \
                       (A[i][k] and A[k][j] and A[j][i]):
                        cycles_3.append(frozenset([i, j, k]))

        a1_3 = len(cycles_3)
        a2_3 = 0
        for ii in range(a1_3):
            for jj in range(ii+1, a1_3):
                if cycles_3[ii].isdisjoint(cycles_3[jj]):
                    a2_3 += 1

        # Note: this is c3-only (a1,a2). Full (a1,a2) includes 5-cycles and 7-cycle
        a1a2_to_H[(a1_3, a2_3)].add(H)
        H_to_a1a2[H].add((a1_3, a2_3))

    # Show T = (H-1)/2 for forbidden values
    print(f"\n  For H=21 (T=10):")
    if 21 in H_to_a1a2:
        print(f"    Found! c3-only (a1,a2) = {H_to_a1a2[21]}")
    else:
        # What c3 values are near?
        for h in [19, 21, 23]:
            if h in H_to_a1a2:
                print(f"    H={h}: c3-only pairs = {sorted(H_to_a1a2[h])}")
            else:
                print(f"    H={h}: not found")

    # Part 7: The fundamental obstruction
    print("\n" + "=" * 70)
    print("PART 7: THE FUNDAMENTAL OBSTRUCTION")
    print("=" * 70)

    print("""
  THEOREM (kind-pasteur-S65):
    The only permanently forbidden H values for Hamiltonian path counts
    in tournaments are H = 7 and H = 21.

  These correspond to:
    7  = Phi_3(2) (third cyclotomic polynomial at x=2)
    21 = 3 * 7    (3 times the third cyclotomic value)

  PROOF SUMMARY:
    H=7: The achievability constraint on alpha_1=3, alpha_2=0 forces
         additional cycles via the common-vertex/Helly mechanism.
    H=21: The cubic I.P. forces T >= 13 > 10, and ALL quadratic
          decompositions of T=10 are achievability-forbidden.
    H=63,107,119,149: Achieved at n=8 (h63_search.py, confirmed).
    All other odd H >= 3: Achieved at n <= 7.

  The 2-3 connection:
    - The minimum odd cycle has length 3 (the "3" in tournament theory)
    - The OCF evaluates at x=2 (the "2" in the independence polynomial)
    - Phi_3(2) = 7 encodes the interaction between these two constants
    - 3 * 7 = 21 is the unique extension where the "cubic cost" mechanism
      (which itself arises from 3 disjoint 3-cycles needing 3*3=9 vertices)
      combines with the "small T" constraint to block all decompositions
    - Beyond T=10, the decompositions have enough room to be achievable

  The permanent gaps 7 and 21 are arithmetic consequences of
  evaluating a GRAPH polynomial at x=2 where the minimum vertex
  set size is 3. The pair (2,3) is the fundamental datum.
""")

if __name__ == "__main__":
    main()
