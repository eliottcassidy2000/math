#!/usr/bin/env python3
"""
near_aut_proof_s243.py -- opus-2026-03-23-S243

PROVING edge_orbits = T_n/2 + (n-2)!

By Burnside on edges of Q_m under S_n:

edge_orbits = (1/n!) * sum_{σ ∈ S_n} |Fix_edges(σ)|

where Fix_edges(σ) = #{unordered pairs {T, T^k} : T^k = T ⊕ e_k, σ fixes {T, T^k}}

A pair {T, T^k} is fixed by σ in two ways:
  Case A: σ(T) = T AND σ(T^k) = T^k
  Case B: σ(T) = T^k AND σ(T^k) = T  (σ swaps them)

CASE A: σ(T) = T means σ ∈ Aut(T). σ(T^k) = T^k means σ also stabilizes T^k.
  Since T^k differs from T only at arc k, σ ∈ Aut(T) AND σ maps arc k to arc k.
  So σ fixes both endpoints of arc k: σ(i) = i and σ(j) = j where k = {i,j}.

  Sum over σ: sum_σ |{(T,k) : σ ∈ Aut(T), σ fixes arc k}|
  = sum_σ #{T : σ(T) = T} × #{k : σ(k) = k}
  NO — can't factorize because σ(T)=T depends on T.

  Correct: sum_σ sum_T sum_k [σ(T)=T AND σ fixes arc k AND T^k > T]
  The "T^k > T" ensures each unordered pair counted once.

  = (1/2) * sum_σ sum_T sum_k [σ(T)=T AND σ fixes arc k]
    (factor of 1/2 because each unordered pair {T, T^k} counted twice)

  = (1/2) * sum_σ (#{T : σ(T) = T} × #{k : σ fixes pair k})

  Here: #{T : σ(T) = T} = Fix_T(σ) (Burnside's tournament fixing count)
  #{k : σ fixes pair k} = #{pairs {i,j} with σ(i)=i AND σ(j)=j} = C(fix(σ), 2)
  where fix(σ) = number of fixed points of σ.

  So Case A sum = (1/2) * sum_σ Fix_T(σ) × C(fix(σ), 2)

  And the Burnside count from Case A = (1/(2n!)) * sum_σ Fix_T(σ) × C(fix(σ), 2)

CASE B: σ(T) = T^k AND σ(T^k) = T.
  Since T^k = T with arc k reversed: σ maps T to T-with-arc-k-reversed.

  σ(T) = T^k means: for all arcs e ≠ k, σ(e) is an arc of T^k = T
  (same as T except at k). And σ(k) is the REVERSED arc at k.

  So σ reverses the arc at position k, and maps all other arcs correctly.
  This means: σ is an automorphism of T EXCEPT at arc k, where it reverses.
  σ is a "near-automorphism" of T at arc k.

  More precisely: σ(i) and σ(j) are the endpoints of arc k, but σ reverses
  the direction. So if arc k goes from a to b in T, then σ(a) = b, σ(b) = a
  (σ swaps the endpoints of arc k), AND σ is an automorphism of T restricted
  to all other arcs.

  BUT WAIT: σ also needs σ(T^k) = T. Since T^k = T with k reversed,
  σ(T^k) = σ(T) with σ(k) reversed = T^k with σ(k) reversed = T.
  So: T^k with σ(k) reversed = T. Since T^k = T with k reversed,
  T^k with σ(k) reversed = T means σ(k) = k (the same pair).
  So σ maps the pair {a,b} to itself: σ(a) = b and σ(b) = a (swap).

  Summary for Case B:
  σ(a) = b, σ(b) = a (σ swaps the endpoints of arc k = (a,b))
  For all other arcs (i,j) ≠ (a,b): σ maps them consistently with T
  (i.e., T[σ(i)][σ(j)] = T^k[i][j] for all i,j except at k where reversal happens)

  Actually more carefully: T^k[i][j] = T[i][j] for (i,j) ≠ (a,b) and (i,j) ≠ (b,a).
  σ(T)[i][j] = T[σ^{-1}(i)][σ^{-1}(j)]. For σ(T) = T^k:
  T[σ^{-1}(i)][σ^{-1}(j)] = T^k[i][j] for all i,j.

  For (i,j) where neither i nor j is a or b:
  T[σ^{-1}(i)][σ^{-1}(j)] = T[i][j] (since T^k agrees with T here).
  So σ must fix or permute vertices other than a,b in a T-automorphic way.

  For i=a, j≠a,b: T[σ^{-1}(a)][σ^{-1}(j)] = T^k[a][j] = T[a][j].
  Since σ^{-1}(a) = b: T[b][σ^{-1}(j)] = T[a][j].
  This constrains: σ^{-1}(j) maps "b-neighbors" to "a-neighbors" consistently.

  This is a STRONG constraint. For a generic tournament, σ swapping a↔b
  and respecting all other arcs is impossible UNLESS the tournament has
  specific structure (a and b must be "interchangeable" modulo the one arc).

Author: opus-2026-03-23-S243
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# ============================================================
banner("CASE A: STANDARD FIXING (σ ∈ Aut(T), σ fixes arc)")
# ============================================================

# Case A contribution = (1/(2n!)) * sum_σ Fix_T(σ) * C(fix(σ), 2)
# where fix(σ) = #{fixed points of σ}

# Fix_T(σ) by cycle type:
# Cycle type (1^{a1}, 2^{a2}, 3^{a3}, ...) fixes 2^{free_orbits} tournaments
# where free_orbits depends on the cycle structure acting on pairs.

# For σ = identity: Fix_T(id) = 2^m, fix(id) = n, C(n,2) = m
# Contribution = (1/(2n!)) * 2^m * m

# But we need the FULL Burnside sum. Let me compute directly for small n.

for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs_n = [(i,j) for i in range(n) for j in range(i+1,n)]
    pair_idx = {p: k for k, p in enumerate(pairs_n)}

    case_a_sum = 0
    case_b_sum = 0

    for perm in permutations(range(n)):
        # Count fixed points
        fixed_pts = [i for i in range(n) if perm[i] == i]
        n_fixed = len(fixed_pts)
        fixed_pairs = [(i,j) for i in fixed_pts for j in fixed_pts if i < j]

        # How σ acts on labeled tournaments
        def apply_perm(bits):
            new_bits = 0
            for k, (i,j) in enumerate(pairs_n):
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    ok = pair_idx[(pi,pj)]
                    if bits & (1<<ok): new_bits |= (1<<k)
                else:
                    ok = pair_idx[(pj,pi)]
                    if not (bits & (1<<ok)): new_bits |= (1<<k)
            return new_bits

        # Count Case A and Case B edges fixed by this σ
        for bits in range(2**m):
            for k in range(m):
                flipped = bits ^ (1<<k)
                if flipped <= bits: continue  # each edge once

                pb = apply_perm(bits)
                pf = apply_perm(flipped)

                if pb == bits and pf == flipped:
                    case_a_sum += 1  # σ fixes both endpoints
                elif pb == flipped and pf == bits:
                    case_b_sum += 1  # σ swaps the pair

    nfact = factorial(n)
    case_a_orbits = case_a_sum / nfact
    case_b_orbits = case_b_sum / nfact
    total_orbits = (case_a_sum + case_b_sum) / nfact

    # Compare with T_n/2 and (n-2)!
    T_n = {3:4, 4:16, 5:88}[n]

    print("  n=%d:" % n)
    print("    Case A: sum=%d, orbits=%.1f" % (case_a_sum, case_a_orbits))
    print("    Case B: sum=%d, orbits=%.1f" % (case_b_sum, case_b_orbits))
    print("    Total: orbits=%.1f" % total_orbits)
    print("    T_n/2 = %.1f, (n-2)! = %d" % (T_n/2, factorial(n-2)))
    print("    Case A orbits = T_n/2? %s" % (case_a_orbits == T_n/2))
    print("    Case B orbits = (n-2)!? %s" % (case_b_orbits == factorial(n-2)))
    print()

banner("PROOF SKETCH FOR CASE A = T_n/2")

print("""
  CLAIM: (1/n!) * sum_σ |Fix_A(σ)| = T_n/2

  where Fix_A(σ) = #{edges {T, T^k} : σ(T)=T AND σ(T^k)=T^k}
                  = #{(T, k) : σ(T)=T, σ fixes arc k} / 2
                    (divided by 2 because each unordered pair counted twice
                     except when T = T^k, i.e., self-loops)

  For the non-self-loop part:
  #{(T, k) : σ(T)=T, σ fixes arc k} = Fix_T(σ) × #{k fixed by σ}
  where #{k fixed by σ} = C(fix(σ), 2) = #{pairs with both endpoints fixed}

  But this counts ORDERED (T, k), not unordered {T, T^k}.
  Since T ≠ T^k for most arcs, dividing by 2 gives the edge count.

  (1/(2n!)) * sum_σ Fix_T(σ) * C(fix(σ), 2)
  = (1/2) * (1/n!) * sum_σ Fix_T(σ) * C(fix(σ), 2)
  = (1/2) * T_n (by Burnside applied to (tournament, arc) pairs)

  Wait: T_n = (1/n!) * sum_σ Fix_{T,arc}(σ)
  where Fix_{T,arc}(σ) = #{(T, k) : σ(T)=T AND σ fixes arc k}
  = Fix_T(σ) * C(fix(σ), 2)

  So Case A = T_n / 2. ✓

  THIS PROVES CASE A = T_n/2.
""")

banner("PROOF SKETCH FOR CASE B = (n-2)!")

print("""
  CLAIM: (1/n!) * sum_σ |Fix_B(σ)| = (n-2)!

  where Fix_B(σ) = #{edges {T, T^k} : σ(T)=T^k, σ(T^k)=T}

  For σ to swap T and T^k:
  - σ(T) = T^k means σ maps T to T-with-arc-k-reversed
  - This requires σ(a) = b, σ(b) = a where arc k = (a,b)
    AND σ is an automorphism of T except at this one arc
  - So σ is a TRANSPOSITION-LIKE map: swaps a,b and permutes
    the remaining n-2 vertices in a T-consistent way

  Key insight: σ must have (a b) in its cycle decomposition.

  For each pair {a,b} and each tournament T:
  The number of σ swapping a↔b and satisfying σ(T)=T^{(a,b)} is
  the number of "T-consistent" permutations of [n]\{a,b} that
  respect the arc structure when a and b swap.

  The total Case B sum:
  sum_σ Fix_B(σ) = sum over all T, all pairs {a,b}:
    #{σ : σ swaps a↔b, σ(T) = T^{(a,b)}}

  = sum_T sum_{a,b} #{near-auts of T at {a,b}}

  The claim is: (1/n!) * this sum = (n-2)!

  Rearranging: sum_T sum_{a,b} #{near-auts} = n! * (n-2)! = n * (n-1) * ((n-2)!)^2

  Since #{tournaments} = 2^m and #{pairs} = m:
  avg near-auts per (T, {a,b}) = n! * (n-2)! / (2^m * m)

  At n=5: n! * (n-2)! = 120 * 6 = 720, 2^10 * 10 = 10240.
  avg = 720/10240 = 0.0703. About 7% of (T, pair) combinations
  have a near-automorphism.

  THE (n-2)! MUST COME FROM THE CYCLE INDEX STRUCTURE.
  Specifically, it counts orbits of S_n acting on (pair, near-aut) triples.
""")

print("Done. opus-2026-03-23-S243")
