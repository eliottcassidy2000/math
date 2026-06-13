#!/usr/bin/env python3
"""
case_b_proof_s243b.py -- opus-2026-03-23-S243b

PROVING Case B = (n-2)!

Strategy: Count the total Fix_B sum and show it equals n! × (n-2)!.

For a permutation σ with σ(a) = b, σ(b) = a:
σ = (a b) ∘ τ where τ ∈ S_{[n]\{a,b}} (permutation of the other n-2 vertices).

For σ to fix the edge {T, T^{(a,b)}}:
σ(T) = T^{(a,b)} means: for all i,j:
  T[σ^{-1}(i)][σ^{-1}(j)] = T^{(a,b)}[i][j]

Since T^{(a,b)} agrees with T everywhere except at {a,b}:
  For (i,j) ≠ (a,b) and (i,j) ≠ (b,a): T[σ^{-1}(i)][σ^{-1}(j)] = T[i][j]
  For (i,j) = (a,b): T[σ^{-1}(a)][σ^{-1}(b)] = 1 - T[a][b]
    i.e., T[b][σ^{-1}(b)] = 1 - T[a][b] (since σ^{-1}(a) = b)

Wait, let me be more careful. σ^{-1}(a) = b, σ^{-1}(b) = a (since σ swaps a,b).

For i=a, j=b: T[σ^{-1}(a)][σ^{-1}(b)] = T[b][a] = 1 - T[a][b].
And T^{(a,b)}[a][b] = 1 - T[a][b].
So: 1 - T[a][b] = 1 - T[a][b]. ✓ AUTOMATICALLY SATISFIED!

For i=a, j=c (c ≠ a,b): T[σ^{-1}(a)][σ^{-1}(c)] = T[b][τ^{-1}(c)] = T[a][c]
So: T[b][τ^{-1}(c)] = T[a][c] for all c ≠ a,b.

For i=b, j=c (c ≠ a,b): T[σ^{-1}(b)][σ^{-1}(c)] = T[a][τ^{-1}(c)] = T[b][c]
So: T[a][τ^{-1}(c)] = T[b][c] for all c ≠ a,b.

For i=c, j=d (c,d ≠ a,b): T[τ^{-1}(c)][τ^{-1}(d)] = T[c][d]
So: τ is an automorphism of T restricted to [n]\{a,b}.

SUMMARY: σ = (a b) ∘ τ maps T to T^{(a,b)} iff:
1. τ is an automorphism of T|_{[n]\{a,b}} (the induced sub-tournament)
2. T[b][τ^{-1}(c)] = T[a][c] for all c ∉ {a,b}
   Equivalently: the out-neighborhoods of a and b, restricted to [n]\{a,b},
   are τ-related: τ maps the out-neighborhood of b to the out-neighborhood of a.

Condition 2 says: τ is a bijection [n]\{a,b} → [n]\{a,b} that maps
{c : b→c in T} to {c : a→c in T}. Combined with condition 1 (τ ∈ Aut(T|_{rest})),
this is a strong constraint.

THE COUNT:
sum over all T, all pairs {a,b}, all τ:
  [τ ∈ Aut(T|_{rest})] × [τ maps N_out(b) to N_out(a)]

Now I'll count this by SWITCHING THE ORDER OF SUMMATION.

Fix a pair {a,b} and a permutation τ ∈ S_{[n]\{a,b}}.
Count the number of tournaments T such that:
1. τ ∈ Aut(T|_{[n]\{a,b}})
2. τ maps N_out(b,T) to N_out(a,T) restricted to [n]\{a,b}

For τ to be an automorphism of T|_{rest}: the sub-tournament on [n]\{a,b}
must be fixed by τ. The number of such sub-tournaments = 2^{free_orbits(τ)}
where free_orbits(τ) is the number of free edge orbits under τ on [n]\{a,b}.

For condition 2: given the sub-tournament (fixed by τ), we need to choose
orientations of arcs from a and b to each c ∈ [n]\{a,b} such that:
  T[b][τ^{-1}(c)] = T[a][c] for all c.

This means: for each τ-orbit {c, τ(c), τ²(c), ...} of [n]\{a,b}:
  T[a][c] = T[b][τ^{-1}(c)] = T[b][τ^{-1}(c)]
  T[a][τ(c)] = T[b][c]
  etc.

So the arc from a to τ^k(c) determines the arc from b to τ^{k-1}(c).

For an orbit of length L: choosing T[a][c] for the first element determines
all T[a][τ^j(c)] and T[b][τ^j(c)] IF the orbit has odd length.
If even length: there's a consistency condition.

Actually, let me just count:
- Arcs between {a} and rest: n-2 arcs from a, n-2 arcs from b. Total 2(n-2).
- Constraint: T[a][c] and T[b][τ^{-1}(c)] linked. So (n-2) free choices.
  Wait: T[a][c] determines T[b][τ^{-1}(c)]. Since τ is a bijection,
  as c ranges over [n]\{a,b}, so does τ^{-1}(c).
  So: T[b][d] = T[a][τ(d)] for all d ∈ [n]\{a,b}.

  This means: the arcs from b are COMPLETELY DETERMINED by the arcs from a
  (via τ). So we have n-2 free binary choices (arcs from a to each other vertex).

- Arc between a and b: 1 arc, free binary choice. But wait — the edge_swap
  condition doesn't constrain the a↔b arc (we showed it's automatically satisfied).
  So 1 more free choice.

- Arcs within [n]\{a,b}: must be fixed by τ. Free choices = 2^{orbits(τ on C(n-2,2) pairs)}.

Total tournaments satisfying both conditions:
  = 2 × 2^{n-2} × 2^{orbits(τ)}
  where the 2 is for the a↔b arc, 2^{n-2} for arcs from a to rest.

Wait: n-2 arcs from a to each of the n-2 other vertices = 2^{n-2} choices.
Then arcs from b are determined (n-2 arcs, no free choice).
Arc a↔b: 2 choices (1 free bit).
Arcs within rest: 2^{orbits(τ)} choices.

Total = 2 × 2^{n-2} × 2^{orbits(τ on (n-2)-vertex pairs)}

Now sum over all T, {a,b}, τ:
= sum_{pairs {a,b}} sum_{τ ∈ S_{n-2}} [2^{1+(n-2)+orbits(τ)}]
= C(n,2) × sum_{τ ∈ S_{n-2}} 2^{n-1+orbits(τ)}

And: sum_{τ ∈ S_{n-2}} 2^{orbits(τ)} = (n-2)! × V_{n-2}
where V_{n-2} = A000568(n-2) = #{iso classes of (n-2)-tournaments}.

Wait, that's Burnside: sum_{τ} 2^{orbits(τ)} = (n-2)! × V_{n-2}.

So: total Fix_B sum = C(n,2) × 2^{n-1} × (n-2)! × V_{n-2}

And Case B = total / n! = C(n,2) × 2^{n-1} × (n-2)! × V_{n-2} / n!
= [n(n-1)/2] × 2^{n-1} × (n-2)! × V_{n-2} / n!
= 2^{n-1} × V_{n-2} × (n-2)! × (n-1) / (2 × (n-1)!)
= 2^{n-1} × V_{n-2} × (n-2)! / (2 × (n-2)!)
= 2^{n-2} × V_{n-2}

So: Case B = 2^{n-2} × V_{n-2}?!

Let me check:
  n=3: 2^1 × V_1 = 2 × 1 = 2. But Case B should be (n-2)! = 1! = 1. FAIL!

Hmm, there's an error. Let me recount.

Actually V_1 = 1 (one tournament on 1 vertex).
2^{n-2} × V_{n-2} = 2 × 1 = 2 ≠ 1.

Let me recheck the counting.

Author: opus-2026-03-23-S243b
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

banner("DETAILED FIX_B COUNTING")

# For each n, for each pair {a,b} and each τ ∈ S_{rest},
# count the number of T such that σ=(ab)∘τ maps T to T^{(a,b)}.

V_tournament = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56}

for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs_n = [(i,j) for i in range(n) for j in range(i+1,n)]
    pair_idx = {p: k for k, p in enumerate(pairs_n)}

    total_count = 0  # Total #{(T, {a,b}, τ) : σ=(ab)τ maps T to T^{(a,b)}}

    for a in range(n):
        for b in range(a+1, n):
            rest = [c for c in range(n) if c != a and c != b]
            nr = len(rest)  # n-2

            # For each τ ∈ S_{rest}
            for tau_perm in permutations(rest):
                tau = {rest[i]: tau_perm[i] for i in range(nr)}
                tau[a] = a; tau[b] = b  # extend to all of [n] (only rest moves)

                # Full σ: swap a↔b, apply τ to rest
                sigma = {i: tau[i] for i in range(n)}
                sigma[a] = b; sigma[b] = a

                # Count T such that σ(T) = T^{(a,b)}
                # σ(T)[i][j] = T[σ^{-1}(i)][σ^{-1}(j)]
                # We need σ(T)[i][j] = T^{(a,b)}[i][j] for all i,j

                # The arcs split into:
                # 1. Arc (a,b): automatic
                # 2. Arcs (a,c) for c in rest: constrained by T[b][τ^{-1}(c)] = T[a][c]
                # 3. Arcs (b,c) for c in rest: determined by #2
                # 4. Arcs (c,d) for c,d in rest: T must be τ-automorphic on rest

                # Count free bits:
                # Arc (a,b): 2 choices (0 or 1)
                # For condition 2: T[a][c] for each c ∈ rest: n-2 free bits
                #   Then T[b][d] = T[a][τ(d)] for all d ∈ rest (determined)
                # For condition 4: T|_{rest} must be fixed by τ
                #   Free bits = #{orbits of τ on pairs of rest}

                # Count orbits of τ on pairs of rest
                rest_pairs = [(rest[i], rest[j]) for i in range(nr) for j in range(i+1, nr)]
                seen = set()
                n_orbits = 0
                for (p, q) in rest_pairs:
                    if (p, q) in seen: continue
                    # Follow the orbit
                    pp, qq = p, q
                    while True:
                        if pp < qq:
                            seen.add((pp, qq))
                        else:
                            seen.add((qq, pp))
                        pp, qq = tau[pp], tau[qq]
                        if pp < qq:
                            key = (pp, qq)
                        else:
                            key = (qq, pp)
                        if key in seen: break
                    n_orbits += 1

                # Total T for this (a,b,τ): 2 × 2^{n-2} × 2^{n_orbits}
                count_T = 2 * (2**(nr)) * (2**n_orbits)
                total_count += count_T

    # Expected: n! × (n-2)! (for Case B orbits = (n-2)!)
    expected = factorial(n) * factorial(n-2)

    # But wait — we're counting each EDGE twice in Case B
    # (once as σ fixing {T, T^k}, once as σ fixing {T^k, T})
    # Actually no: each (T, {a,b}, τ) triple is counted once.
    # The Case B sum is: sum_σ |Fix_B(σ)|
    # = sum_σ #{edges {T, T^k} σ swaps}
    # For each σ = (ab)τ: we need to count edges {T, T^{(a,b)}}
    # where σ(T) = T^{(a,b)}. Each such T gives one edge.
    # But if T = T^{(a,b)} (self-loop), it's NOT an edge to swap.

    # Let's verify: the total_count should include self-loops.
    # We need to subtract self-loop cases: T = T^{(a,b)}.

    # For now, let me just compare.

    print("  n=%d: total_count = %d, n!*(n-2)! = %d, ratio = %.4f" %
          (n, total_count, expected, total_count / expected))

    # Hmm, let me also verify by direct edge counting
    # Fix_B(σ) = #{edges {T, T^k} : σ swaps them, T ≠ T^k}

    # Direct computation:
    direct_fix_b = 0
    for perm_tuple in permutations(range(n)):
        perm = {i: perm_tuple[i] for i in range(n)}
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

        for bits in range(2**m):
            for k in range(m):
                flipped = bits ^ (1<<k)
                if flipped <= bits: continue
                pb = apply_perm(bits)
                pf = apply_perm(flipped)
                if pb == flipped and pf == bits:
                    direct_fix_b += 1

    print("    Direct Fix_B sum = %d" % direct_fix_b)
    print("    Expected (n!*(n-2)!) = %d" % expected)
    print("    total_count / direct = %.4f" % (total_count / direct_fix_b if direct_fix_b > 0 else 0))
    print()

banner("THE PROOF STATUS")

print("""
  CASE A: PROVED. Case A = T_n/2 by Burnside on (T, arc) pairs.

  CASE B: The formula Case B = (n-2)! is VERIFIED at n=3,4,5.
  The proof would show: sum_σ Fix_B(σ) = n! × (n-2)!

  The counting argument above gives a product formula for the
  number of tournaments T compatible with each (a,b,τ).
  Summing over all (a,b) and τ should give n! × (n-2)!.

  If the total from the counting formula matches direct computation,
  the proof reduces to summing 2^{1+(n-2)+orbits(τ)} over all
  pairs and τ, which is a cycle-index computation.
""")

print("Done. opus-2026-03-23-S243b")
