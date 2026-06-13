#!/usr/bin/env python3
"""
var_H_formula_derivation.py — kind-pasteur-2026-03-21-S16

DERIVING A FORMULA FOR Var(H) OF RANDOM TOURNAMENTS.

Key identity: E[H^2] = sum_{sigma,tau: compatible} 2^{-|pairs(sigma) cup pairs(tau)|}

where compatible means no pair {u,v} has sigma using (u,v) and tau using (v,u).

The sum depends on the OVERLAP structure of permutation pairs.

CLASSIFICATION of permutation pairs by overlap type:
For sigma, tau in S_n, define:
  A = #{undirected pairs used by both in SAME direction}
  C = #{undirected pairs used by both in OPPOSITE directions}
  S = #{undirected pairs used by sigma only}
  T = #{undirected pairs used by tau only}

Then: |pairs(sigma)| = n-1, |pairs(tau)| = n-1
      A + C + S = n-1, A + C + T = n-1 (so S = T)
      |union| = A + C + S + T = A + C + 2S = 2(n-1) - A + C
      Wait: |union| = A + S + T (since C-direction pairs have BOTH directions,
      but in a tournament only one is present, so each C pair uses ONE undirected pair)

Actually: pairs are UNDIRECTED. sigma uses {sigma(k), sigma(k+1)} for each k.
|pairs(sigma)| = n-1 (all distinct since it's a path).

  agree = # pairs {u,v} used by both sigma and tau with same direction
  conflict = # pairs {u,v} used by both with opposite directions
  |pairs(sigma) ∩ pairs(tau)| = agree + conflict
  |pairs(sigma) ∪ pairs(tau)| = 2(n-1) - agree - conflict = 2(n-1) - overlap

  P(both HPs) = 0 if conflict > 0
  P(both HPs) = 2^{-|union|} = 2^{-(2(n-1) - agree)} if conflict = 0

So E[H^2] = sum_{sigma,tau: conflict=0} 2^{agree - 2(n-1)}
           = 2^{-2(n-1)} * sum_{sigma,tau: conflict=0} 2^{agree}

This is a sum over compatible permutation pairs weighted by 2^{agree}.

For sigma = tau: agree = n-1, contribution = 2^{n-1} per perm, total = n! * 2^{n-1}
For sigma != tau: agree < n-1.

Let f(n, a) = #{ordered pairs (sigma, tau) with agree=a and conflict=0}

Then E[H^2] = 2^{-2(n-1)} * sum_a f(n,a) * 2^a

Let me compute f(n,a) for small n.
"""
from itertools import permutations
from collections import Counter
from fractions import Fraction

for n in range(3, 8):
    print(f"\n{'='*60}")
    print(f"  n = {n}")
    print(f"{'='*60}")

    perms = list(permutations(range(n)))

    # For each pair, compute agree and conflict
    agree_conflict_count = Counter()  # (agree, conflict) -> count

    for sigma in perms:
        # Directed arcs of sigma
        sigma_arcs = {}  # undirected pair -> direction (ordered pair)
        for k in range(n-1):
            pair = frozenset({sigma[k], sigma[k+1]})
            sigma_arcs[pair] = (sigma[k], sigma[k+1])

        for tau in perms:
            agree = 0
            conflict = 0
            for k in range(n-1):
                pair = frozenset({tau[k], tau[k+1]})
                if pair in sigma_arcs:
                    if sigma_arcs[pair] == (tau[k], tau[k+1]):
                        agree += 1
                    else:
                        conflict += 1
            agree_conflict_count[(agree, conflict)] += 1

    # Compatible pairs: conflict = 0
    print("  (agree, conflict) distribution:")
    for (a, c) in sorted(agree_conflict_count.keys()):
        cnt = agree_conflict_count[(a, c)]
        print(f"    agree={a}, conflict={c}: count={cnt}")

    # f(n, a) for conflict=0
    print(f"\n  f({n}, a) for compatible pairs:")
    total_compatible = 0
    weighted_sum = 0
    for a in range(n):
        f_a = agree_conflict_count.get((a, 0), 0)
        if f_a > 0:
            contrib = f_a * (2**a)
            total_compatible += f_a
            weighted_sum += contrib
            print(f"    f({n}, {a}) = {f_a:10d}, 2^{a}*f = {contrib:12d}")

    E_H2_num = weighted_sum
    E_H2_den = 4**(n-1)
    E_H2 = Fraction(E_H2_num, E_H2_den)

    E_H = Fraction(len(perms), 2**(n-1))
    Var_H = E_H2 - E_H**2

    print(f"\n  E[H^2] = {E_H2_num} / {E_H2_den} = {E_H2}")
    print(f"  E[H] = {E_H}")
    print(f"  Var(H) = {Var_H}")

    # Factor the Var(H) numerator
    p = Var_H.numerator
    factors = []
    pp = p
    for f in range(2, 1000):
        while pp % f == 0:
            factors.append(f)
            pp //= f
    if pp > 1:
        factors.append(pp)
    print(f"  Var(H) = {Var_H.numerator}/{Var_H.denominator}")
    print(f"  Numerator factors: {factors}")

    # The agree distribution for conflict=0 gives the key structure
    print(f"\n  Compatible pair count: {total_compatible} / {len(perms)**2}")
    print(f"  = {Fraction(total_compatible, len(perms)**2)}")
