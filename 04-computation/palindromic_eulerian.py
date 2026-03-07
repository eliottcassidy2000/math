#!/usr/bin/env python3
"""
Tournament Eulerian distribution palindromy — is this a new theorem?

For a tournament T, define a_k(T) = #{perms P : exactly k forward edges in P}.
"Forward edge" means A[P(i), P(i+1)] = 1.

OBSERVATION: a_k(T) = a_{n-1-k}(T) for ALL tournaments T at n=5.

This is equivalent to: the number of perms with k forward edges equals
the number with n-1-k forward edges.

Proof attempt: Given a perm P = (p_0, ..., p_{n-1}), its REVERSE
P' = (p_{n-1}, ..., p_0) has the property that step i in P' uses edge
(p_{n-1-i}, p_{n-i}), which is the reverse of step (n-1-i) in P.

In P: step j uses edge (p_j, p_{j+1}), which is forward iff A[p_j, p_{j+1}] = 1.
In P': step i uses edge (p_{n-1-i}, p_{n-2-i}), which is forward iff
A[p_{n-1-i}, p_{n-2-i}] = 1, i.e., iff step j=n-2-i in P is BACKWARD.

So #forward(P') = #backward(P) = (n-1) - #forward(P).

This means reversal maps k-forward perms to (n-1-k)-forward perms,
and it's a bijection (reversal is an involution on S_n).

Therefore a_k(T) = a_{n-1-k}(T) for ALL tournaments and ALL n.
THIS IS A THEOREM (elementary!).

Let's verify at n=4,6,7 and derive consequences.

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations, combinations
from collections import defaultdict
from fractions import Fraction

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

import random
random.seed(42)

for n in [3, 4, 5, 6]:
    m = num_tiling_bits(n)
    total = 2**m
    test_set = range(total) if total <= 1024 else random.sample(range(total), 200)

    palindromic = 0
    non_palindromic = 0
    for bits in test_set:
        A = tournament_from_tiling(n, bits)
        forward_counts = defaultdict(int)
        for perm in permutations(range(n)):
            fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
            forward_counts[fwd] += 1

        is_pal = all(forward_counts.get(k, 0) == forward_counts.get(n-1-k, 0) for k in range(n))
        if is_pal:
            palindromic += 1
        else:
            non_palindromic += 1
            print(f"  FAILURE at n={n}, bits={bits}: {dict(sorted(forward_counts.items()))}")

    tested = len(list(test_set)) if total <= 1024 else 200
    print(f"n={n}: {palindromic}/{tested} palindromic, {non_palindromic} failures")

# Now: consequences of palindromy for W(r)
print(f"\n{'='*70}")
print("CONSEQUENCES OF PALINDROMY")
print("="*70)
print("""
W(r) = sum_k a_k * (r+1/2)^k * (r-1/2)^{n-1-k}

With a_k = a_{n-1-k}, let p = r+1/2, q = r-1/2. Then:
W(r) = sum_k a_k * p^k * q^{n-1-k}

Since a_k = a_{n-1-k}:
W(r) = sum_k a_{n-1-k} * p^k * q^{n-1-k}

Substituting j = n-1-k:
= sum_j a_j * p^{n-1-j} * q^j

Compare with:
W(-r) = sum_k a_k * (-r+1/2)^k * (-r-1/2)^{n-1-k}
      = sum_k a_k * (-q)^k * (-p)^{n-1-k}
      = (-1)^{n-1} sum_k a_k * q^k * p^{n-1-k}
      = (-1)^{n-1} sum_k a_k * p^{n-1-k} * q^k

So W(-r) = (-1)^{n-1} * sum_k a_k * p^{n-1-k} * q^k

And from palindromy: sum_k a_k * p^{n-1-k} * q^k = W(r)
(by the substitution above with a_k = a_{n-1-k}).

Wait, that gives W(-r) = (-1)^{n-1} * W(r).

THIS IS THE r-PARITY OF W(r)! (THM-059 property (iii))

So the palindromy of the forward-edge distribution IS EQUIVALENT to the
r-parity of W(r)! Both follow from the same reversal bijection.

The reversal proof is actually the SIMPLEST proof of the r-parity property.
""")

# Average forward edge distribution
print(f"\n{'='*70}")
print("AVERAGE FORWARD EDGE DISTRIBUTION")
print("="*70)
for n in [3, 4, 5]:
    m = num_tiling_bits(n)
    avg_dist = defaultdict(lambda: Fraction(0))
    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        for perm in permutations(range(n)):
            fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
            avg_dist[fwd] += Fraction(1, 2**m)

    print(f"  n={n}: avg forward edge distribution (per tournament):")
    for k in range(n):
        print(f"    k={k}: {avg_dist[k]} = {float(avg_dist[k]):.4f}")
    print(f"    Sum = {sum(avg_dist.values())} = {n}!")

# What IS the average W polynomial?
print(f"\n{'='*70}")
print("AVERAGE W(r) = E_T[W_T(r)]")
print("="*70)
for n in [3, 4, 5]:
    m = num_tiling_bits(n)
    avg_W = defaultdict(lambda: Fraction(0))
    for bits in range(2**m):
        A = tournament_from_tiling(n, bits)
        result = defaultdict(lambda: Fraction(0))
        for perm in permutations(range(n)):
            poly = {0: Fraction(1)}
            for i in range(n-1):
                s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
                new_poly = {}
                for power, coeff in poly.items():
                    new_poly[power+1] = new_poly.get(power+1, Fraction(0)) + coeff
                    new_poly[power] = new_poly.get(power, Fraction(0)) + coeff * s
                poly = new_poly
            for power, coeff in poly.items():
                result[power] += coeff
        for p, c in result.items():
            avg_W[p] += c / 2**m

    terms = []
    for p in sorted(avg_W.keys(), reverse=True):
        if avg_W[p] != 0:
            terms.append(f"{avg_W[p]}*r^{p}")
    print(f"  n={n}: avg W(r) = {' + '.join(terms)}")

    # Compare with F_{n-1}(r) + invariant contributions
    # Average t3 = C(n,3)/4 (each triple has 1/4 chance of each direction)
    # Wait: average t3 = C(n,3)/4 since for random tournament, each triple
    # has prob 1/4 of being a 3-cycle (either direction)
    from math import comb
    avg_t3 = Fraction(comb(n, 3), 4)
    print(f"    avg t3 = C({n},3)/4 = {avg_t3}")
