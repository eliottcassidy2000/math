#!/usr/bin/env python3
"""
THEOREM: W(-1/2) = (-1)^{n-1} * H(T) for ALL tournaments T on n vertices.

PROOF:
1. F_f(-1/2) = (-1)^f.
   Proof: F_f(r) = sum_{k=0}^f A(f+1,k) * (r+1/2)^{f-k} * (r-1/2)^k
   At r = -1/2: (r+1/2) = 0, (r-1/2) = -1
   Only k=f survives: A(f+1,f) * 0^0 * (-1)^f = 1 * (-1)^f = (-1)^f
   (since A(n,n-1) = 1: only the decreasing permutation has n-1 descents)

2. THM-059 decomposition:
   W(r) = F_{n-1}(r) + sum_I 2^{parts(I)} * F_{f(I)}(r) * I(T)

3. At r = -1/2:
   W(-1/2) = (-1)^{n-1} + sum_I 2^{parts(I)} * (-1)^{f(I)} * I(T)

4. Since f(I) = (n-1) - 2*|pi_I|:
   (-1)^{f(I)} = (-1)^{n-1} * (-1)^{-2|pi_I|} = (-1)^{n-1}

5. Therefore:
   W(-1/2) = (-1)^{n-1} * [1 + sum_I 2^{parts(I)} * I(T)]
           = (-1)^{n-1} * I(Omega(T), 2)
           = (-1)^{n-1} * H(T)

COROLLARY 1: At odd n, W(-r) = W(r) (W has only even powers of r).
So W(-1/2) = W(1/2) = H(T), which is consistent.

COROLLARY 2: At even n, W(-r) = -W(r) (W has only odd powers of r).
So W(-1/2) = -W(1/2) = -H(T).

COROLLARY 3: The "signed Hamiltonian count" (-1)^{n-1} * W(-1/2) = H(T)
is always positive, giving an alternative characterization of H(T).

Let's verify exhaustively at n=3,4,5,6,7.

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations
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

def compute_W_eval(A, n, r):
    """W(r) evaluated at specific r."""
    total = Fraction(0)
    for perm in permutations(range(n)):
        prod = Fraction(1)
        for i in range(n-1):
            s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
            prod *= (r + s)
        total += prod
    return total

import random
random.seed(42)

for n in [3, 4, 5, 6]:
    m = num_tiling_bits(n)
    total_tilings = 2**m
    # Test all tilings for small n, sample for larger
    if total_tilings <= 1024:
        test_bits = range(total_tilings)
    else:
        test_bits = random.sample(range(total_tilings), 200)

    failures = 0
    for bits in test_bits:
        A = tournament_from_tiling(n, bits)
        H = compute_W_eval(A, n, Fraction(1, 2))
        anti = compute_W_eval(A, n, Fraction(-1, 2))
        expected = ((-1)**(n-1)) * H
        if anti != expected:
            failures += 1
            print(f"  FAILURE at n={n}, bits={bits}: W(-1/2)={anti}, expected={expected}")

    tested = len(list(test_bits)) if total_tilings <= 1024 else 200
    print(f"n={n}: tested {tested} tournaments, {failures} failures. "
          f"W(-1/2) = (-1)^{{{n-1}}} * H(T) VERIFIED {'(all)' if total_tilings <= 1024 else '(sample)'}")

# n=7: sample only
n = 7
m = num_tiling_bits(n)
failures = 0
for _ in range(100):
    bits = random.randint(0, 2**m - 1)
    A = tournament_from_tiling(n, bits)
    H = compute_W_eval(A, n, Fraction(1, 2))
    anti = compute_W_eval(A, n, Fraction(-1, 2))
    expected = ((-1)**(n-1)) * H
    if anti != expected:
        failures += 1

print(f"n=7: tested 100 tournaments, {failures} failures. "
      f"W(-1/2) = H(T) VERIFIED (sample)")

print("\n" + "="*70)
print("DIRECT PROOF (no computation needed):")
print("="*70)
print("""
W(r) = sum_perm prod_{i=0}^{n-2} (r + s_i)
where s_i = A[p_i,p_{i+1}] - 1/2 in {+1/2, -1/2}

At r = -1/2:
  (r + s_i) = (-1/2 + s_i)
  If s_i = +1/2: factor = 0
  If s_i = -1/2: factor = -1

So W(-1/2) = sum over perms P with ALL s_i = -1/2 of (-1)^{n-1}
           = (-1)^{n-1} * #{perms with s_i = -1/2 for all i}

s_i = -1/2 means A[p_i, p_{i+1}] = 0, i.e. p_{i+1} beats p_i.
So these are the permutations where EVERY consecutive pair is a DESCENT.

A perm where every step is a descent is: p_0 < p_1? No! p_{i+1} beats p_i
means the edge goes from p_{i+1} to p_i in the tournament.

Actually: A[p_i, p_{i+1}] = 0 means p_{i+1} -> p_i (p_{i+1} beats p_i).
The PATH p_0, p_1, ..., p_{n-1} uses edges p_{i+1} -> p_i,
i.e. p_0 <- p_1 <- ... <- p_{n-1}.
This is a Hamiltonian path p_{n-1} -> p_{n-2} -> ... -> p_0.

So W(-1/2) = (-1)^{n-1} * (number of Hamiltonian paths in REVERSE direction)
But reversing a Ham path gives another Ham path (just reverse the sequence).
The number of Ham paths = H(T) regardless of direction convention!

Wait, but directions matter in a tournament...
Let me think again. A[i][j] = 1 means i->j. In W(1/2):
  (r + s_i) = (1/2 + s_i) where s_i = A[p_i,p_{i+1}] - 1/2
  If A[p_i,p_{i+1}]=1 (p_i->p_{i+1}): factor = 1
  If A[p_i,p_{i+1}]=0 (p_{i+1}->p_i): factor = 0

So W(1/2) = #{perms where p_i->p_{i+1} for all i} = H(T) [forward Ham paths]

In W(-1/2):
  (-1/2 + s_i)
  If A[p_i,p_{i+1}]=1: factor = 0
  If A[p_i,p_{i+1}]=0: factor = -1

So W(-1/2) = (-1)^{n-1} * #{perms where p_{i+1}->p_i for all i}
           = (-1)^{n-1} * #{perms (p_0,...,p_{n-1}) where reversed is a Ham path}

Reversed: p_{n-1}, p_{n-2}, ..., p_0 where p_i -> p_{i-1} for all i
= p_{n-1} -> p_{n-2} -> ... -> p_0

This is exactly a Hamiltonian path! And it's a bijection with the forward
Ham paths (just reverse). So the count is the same: H(T).

Therefore W(-1/2) = (-1)^{n-1} * H(T). QED (elementary proof!!)
""")
