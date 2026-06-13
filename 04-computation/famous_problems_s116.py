#!/usr/bin/env python3
"""famous_problems_s116.py — Famous problems through the Cayley lens"""
from math import sqrt, log, pi

phi = (1+sqrt(5))/2

print("FAMOUS PROBLEMS THROUGH THE CAYLEY LENS")
print("="*60)

print("""
1. RIEMANN HYPOTHESIS
  zeta(s) = sum Q(x_n)^{-s} over Cayley addresses.
  RH: Re(s)=1/2 means amplitude at each address = exp(-arctanh(x_n)) = n^{-1/2}.
  The critical line IS the hyperbolic midpoint: where amplitude
  equals the geometric mean of the Cayley transform.

  DOES IT HELP? The geometric picture is clean but probably not
  sufficient for a proof. The RH is about ZEROS of an entire function,
  and our tools are about COEFFICIENTS and WEIGHTS.

2. TWIN PRIME CONJECTURE
  Twin primes (p, p+2) have arctanh gap ln(1+2/p)/2 ~ 1/p -> 0.
  The conjecture: infinitely many primes with another prime within
  arctanh distance 1/p (vanishingly close in hyperbolic space).

  Compare: Bertrand gives primes within distance ln(2)/2 = 0.347 (constant).
  Twin primes: distance ~ 1/p (vanishing). MUCH stronger.
""")

# Twin prime arctanh gaps
print("  Twin prime arctanh gaps:")
for p in [3, 5, 11, 29, 41, 101, 1049, 10007]:
    p2 = p + 2
    is_prime_p = all(p % d != 0 for d in range(2, int(sqrt(p))+1))
    is_prime_p2 = all(p2 % d != 0 for d in range(2, int(sqrt(p2))+1))
    if is_prime_p and is_prime_p2:
        gap = log((p+2)/p)/2
        print(f"    ({p}, {p+2}): arctanh gap = {gap:.6f} ~ 1/{p} = {1/p:.6f}")

print("""
3. GOLDBACH CONJECTURE
  Goldbach: for each even 2k, there exist primes p,q with p+q=2k.
  In Cayley: Q(x_p) + Q(x_q) = 2k.
  This is a DIOPHANTINE equation in the Cayley transform of prime addresses.

  The difficulty: addition is NON-HYPERBOLIC on the Cayley line.
  The line is built for multiplication. Goldbach is additive.
  This incommensurability may be the geometric reason Goldbach is hard.

4. P vs NP (the most promising connection)
  Computing H(T) is #P-hard.
  Our formula: H(T)/E[H] ~ 1 + sum c_k / n^k (Cayley-Delannoy expansion).

  Truncating at level K gives H(T) with error O(1/n^K).
  Each level k requires O(n^{2k}) time to compute.
  For fixed K: polynomial time, O(1/n^K) multiplicative error.

  CONJECTURE: For MOST tournaments, truncating the Cayley-Delannoy
  expansion at level K gives H(T) to within (1+O(1/n^K)) factor.
  This would be a PTAS for a #P-hard problem.

  CAVEAT: The expansion is for the VARIANCE over random tournaments.
  It does NOT directly give a per-tournament approximation.
  The per-tournament version would need concentration inequalities
  showing individual H(T) is close to the ensemble average
  at each Fourier level. This is plausible but UNPROVED.

5. COLLATZ CONJECTURE
  Already explored: Collatz as random walk on Cayley line.
  Even step: -log(2)/2. Odd step: +log(3)/2.
  Drift: -0.048/step (leftward). But drift is not enough for a proof.

  The Cayley perspective gives EXACTLY the well-known heuristic.
  It does not go beyond it.

6. ABC CONJECTURE (proved by Mochizuki, controversial)
  For a+b=c with gcd(a,b)=1: rad(abc)^{1+eps} > c for all eps>0.
  rad(n) = product of distinct prime factors of n.

  On the Cayley line: rad(n) corresponds to the PRODUCT of prime
  addresses appearing in n (without multiplicity). The abc conjecture
  says: the sum a+b=c (additive) is CONTROLLED by the multiplicative
  structure (rad = prime address product).

  This is EXACTLY the tension between addition and multiplication
  on the Cayley line that we identified. ABC says: addition cannot
  stray too far from what multiplication allows.
""")

print("="*60)
print("THE HONEST VERDICT")
print("="*60)
print("""
The Cayley framework gives CLEAN REFORMULATIONS but does not
directly SOLVE any famous problem. The three most promising:

1. P vs NP / #P: The Cayley-Delannoy expansion might give a
   polynomial-time approximation scheme for H(T) counting.
   Needs: per-tournament concentration at each Fourier level.
   Status: CONJECTURE, testable computationally.

2. Twin primes: The hyperbolic gap formulation is clean:
   can prime addresses be arbitrarily close in arctanh-space?
   But this is just a restatement, not a proof technique.

3. ABC: The addition/multiplication tension on the Cayley line
   is EXACTLY what ABC is about. But Cayley tools are multiplicative
   and ABC needs additive control, which is the hard part.

The Cayley framework is a LANGUAGE for stating problems, not
a METHOD for solving them. But sometimes the right language
reveals the right approach. The per-tournament PTAS conjecture
is the most likely to yield something concrete.
""")
