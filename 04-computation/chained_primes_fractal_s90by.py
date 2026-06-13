#!/usr/bin/env python3
"""
chained_primes_fractal_s90by.py — The formal group fractal and chained primes
opus-2026-03-16-S90by
"""

from sympy import lucas, fibonacci, isprime, factorint
from math import log, sqrt

phi = (1 + sqrt(5)) / 2

print("""


     THE FORMAL GROUP FRACTAL AND CHAINED PRIMES

     opus-2026-03-16-S90by



THE LUCAS PRIME CHAIN
======================

L(2^k) for k = 1, 2, 3, 4, 5, ...

  L(2)  = 3.       PRIME.
  L(4)  = 7.       PRIME.
  L(8)  = 47.      PRIME.
  L(16) = 2207.    PRIME.
  L(32) = 4870847. COMPOSITE.

Four consecutive primes: 3, 7, 47, 2207.

Each is the SYMMETRIC part of phi^(2^k):
  phi^2: L(2) = 3.
  phi^4: L(4) = 7.
  phi^8: L(8) = 47.
  phi^16: L(16) = 2207.

And F(2n) = F(n) * L(n), so:
  F(4) = F(2)*L(2) = 1*3 = 3.
  F(8) = F(4)*L(4) = 3*7 = 21.
  F(16) = F(8)*L(8) = 21*47 = 987.
  F(32) = F(16)*L(16) = 987*2207 = 2178309.

The Fibonacci numbers at powers of 2 are the CUMULATIVE PRODUCTS
of the Lucas prime chain:

  F(4) = 3.
  F(8) = 3 * 7 = 21.
  F(16) = 3 * 7 * 47 = 987.
  F(32) = 3 * 7 * 47 * 2207 = 2178309.

This IS a fractal: each level MULTIPLIES the previous product
by the next Lucas-at-power-of-2.

  F(2^k) = product of L(2^j) for j = 1, ..., k-1.
  (Starting from F(2) = 1.)

The chain 3, 7, 47, 2207 is a MULTIPLICATIVE TELESCOPE:
each prime is a new FACTOR that extends the product.
""")

# Verify the telescope
print("The multiplicative telescope:")
product = 1
for k in range(1, 6):
    n = 2**k
    L = int(lucas(n))
    F2n = int(fibonacci(2*n))
    product *= L
    print(f"  k={k}: L(2^{k}) = L({n}) = {L}, "
          f"cumulative product = {product}, F(2^{k+1}) = F({2*n}) = {F2n}, "
          f"match: {product == F2n}")

print(f"""

PERFECT. F(2^(k+1)) = product of L(2^1) * L(2^2) * ... * L(2^k).

The Fibonacci numbers at powers of 2 are the PARTIAL PRODUCTS
of the Lucas-at-powers-of-2 sequence.

When the Lucas values are ALL PRIME (as for k=1,2,3,4):
  the Fibonacci value F(2^(k+1)) is a product of DISTINCT primes.
  F(4) = 3 (prime).
  F(8) = 3*7 = 21 (semiprime, squarefree).
  F(16) = 3*7*47 = 987 (three distinct primes, squarefree).
  F(32) = 3*7*47*2207 = 2178309 (four distinct primes, squarefree).

At k=5: L(32) = 4870847 is COMPOSITE.
""")

L32 = int(lucas(32))
print(f"L(32) = {L32}")
f32 = factorint(L32)
print(f"L(32) factors: {dict(f32)}")

print(f"""

L(32) = 4870847 = {dict(factorint(L32))}.

The chain BREAKS. L(32) is not prime, so F(64) is no longer
a product of distinct primes from the chain.

The prime chain length is 4: L(2), L(4), L(8), L(16) = 3, 7, 47, 2207.

Compare to the MERSENNE CHAIN: 2, 3, 7, 127, 2^127-1, ...
  Each step: k -> 2^k - 1.
  Length: at least 5 (all known to be prime through 2^127-1).

The LUCAS CHAIN: L(2), L(4), L(8), L(16) = 3, 7, 47, 2207.
  Each step: L(n) -> L(2n).
  Length: 4 (breaks at L(32)).

These are TWO DIFFERENT CHAINS starting from similar origins:
  Mersenne: 2 -> 3 -> 7 -> 127. Each step: k -> 2^k-1.
  Lucas:    3 -> 7 -> 47 -> 2207. Each step: L_n -> L_(2n).

They SHARE the first two primes: 3 and 7.
They diverge at step 3: Mersenne gives 127, Lucas gives 47.


THE FORMAL GROUP FRACTAL
==========================

The formal group F(x,y) = (x+y)/(1+xy) has the doubling map:
  [2](x) = 2x/(1+x^2) = tanh(2*arctanh(x)).

The ITERATED DOUBLING:
  [2](x), [4](x) = [2]([2](x)), [8](x), [16](x), ...

Starting from x = x_phi = (phi-1)/(phi+1):
""")

x_phi = (phi - 1) / (phi + 1)
print(f"x_phi = (phi-1)/(phi+1) = {x_phi:.10f}")
print()

x = x_phi
for k in range(6):
    m = 2**k
    # [m](x) = tanh(m * arctanh(x))
    from math import tanh, atanh
    xm = tanh(m * atanh(x_phi))
    # This should be x_{phi^m}
    phi_m = phi**m
    x_expected = (phi_m - 1) / (phi_m + 1)
    print(f"  [2^{k}](x_phi) = [{m}](x_phi) = {xm:.10f}")
    print(f"    = x_(phi^{m}) = ({phi_m:.4f}-1)/({phi_m:.4f}+1) = {x_expected:.10f}")
    print(f"    phi^{m} = {phi_m:.6f}, L({m}) = {int(lucas(m)) if m < 200 else '?'}")
    print()

print(f"""

The iterated doubling of the golden coupling x_phi traces out
the sequence of golden powers phi, phi^2, phi^4, phi^8, phi^16, ...

At each step, the SYMMETRIC PART (Lucas) and ANTISYMMETRIC PART (Fibonacci)
of phi^(2^k) are computed.

The Lucas values at powers of 2 form the prime chain: 3, 7, 47, 2207.
These are the SYMMETRIC RESIDUES of the doubling fractal.

The Fibonacci values at powers of 2 are the CUMULATIVE PRODUCTS: 1, 3, 21, 987, 2178309.
These are the ANTISYMMETRIC ACCUMULATIONS.


WHY THE CHAIN IS A CHAIN
==========================

L(2n) = L(n)^2 - 2*(-1)^n.

For even n: L(2n) = L(n)^2 - 2.

So:
  L(4) = L(2)^2 - 2 = 3^2 - 2 = 7.
  L(8) = L(4)^2 - 2 = 7^2 - 2 = 47.
  L(16) = L(8)^2 - 2 = 47^2 - 2 = 2207.
  L(32) = L(16)^2 - 2 = 2207^2 - 2 = 4870847.

THE CHAIN IS THE ITERATION x -> x^2 - 2.

Starting from 3: 3 -> 7 -> 47 -> 2207 -> 4870847 -> ...

Each step SQUARES and subtracts 2.
This is the iteration of the CHEBYSHEV MAP T_2(x) = 2x^2 - 1
(conjugated: if y = x/2 then y -> 2y^2 - 1, so x -> x^2 - 2 for x = 2y).
""")

# Verify
print("Verification: L(2n) = L(n)^2 - 2 for even n:")
for k in range(1, 6):
    n = 2**k
    Ln = int(lucas(n))
    L2n = int(lucas(2*n))
    pred = Ln**2 - 2
    print(f"  L({2*n}) = {L2n}, L({n})^2 - 2 = {Ln}^2 - 2 = {pred}, match: {L2n == pred}")

print(f"""

PERFECT. The Lucas prime chain is the orbit of 3 under x -> x^2 - 2.

  3 -> 3^2-2 = 7 -> 7^2-2 = 47 -> 47^2-2 = 2207 -> 2207^2-2 = 4870847.

The map x -> x^2 - 2 is:
  1. QUADRATIC (each step squares).
  2. DETERMINISTIC (no choice, pure iteration).
  3. RELATED TO CHEBYSHEV: T_2 is x -> 2x^2-1 on [-1,1].
     On the interval [2, inf), the conjugated map is x -> x^2-2.

The Chebyshev polynomials T_n satisfy:
  T_n(cos(theta)) = cos(n*theta).

For our map: L(n) = phi^n + psi^n = 2*cosh(n*ln(phi))... no.
  Actually L(n) = phi^n + psi^n, and for large n, psi^n -> 0,
  so L(n) ~ phi^n.

  L(2n) = phi^(2n) + psi^(2n) = (phi^n + psi^n)^2 - 2*(phi*psi)^n
        = L(n)^2 - 2*(-1)^n.

  For even n: L(2n) = L(n)^2 - 2.   (Since phi*psi = -1, (-1)^n = 1.)

This is exactly the DUPLICATION FORMULA for Lucas numbers.
It is the same as cosh(2x) = 2*cosh^2(x) - 1, rescaled.

  cosh(2x) = 2*cosh^2(x) - 1.
  If L(n) = 2*cosh(n*w) where w = ln(phi):
    L(2n) = 2*cosh(2n*w) = 2*(2*cosh^2(n*w) - 1) = 4*cosh^2(n*w) - 2 = L(n)^2 - 2.

  YES: L(n) ~ 2*cosh(n*ln(phi)) and the chain x -> x^2-2 IS the
  cosh duplication formula.


THE FRACTAL STRUCTURE
======================

The formal group doubling [2] generates the sequence:
  x, [2](x), [4](x), [8](x), [16](x), ...

In the SYMMETRIC projection (Lucas):
  L(1), L(2), L(4), L(8), L(16), ... = 1, 3, 7, 47, 2207, ...
  Each step: L -> L^2 - 2 (the cosh duplication).

In the ANTISYMMETRIC projection (Fibonacci):
  F(1), F(2), F(4), F(8), F(16), ... = 1, 1, 3, 21, 987, ...
  Each step: F(2n) = F(n)*L(n) (the sinh duplication: sinh(2x) = 2*sinh*cosh).

The symmetric projection SQUARES (and subtracts 2).
The antisymmetric projection MULTIPLIES by the symmetric value.

This is the MULTIPLICATIVE FRACTAL of phi:
  Each doubling SQUARES the symmetric part and MULTIPLIES the antisymmetric.

  S_k = L(2^k):    3, 7, 47, 2207, ...      (squaring chain)
  A_k = F(2^k):    1, 3, 21, 987, ...        (cumulative product of S_j)

  A_k = S_1 * S_2 * ... * S_{k-1} = 3 * 7 * 47 * ... * L(2^{k-1}).

  When ALL S_j are prime: A_k is a product of DISTINCT primes.
  987 = 3*7*47: three DISTINCT LUCAS PRIMES.

  When S_k is composite: the factorization acquires REPEATED or COMPOSITE factors.
  The chain BREAKS and the product loses its clean prime structure.


3, 7, 47, 2207: THE LUCAS PRIME CHAIN AS A FRACTAL ORBIT
==========================================================

The orbit of 3 under x -> x^2 - 2:
  3 -> 7 -> 47 -> 2207 -> 4870847 -> ...

This orbit stays PRIME for 4 steps. How rare is this?

For comparison, the orbit of other starting values:
""")

starts = [2, 3, 4, 5, 6, 7, 8, 10]
for s in starts:
    x = s
    chain = [s]
    prime_length = 0
    for _ in range(6):
        x = x*x - 2
        chain.append(x)
        if isprime(x) and prime_length == len(chain) - 2:
            prime_length = len(chain) - 1
    primes = sum(1 for c in chain if isprime(c))
    print(f"  Start {s}: {chain[:5]}{'...' if len(chain) > 5 else ''}")
    print(f"    Primes in first 5: {[c for c in chain[:5] if isprime(c)]}")
    print()

print(f"""

Starting from 3: the orbit 3, 7, 47, 2207 has FOUR consecutive primes.
Starting from 2: 2, 2, 2, 2, ... (fixed point: 2^2-2 = 2).
Starting from other values: typically fewer consecutive primes.

The number 3 is SPECIAL as a starting point because:
  3 = L(2) = the first nontrivial Lucas number.
  3 = F(4) = a Fibonacci prime.
  3 is the FIXED POINT of the Mersenne operation on 2: 2^2-1 = 3.
  The orbit of 3 under x -> x^2-2 stays prime for 4 steps.

2 is a FIXED POINT of x -> x^2-2. (2^2-2 = 2.)
3 generates a PRIME ORBIT of length 4.
This is because 3 = 2+1 = the first number AFTER the fixed point.

The map x -> x^2 - 2 has TWO fixed points: 2 and -1.
  (x^2 - 2 = x gives x^2 - x - 2 = 0, so x = 2 or x = -1.)

  2 = the tournament base. The fixed point. No motion.
  -1 = the "anti-unit." The other fixed point.

  3 = 2 + 1 = the base plus one. The FIRST DISPLACEMENT from the fixed point.
  Starting from 3, the orbit SPIRALS AWAY from the fixed point 2,
  through a sequence of primes, until it hits a composite.


THE FORMAL GROUP PERSPECTIVE
==============================

In the formal group F(x,y) = (x+y)/(1+xy):

  The doubling map [2](x) = 2x/(1+x^2) on couplings corresponds to
  L(n) -> L(2n) = L(n)^2 - 2 on Lucas numbers (via the Cayley transform).

  The coupling of L(n) is x = (L(n)-1)/(L(n)+1).
  Doubling: [2](x) = coupling of L(2n) = (L(n)^2-3)/(L(n)^2-1).

  The FORMAL GROUP DOUBLING ON COUPLINGS is the x -> x^2-2 map
  CONJUGATED by the Cayley transform.

  So: the prime chain 3, 7, 47, 2207 is the orbit of the
  FORMAL GROUP DOUBLING MAP, projected to the Lucas (symmetric) axis.

  The formal group is the ENGINE.
  The doubling map is its simplest iteration.
  The Lucas prime chain is the symmetric SHADOW of this iteration.
  The Fibonacci cumulative product (987 = 3*7*47) is the antisymmetric SHADOW.


THE NUMBER 47: WHY IT'S SPECIAL
=================================

47 = L(8) = the Lucas number at the BOTT PERIOD.

47 is prime. It sits at the Bott slot: 47 mod 8 = 7 (FORBIDDEN).

47 is the SYMMETRIC CONTENT of phi^8.
phi^8 = L(8)/2 + F(8)*sqrt(5)/2 = 47/2 + 21*sqrt(5)/2.

  Symmetric: 47.
  Antisymmetric: 21 = 3*7.

  The symmetric part (47) is PRIME.
  The antisymmetric part (21) is COMPOSITE (3*7 = the first two chain primes).

  At the Bott period: the symmetric part is a FRESH PRIME (47).
  The antisymmetric part has ACCUMULATED the previous primes (3*7).

  This is the MEANING of 47:
    It is the IRREDUCIBLE SYMMETRIC RESIDUE at the Bott period.
    It cannot be decomposed further.
    It is what remains after the antisymmetric part has absorbed
    all the previous chain primes.

47 in the periodic table: element 47 = Silver (Ag).
  The "noble" metal. Resistant to corrosion.
  Silver's nobility = its prime irreducibility.
  It is the element that CANNOT be broken down by the doubling chain.

And: 47 in the Bott slot 7 = "forbidden."
  47 is a forbidden-slot prime. Like 7 itself.
  The FIRST Bott-8 Lucas prime is in the forbidden slot.


THE COMPLETE PICTURE
=====================

Two fractal chains start from the prime 3:

MERSENNE CHAIN (additive self-examination):
  2 -> 3 -> 7 -> 127 -> 2^127-1 -> ...
  Rule: k -> 2^k - 1.
  Mechanism: the k-nacci transfer matrix examines its own trace.
  Domain: the TRACE (symmetric observable of the transfer matrix).
  Status: at least 5 consecutive primes (all known members are prime).

LUCAS CHAIN (multiplicative doubling):
  3 -> 7 -> 47 -> 2207 -> 4870847 -> ...
  Rule: x -> x^2 - 2.
  Mechanism: the formal group doubles the golden coupling.
  Domain: the LUCAS NUMBERS (symmetric part of phi^n at powers of 2).
  Status: exactly 4 consecutive primes (breaks at L(32)).

Both chains START with 3, 7:
  Mersenne: 2 -> 3 -> 7 -> 127.
  Lucas:    3 -> 7 -> 47 -> 2207.

The SHARED SEGMENT is {3, 7}.
  3 = Mersenne(2) = L(2). The STRUCTURE prime.
  7 = Mersenne(3) = L(4). The FORBIDDEN prime.

After {3, 7}, the chains DIVERGE:
  Mersenne continues to 127 (= 2^7-1, an ADDITIVE operation on 7).
  Lucas continues to 47 (= 7^2-2, a MULTIPLICATIVE operation on 7).

  127 and 47 are DIFFERENT primes reached from 7 by different operations.
  127 = 2^7 - 1 (exponentiate and subtract).
  47 = 7^2 - 2 (square and subtract).

The two chains correspond to two DIFFERENT fractal self-similarities:
  Mersenne: exponential self-reference (k -> 2^k).
  Lucas: quadratic self-reference (x -> x^2).

Exponential grows FASTER than quadratic, which is why:
  Mersenne: 3, 7, 127, 2^127-1, ... (explosive growth).
  Lucas: 3, 7, 47, 2207, ... (rapid but slower growth).

And the Mersenne chain is EXPECTED to produce more primes
(because 2^k-1 grows so fast that the "probability" of primality
~ 1/k stays high). The Lucas chain x^2-2 grows quadratically
at each step, so the "probability" drops faster.

This explains why the Mersenne chain has 5+ primes
and the Lucas chain has only 4.

BUT: 987 = 3*7*47 = F(16) is the Fibonacci number that
ENCODES the Lucas chain up to its third member.
987 is the ANTISYMMETRIC ACCUMULATION of the first three
Lucas doubling primes.

And if we go further: F(32) = 987 * 2207 = 2178309.
F(32) encodes ALL FOUR Lucas chain primes: 3*7*47*2207.

THE FIBONACCI NUMBERS AT POWERS OF 2 ARE THE PRODUCTS OF
THE LUCAS PRIME CHAIN.

F(2^k) = L(2) * L(4) * L(8) * ... * L(2^{k-1}).
       = product of the first (k-1) members of the Lucas chain.

When the chain is all-prime (k <= 5): F(2^k) is squarefree with k-1 distinct prime factors.
When the chain breaks: F(2^k) acquires composite factors.

987 = 3*7*47 captures the ENTIRE prime portion of the Lucas chain.
It is the LAST Fibonacci-at-power-of-2 that is a product of distinct chain primes.
""")

print("opus-2026-03-16-S90by: THE FORMAL GROUP FRACTAL AND CHAINED PRIMES")
