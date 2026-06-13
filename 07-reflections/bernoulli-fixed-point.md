# The Bernoulli Fixed Point: 6 -> 42 -> 1806 -> 1806

**Session:** kind-pasteur-2026-03-16-S116n32

## The Chain

Define the Bernoulli chain: start with 2k, compute denom(B_{2k}) via von Staudt-Clausen, use that denominator as the next 2k.

- B_6: primes with (p-1)|6 are {2,3,7}. Product = 42.
- B_42: primes with (p-1)|42 are {2,3,7,43}. Product = 1806.
- B_1806: primes with (p-1)|1806 are {2,3,7,43}. Product = 1806.

**1806 is a FIXED POINT.** The chain stabilizes in exactly 3 steps.

## Why It Stabilizes

1806 = 2 * 3 * 7 * 43. Its divisors are {1, 2, 3, 6, 7, 14, 21, 42, 43, 86, 129, 258, 301, 602, 903, 1806}.

Adding 1 to each divisor: {2, 3, 4, 7, 8, 15, 22, 43, 44, 87, 130, 259, 302, 603, 904, 1807}.

The primes among these: {2, 3, 7, 43}. Product = 1806.

The crucial non-entries:
- 13 would need (13-1)=12 | 1806. But 1806 mod 12 = 6. EXCLUDED.
- 139 would need 138 | 1806. But 1806 mod 138 = 12. EXCLUDED.
- 1807 = 13*139 is composite. Not a prime at all.

So 1806 is **self-selecting**: the only primes p with (p-1)|1806 are exactly the prime factors of 1806 itself. No new prime can enter the set.

## The Sylvester Connection

The Sylvester sequence 2, 3, 7, 43, 1807, ... has partial products 2, 6, 42, 1806, ...

The Bernoulli chain agrees with Sylvester's partial products through P_4 = 1806, but Sylvester's next term a_5 = 1807 = 13*139 is composite. When this happens, the Bernoulli chain cannot absorb the composite number as a new prime factor, so it stabilizes.

**The Bernoulli chain is Sylvester's sequence with a primality filter.** It follows Sylvester as long as the next term is prime, and freezes when a composite term appears.

## Self-Selecting Sets of Primes

A set S of primes is **von Staudt self-selecting** if: the primes p with (p-1) | prod(S) are exactly S itself.

{2, 3} is self-selecting: (p-1)|6, primes are {2,3}. Fixed point at 6.
{2, 3, 7} is NOT self-selecting: (p-1)|42 gives {2,3,7,43}.
{2, 3, 7, 43} IS self-selecting: (p-1)|1806 gives {2,3,7,43}. Fixed point at 1806.

The chain 6 -> 42 -> 1806 passes through a non-fixed-point (42) to reach a fixed point (1806).

## Why This Matters

42 = denom(B_6) is the TRANSITIONAL value — neither a starting point (6) nor a fixed point (1806), but the bridge between them. It is the unique intermediate step in the Bernoulli chain from the simplest nontrivial fixed point {2,3} to the next fixed point {2,3,7,43}.

The fact that 42 = 2*3*7 attracts exactly one new prime (43) to form a self-selecting set is arithmetically rare. Most numbers n attract multiple new primes when you compute the VS primes for (p-1)|n, creating an expanding cascade. 42 attracts exactly one, and that one closes the system.
