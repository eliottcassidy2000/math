# The Bernoulli Fixed Point: 6 -> 42 -> 1806 -> 1806

**Historical note, corrected 2026-10-29.** The displayed chain and fixed
point1806 are valid. The original claim that {2,3} was already fixed at6
was false: 7-1 divides6. The suggested universal Sylvester/primality-filter
rule was also false; n=8 has composite n+1=9 but maps to30. See the
[exact repair and connection atlas](../05-knowledge/results/twelve_twenty_two_trace_20261029/ARCHAEOLOGY.md)
and its `check_precision.py` receipt. This reflection is provenance, not a
general classification of Bernoulli-denominator dynamics.

**Session:** kind-pasteur-2026-03-16-S116n32

## The Chain

Define the Bernoulli chain: start with 2k, compute denom(B_{2k}) via von Staudt-Clausen, use that denominator as the next 2k.

- B_6: primes with (p-1)|6 are {2,3,7}. Product = 42.
- B_42: primes with (p-1)|42 are {2,3,7,43}. Product = 1806.
- B_1806: primes with (p-1)|1806 are {2,3,7,43}. Product = 1806.

**1806 is a FIXED POINT.** Starting at6, two updates reach1806; the third
displayed update confirms that it remains there.

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

This particular Bernoulli chain agrees with the displayed Sylvester partial
products through1806. Its fixed-point proof uses the COMPLETE divisor list
above: the prime values d+1 are exactly2,3,7,43. The compositeness of1807
alone would not establish the fixed point; other divisors can introduce
new primes.

The historical sentence identifying every Bernoulli chain with a
primality-filtered Sylvester sequence is retracted. At n=8, n+1=9 is
composite, but divisors1,2,4 contribute primes2,3,5 and the next value is30.

## Self-Selecting Sets of Primes

A set S of primes is **von Staudt self-selecting** if: the primes p with (p-1) | prod(S) are exactly S itself.

{2, 3} is NOT self-selecting: (p-1)|6 gives primes {2,3,7}, so6 maps to42.
{2, 3, 7} is NOT self-selecting: (p-1)|42 gives {2,3,7,43}.
{2, 3, 7, 43} IS self-selecting: (p-1)|1806 gives {2,3,7,43}. Fixed point at 1806.

The chain 6 -> 42 -> 1806 passes through a non-fixed-point (42) to reach a fixed point (1806).

## Why This Matters

42 = denom(B_6) is an intermediate value in the exact chain starting at6
and ending at the fixed point1806. The starting set {2,3} is not fixed.

42 adds the prime43, and the complete divisor test at1806 then closes
this particular chain. No density estimate or assertion that most inputs
generate expanding cascades is established by this example.
