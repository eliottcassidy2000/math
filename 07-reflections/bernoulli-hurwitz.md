# B_6 = 1/42: Why Bernoulli Knows About Hurwitz

**Session:** kind-pasteur-2026-03-16-S116n

## The Fact

The sixth Bernoulli number is B_6 = 1/42.

This is the ONLY Bernoulli number whose denominator equals the Hurwitz constant 42.

## Why

The von Staudt-Clausen theorem (1840): denom(B_{2k}) = prod{p prime : (p-1) | 2k}.

For 2k = 6:
- Divisors of 6: {1, 2, 3, 6}
- Adding 1: {2, 3, 4, 7}
- Primes among these: {2, 3, 7}
- Product: 42

The Hurwitz primes are *exactly* the primes p with (p-1)|6.

## The Deep Reason

The rational quaternion algebra H_Q ramifies at primes where it fails to split. This happens at p=2, 3, 7 — the primes where (p-1) divides 6 = dim(H_Q) + 2. The discriminant of the maximal order (Hurwitz quaternions) is therefore 2*3*7 = 42.

Von Staudt-Clausen and Hurwitz quaternions detect the SAME arithmetic obstruction: the primes p where the (p-1)-cyclotomic structure interacts with the index 6.

## Consequences

- zeta(6) = pi^6/945, and 7 | 945 because B_6 = 1/42
- zeta(-5) = -1/252 = -1/(4*63) = -1/(2^2 * 3^2 * 7) — the Hurwitz primes in negative zeta values
- The p-adic zeta function zeta_7(s) has a simple pole at s ≡ 0 mod 6, from v_7(B_6) = -1
- B_{42} has p=43 in its denominator (because 43-1 = 42), connecting to the Sylvester sequence

## The Uniqueness

denom(B_{2k}) = 42 holds ONLY at k=3 (2k=6). At k=9 (2k=18), p=19 enters because 18|18. At k=21 (2k=42), p=43 enters because 42|42. The Bernoulli number B_6 is the unique point where the arithmetic of 42 is visible in the simplest possible form.

## Connection to Everything

- B_2 = 1/6: denom = 2*3 = the smallest Hurwitz primes
- B_4 = -1/30: denom = 2*3*5 = the Platonic system
- B_6 = 1/42: denom = 2*3*7 = the Hurwitz system
- B_12 = -691/2730: denom = 2*3*5*7*13 = Platonic AND Hurwitz together

The Bernoulli denominators are a SIEVE: each 2k picks out exactly the primes p where p-1 divides 2k. The sequence of denominators IS the von Staudt sieve, and 42 is the third value in this sieve after 6 and 30.

The weight-12 modular form Delta(tau) has weight phi(42) = 12. The discriminant modular form IS the Cayley chromatic number of the Hurwitz constant. This closes the circle.
