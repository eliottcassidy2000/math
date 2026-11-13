# Independent review of joint valuation refinement

**Status: VERIFIED; no unresolved defect found.** This review covers
[JOINT_REFINEMENT.md](arithmetic/JOINT_REFINEMENT.md) and its stated incoming
dependency, [Session292, arithmetic section 3](../session292_triple_bridges_20261113/arithmetic/PROOF.md),
commit `8e0debd6b2`.

The incoming theorem applies to all three displayed sequences: their
coefficients and constants are positive, the base is 2, and each constant
is odd. Consequently every prime divisor is odd, so all multiplicative
orders modulo the required prime powers exist. No assumption of prime
values or of coprimality between the coefficient and constant is needed.

For a finite list of already forced primes, choosing a witness exponent
and retaining its residue modulo each `ord_(p^(e+1))(2)` preserves the
**exact** valuation `e`. Together with the old exponent modulus these
constraints are compatible: that same witness solves all of them. The
refined progression is infinite. Since the positive exponential tends to
infinity while the product of the fixed prime powers is constant, a larger
term must have a new prime divisor. Pinning that divisor and repeating
proves the induction. This explains why merely preserving divisibility,
without exact valuation control, would have been insufficient.

Applying the argument to A, then B, then C is valid. Every later refinement
retains every earlier congruence, and all witnesses are chosen from the
current progression. Thus this is a simultaneous construction using one
exponent, not three independently chosen CRT phases. It works in each
accepted joint-wheel residue and for every prescribed nonnegative triple
of factor-count lower bounds. Acceptance excludes every prime at most 97.

The three gcd identities in the note are exact because A, B, C are odd
and B=3A+2, C=3A+4. Therefore the forced prime lists are disjoint, even
before imposing the sieve. The claim is a lower bound on distinct factors;
it does not claim an exact factorization, a uniform bound on the required
exponent, or a new prime-producing construction.

The square-plus-one exclusion also checks: A and B are 3 modulo 4 for
`k>=2`; C-1 is 8 modulo 16 for `k>=4`. The five separately listed head
values 220, 664, 666, 1260 and 2448 are nonsquares. Therefore the exclusion
is for every entry of the formal triple and every positive exponent,
independently of primality.

For the final A135244 consequence, choosing all three factor lower bounds
at least 2 immediately gives simultaneous strict inequalities. Alternatively,
forcing one fixed divisor in each coordinate and discarding finitely many
initial terms makes all three divisors proper. A coordinate with requested
factor bound zero, by itself, supplies no compositeness conclusion. The
note's reference to forced-composite subprogressions is read with this
necessary guard, consistent with its earlier explicit progression.

Cheap independent controls checked pairwise coprimality and square exclusion
for every `k=1,...,256` using exact integers. Complete trial division at
accepted exponents 6 and 30 additionally gave

```text
k=6:  (6359, 19079, 19081), each prime.
k=30: A = 3803 * 27951733,
      B = 318901321799 (prime),
      C = 163 * 1956449827.
```

All displayed factors are prime and exceed 97, with pairwise disjoint
lists. These are finite corroborations, not the basis for the arbitrary
factor-count induction. The proof's new simultaneous conclusion and the
incoming single-sequence theorem remain separately attributed.
