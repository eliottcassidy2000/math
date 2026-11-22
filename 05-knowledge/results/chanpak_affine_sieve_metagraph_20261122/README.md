# Chan–Pak: the affine sieve, exact selectors, and the metagraph

**2026-11-22 · PROVED in the stated scopes. Chan–Pak OPEN.**
This session follows the [branching inverse packet](../chanpak_branching_inverse_20261122/README.md).
It targets the prime and factor-seven lanes, not just a larger finite census.

## Main absolute and arithmetic results

- [A positive-alphabet affine sieve](AFFINE-SIEVE.md). Three fixed legal
  port moves give a uniform congruence spectral gap after comparison with
  two integral shears. A sieve of the actual output words gives infinitely
  many distinct strong path counts with a bounded number of prime factors,
  even in each fixed reduced class modulo a modulus prime to 6. Prescribing
  valuation one at 7 gives an almost-prime version of the factor-seven
  lane. Neither conclusion says the remaining factor is prime.
- [An exact scalar inverse for fixed imbalances](DIGIT-INVERSE.md).
  The large reverse-digit term cancels identically at zero imbalance;
  separated height bands and integer square-root rounding give a complete
  inverse, including that previously excluded chamber. The scope is the
  physical fixed-imbalance slice, not all target counts.
- [A one-apex Monsky character](MONSKY.md). An alternating contact-bit sum
  modulo 3 determines the binary nullity. The existing arithmetic realization
  gives an actual pure 2-Selmer interpretation, not an elliptic rank
  calculation. Two strong hosts of the same order and the same count 9
  have different nullities.

## How the number-theoretic ideas were used

The [Lean and Mathlib joint statement](https://leanprover-community.github.io/blog/posts/riemann-hypothesis/)
certifies Alpöge's classical RH proof, and the
[MIT verification](https://math.mit.edu/~drew/grh/)
also states the verified Dirichlet-GRH scope. These are distinct from
what the [Clay listing](https://www.claymath.org/millennium/riemann-hypothesis/)
records for the prize. The session does not assume an analogous theorem
for arbitrary Hecke, Dedekind, or automorphic functions.

The [proof announcement](https://alpo.ge/remark.pdf) gives a useful
design pattern: get a definite-sign expression and iterate a gain
that can pay for an arithmetic remainder. In this session the exact
positive expression is finite variance of two-step port moves. The
unconditional comparison uses Bourgain–Varjú, not an RH estimate.
RH itself supplies a clean parity-barrier example: two integer
weights have essentially the same divisibility data even though one
is zero on all primes. This explains the distinction between sieving
for bounded almost-primes and finding the precise missing prime targets.

The [BSD](https://www.claymath.org/millennium/birch-and-swinnerton-dyer-conjecture/)
lesson is to distinguish a finite local invariant, a height count, and
an actual rank or physical family. The Monsky theorem is an arithmetic
transport with an explicit matrix; neither the analytic nor algebraic
rank of the associated elliptic curve is asserted. The
[Hodge](https://www.claymath.org/millennium/hodge-conjecture/)
lesson is to keep the actual cycles and their allowed local fillings,
rather than equating a scalar invariant with a lift. The selected
tournament faces are combinatorial objects, not algebraic cycles on
a projective variety.

Exact local controls for the new formulas are included. The proofs,
not the finite checks, supply the unbounded assertions.
