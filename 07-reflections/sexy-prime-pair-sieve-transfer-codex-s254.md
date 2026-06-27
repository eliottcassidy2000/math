# Sexy prime pair sieve transfer (codex-2026-06-27-S254)

The repo relates to the sexy prime conjecture most cleanly through the old
pair-lens work, not through a direct LRC theorem. Sexy primes are the fixed
gap-6 row `(p,p+6)=(m-3,m+3)`, so they are the `h=3` ray in the midpoint-gap
surface from HYP-1966 and HYP-1965.

The local story is already sharp. For every odd prime `q`, the pair is killed
when `m=+h` or `m=-h mod q`; for `h=3` those two residues collapse modulo `3`.
That gives the Hardy-Littlewood chord factor `2` relative to twins. The S254
scout through `10^6` sees `16386` gap-6 pairs versus `8169` twin pairs, ratio
`2.006`, matching the predicted local shape.

That is evidence for the correct coordinate system, not a proof. The missing
object is still a lower-bound fixed-gap sieve that detects two primes and
breaks the parity barrier. Counts, local admissibility, and singular-series
ratios explain why the conjecture is expected; they do not prove infinitely
many pairs.

The useful transfer from LRC14 is the finite-address discipline. Before a proof
forgets a coordinate, it should say which predicate is preserved and what debt
is paid. For sexy primes the preserved predicate is infinitely many `m` with
`m-3` and `m+3` prime. The destroyed data are endpoint identities, local
residue channels, prime-power sidecars, exceptional moduli, weight choices,
and parity/distribution debt.

The incoming level-7 and gK8 work sharpens the analogy but not the theorem.
THM-573 is a seven-lift pigeonhole sieve for LRC covering rows, reducing the
hard residual to at most six multiples of `7`. HYP-3084 now says the gK8
covariance route is not a literal Clebsch design Gram; the surviving certificate
is a reflection-symmetric pairwise co-emptiness matrix with a dominant Perron
mode and a `3x3` symmetric-block target. That is the same style of local gate,
residual ledger, and low-order certificate a sexy-prime proof would need, but
it is not a prime-pair distribution theorem.

The LRC singular series is not the Hardy-Littlewood Euler product; one is a
geometric/archimedean proof object in the LRC stack, the other is a local
congruence product for prime pairs.

The right next move is therefore not to claim the repo proves sexy primes. It
is to build a `sexy_prime_pair_ledger` that makes the unpaid obligations visible
after the easy local work is complete: lower-bound sieve, parity breaking,
distribution in arithmetic progressions, prime-power sidecars, and named
exceptional exits.

Tournament Analysis assumption challenge: the vertices should not default to
prime endpoints or raw pair counts. Better vertices are proof obligations,
fixed-gap rays, residue channels, and weight sidecars. This preserves the exact
gap-6 infinitude predicate while recording exactly what each quotient destroys.
