# HYP-2448 - Irreducible-prime certificate-state tournament

**Status:** OPEN synthesis; finite atlas confirms the carrier shape.
**Source:** codex-2026-06-12.
**External sources:** Jitender Singh, arXiv:2411.18366,
https://arxiv.org/abs/2411.18366; Shahriar Iravanian, arXiv:2410.15880,
https://arxiv.org/abs/2410.15880.
**Companions:** HYP-2227, HYP-2425, HYP-2430, HYP-2444, HYP-2445,
HYP-2447,
OPEN-Q-070.
**Computation:** `04-computation/irreducible_prime_carrier_tournament_codex.py`;
stored output `05-knowledge/results/irreducible_prime_carrier_tournament_codex.out`.

## Statement

The repo's recurring irreducibility theme can be organized as a bidirectional
atom carrier:

```text
atoms in Z[x]  --evaluation-->  atoms in Z
irreducible primitive f         prime values f(n)

atoms in Z     --certificates--> atoms/factor bounds in Z[x]
prime or low-Omega values        irreducibility or factor-count bounds
```

Bunyakovsky/Buniakowski is the forward direction: if `f in Z[x]` is
irreducible and the fixed divisor of `f(N)` is `1`, then `f` should take
infinitely many prime values.

Singh/Murty and Cohn give the reverse direction: if a sufficiently large value
`f(m)` has few prime factors, or if the base-`b` digit value `N=f(b)` has few
prime factors, then the number of irreducible factors of `f` is bounded.  In
particular, prime values can certify irreducibility.

Iravanian supplies the algorithmic middle: factor over `R`, then recombine
linear/quadratic real factors by subset-sum/integer-trace constraints.  In
that lens, irreducibility is absence of a nontrivial subset that survives all
integer coefficient tests.

Thus the clean carrier is not "primes are tournaments" or "polynomials are
tournaments".  The clean carrier is the finite/infinite tournament of
certificates:

```text
vertices = proof channels or candidate factor subsets;
edge A -> B if A separates atoms/factors with less lost support than B.
```

## Exact Finite Evidence

The stored atlas confirms the basic split on small rows.

1. `x^2+x+2` is irreducible but has fixed divisor `2`, so it is excluded by
   Bunyakovsky admissibility.  It has only the trivial prime value `2` in the
   stored window.
2. `x^2-x+41`, `x^2+x+1`, `x^2+1`, and `x^4-10x^2+1` are primitive
   irreducible rows with fixed divisor `1`.  Each receives a Singh/Murty
   one-prime certificate in the tested range.
3. Cohn's digit mechanism is visible exactly: `101` in base `10` gives
   `x^2+1` with `Omega(101)=1`; `9841` in base `3` gives
   `1+x+...+x^8`, with factor degrees `[2,6]` and
   `Omega(9841)=2`.
4. Singh's reducible degree-8 example `x^8+6x^7+5` has factor degrees
   `[1,7]`, and the stored search finds a tight plain value certificate at
   `m=12`.
5. Iravanian-style first-trace recombination is necessary but not sufficient:
   the irreducible quartic `x^4-10x^2+1` has two false integer-trace subset
   candidates.  The reducible product
   `(x^2+x+1)(x^2+3x+1)` has exactly the two nontrivial trace candidates
   expected from its two factors.

As an addendum to HYP-2447's four-face prism, this is the first useful
computational form of the user's two-way sentence:
irreducible polynomials can produce primes, and primes can produce
irreducible-polynomial certificates.

## Tournament Analysis

The atlas uses proof carriers as vertices:

```text
bunyakovsky_forward_atom
singh_murty_value_factor_bound
cohn_digit_prime_certificate
iravanian_real_trace_recombination
finite_field_hensel_recombination
newton_nonarch_bivariate_gate
local_fixed_divisor_sieve
repo_support_gate_LRC72
```

Pairwise observable:

```text
majority(forward_primes, reverse_irreducibility, recombination,
         support_retention, repo_bridge)
```

Tie path: lower risk, then declaration/name.

Fingerprint of the stored run:

```text
score_hist = {0:1, 2:2, 3:1, 4:2, 6:1, 7:1}
directed_3cycles = 3
SCC sizes = [5,1,1,1]
Hamiltonian paths = 9
edge_flips_vs_reverse_irreducibility_only = 8
leader = singh_murty_value_factor_bound
runner-up = repo_support_gate_LRC72
```

The nontransitivity is useful.  It says no single scalar ranking controls the
field: forward prime production, reverse certificates, recombination
algorithmics, and retained support channels trade off.

## Transfer to LRC14 and the 72-Code Gate

This hypothesis sharpens the support-gate program from HYP-2445.

For LRC14, "blocked denominator" is analogous to "composite value" only after
the fixed-divisor/local-obstruction channel is retained.  The next resource
lemma should not merely count blocked `q`; it should track an irreducible
resource vector over shell-27 classes, the 13-clock, divisor fibers, and
Bprime/owner-deletion exceptions.

For the extremal Type II `[72,36,16]` problem, the scalar weight enumerator is
the analogue of a value of a polynomial: it is necessary but lossy.  The
irreducible object is the support/matroid/design realization.  HYP-2430's
Tutte support gate is the code-side version of asking whether the apparent
atom in the scalar quotient genuinely lifts to an atom in the richer category.

The slogan becomes:

```text
prime value : irreducible polynomial
scalar enumerator : irreducible support design
blocked q scalar : irreducible LRC resource vector
```

In all three cases, the proof must retain the channel that the scalar quotient
forgets.

## Infinite-Tournament Program

Define certificate states for primitive integer polynomials:

```text
C(f;X) = (fixed divisor, local residue obstructions up to X,
          least Singh/Cohn certificate depth up to X,
          trace-subset survivor profile up to X,
          Newton/non-Archimedean support data when present).
```

Orient `C(f;X) -> C(g;X)` if `f` has smaller normalized unresolved
factorization ambiguity than `g`; orient ties by the richer retained side
channel.  Then study edge flips as `X` grows.

Bunyakovsky becomes the assertion that admissible irreducible vertices keep
receiving forward prime-hit witnesses infinitely often.  Singh/Cohn/Iravanian
give reverse and algorithmic edges in the same infinite tournament.

## Assumption Challenge

Alternate vertex sets considered: polynomials, integer values, primes,
residue classes, roots, real/complex conjugate blocks, factor subsets, digit
positions, Newton polygon edges, bivariate specialization channels, LRC
runners, gaps, fixed circle sections, section boundaries, wall-crossing
events, codewords, matroid circuits, design blocks, and proof obligations.

Chosen vertex set: proof carriers/certificate states.

Preserved: the irreducible-prime two-way relationship, fixed-divisor/local
obstruction data, reverse certificate strength, recombination ambiguity,
support retention, and transfer potential to LRC14 and `[72,36,16]`.

Destroyed: actual asymptotic prime production, detailed root geometry beyond
the first recombination trace, individual code support realizability, and the
full arithmetic of each polynomial.

Challenged assumption: a Tournament Analysis vertex need not be a runner, arc,
prime, or polynomial.  For this session, the predicate-preserving quotient is
the certificate channel.

## Next Moves

1. Replace the floating recombination scout by exact algebraic trace lattices
   and count all false trace survivors for standard hard irreducible families.
2. Build `C(f;X)` for a larger polynomial family and measure edge flips as
   `X` grows; look for stable SCCs that correspond to known irreducibility
   criteria.
3. Translate the same carrier to LRC14: certificate states become Q27 resource
   vectors and Bprime/owner-deletion exceptions.
4. Translate the same carrier to `[72,36,16]`: certificate states become
   support/matroid/design realization moves rather than scalar weight
   enumerators.
