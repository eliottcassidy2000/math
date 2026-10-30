# Session 110 — symmetry, prime composition, arithmetic information

Completed research packet, 30 October 2026.

## Main results

1. **Logarithm-free prime tournament construction.** In the simply ramified quadratic family, with N=p^(2k), compatible balanced lifts exist with skew adjacency norm at most 200√(pN). Every such balanced lift is prime for k≥2. Every bank member has norm at least √(pN), so the existence bound has the correct order in both variables. This is an existence theorem using discrepancy theory, not an implemented optimizer. See [proof](LOGFREE_PRIME_TOURNAMENTS.md).

2. **Prime composition with the actual action retained.** For a fixed odd transitive action G/K, the prime count is the incidence inverse of the double-coset orientation count:
A(G,K)=2^{(|K\\G/K|−1)/2},
P(G,K)=A(G,K)−Σ_{K<H<G}P(G,H)A(H,K).
This covers action presentations, not a disjoint enumeration of graph classes. Direct enumeration gives992 prime Cayley states for C21 and980 for C7⋊C3, out of1024 in either case. See [checkpoint2](CHECKPOINT_02.md).

3. **A precise small ramified atlas, with hostile higher-depth controls.** The order125 bank has408 prime graph classes and one composite graph class. At order625 the corrected normalizer has order5000, giving83,232 prime full-affine classes. Full graph classification at625 remains open. At3125, an additional order-five scalar fixes17,045,913,600 prime states, defeating extrapolation of the smaller stabilizer pattern. See [checkpoint2](CHECKPOINT_02.md) and the portable normalizer, graph and depth-five checks.

4. **Tournament enumeration via invariant subspaces.** Palette symmetry can be replaced by row-space data, with Gaussian-binomial formulas in semisimple cases and primary-module formulas in modular characteristic. Strong-component inversion must occur before averaging. The binary implementation reproduces A000568 through order12 and passes280 termwise comparisons with explicit palette groups. This extends the enumeration framework; it is not a claim that unrestricted prime graph classes have been classified. See [module enumeration](MODULE_ENUMERATION.md).

5. **New exact arithmetic information boundaries.** For a fixed prime power p^t, the pair s(n),s(p^t n) uniquely recovers n. Composite multipliers behave differently: complete two-source fibres occur for10 and15. The square-class controls exhibit a projected Fano pattern with16 apparent square selections but only one genuine selection. Full hidden rank, not a selected prime projection, determines the missing information. See [checkpoint1](CHECKPOINT_01.md) and [arithmetic checks](portable_checks/arithmetic_probes.js).

6. **Quinary path and augmentation recurrences.** The orbit-based Hamilton-path carrier for C5[C5] has65,536 states and a24,297-state reachable weighted level-preserving quotient. Its exact path count reproduces the earlier October21 result. The corresponding regular-parent family has6,559 unmarked order26 extensions, including1,554 rigid ones. These extension counts are shared with concurrent work and are not counts of all order26 tournaments. See [quinary carriers](QUINARY_CARRIERS.md).

7. **Two further prime-existence results.** A finite odd abelian group has a prime Cayley tournament unless it is C3². In the nonmaximal dual-number ring, three essential cuts yield an explicit inclusion–exclusion prime formula. At p=5 this gives4,244,832,000 prime labelled states out of2^32. Priority of the general abelian existence statement has not been established. See [additional proofs](ADDITIONAL_PROOFS.md).

8. **A source-preserving limit for Collatz certificates.** Exact signed residual identities retain source weights and negative credit. Finite modifications leave the reciprocal tail coefficient equal to one; an infinite reset construction can have a smaller coefficient without contradicting convergence in ℓ1. Supersolution endpoint data impose a Bellman no-gain boundary on finite signed optimization. See [additional proofs](ADDITIONAL_PROOFS.md). No Collatz convergence theorem is claimed.

## What unifies the results?

A useful organizing principle is to retain the object that composition actually acts on:

- a subgroup interval, rather than just a group order;
- a joint module/quotient type, rather than a marginal dimension;
- an embedded action, rather than only an abstract automorphism group;
- a full arithmetic syndrome, rather than its visible prime coordinates;
- original source weights and native clocks, rather than an averaged surrogate;
- a specified family of observations, when asserting minimal state counts.

This is a proposed research framework, not one new universal theorem. The precise theorems above show where it succeeds; the higher-depth stabilizer and false-Fano controls show why coarser versions fail.

A useful distinction is between counting, recognition and continued composition. Acceptance-only linear rank25 in the quinary path problem does not yield a free25-state generator: obtaining the completion coordinates already requires the path computation, and that representation does not retain all level observations.

## Verification and reproducibility

The [machine-readable receipt](portable_checks/RESULTS.json) records exact checks and their scope.

Published portable JavaScript was fetched at immutable commits and replayed in the parent execution isolate:
- first four checks: direct order21 graph census,625/125 normalizers,125 graph census, depth-five hostile control;
- seven further checks: binary enumeration, ramified cyclotomic period identities, dual-number ideals, corrected quinary paths, arithmetic inverses, regular augmentation, and abelian existence controls.

The seven later sources were fetched from f31384340adb9f957ca4299a8cab4f7755807072. The receipt includes their blob hashes and complete outputs. Scripts use no external packages and can be run separately with a JavaScript runtime supporting BigInt.

**Failed consolidated replay:** the parent Python batch terminated with a240-second timeout on its first script during execution-service disruption. It is not a passing14-run verification. Earlier individual producer and independent-auditor checks are distinct evidence. No timing speedup is claimed.

A validation-only dangling-else error in the initial portable quinary checker was repaired in a039df2bacadc46f3aaece7c5188d0474681c3f5. The corrected run explicitly covers2,304 ternary raw-transition states and1,000 quinary controls. The Hamilton count and partition counts were unchanged. The central mistakes ledger records this and the normalizer corrections.

Finite computations support the stated finite ranges. The all-order conclusions rest on their accompanying proofs and cited classical results, not extrapolation from these checks.

## Concurrent work and attribution

The packet incorporates, rather than reclaims, earlier and concurrent results:
- the exact25-vertex Hamilton count from lattice_bridges_20261021;
- regular-parent symmetry and the shared quinary extension inventory;
- the parallel fixed-p logarithm-free balanced-lift construction;
- source-cutoff, reset-tail and Bellman analyses;
- classical discrepancy, invariant-subspace, Hall/Birkhoff and modular-decomposition inputs.

Detailed links and distinctions appear in the individual notes. Numerical coincidences with Sha, the Monster or Collatz are not promoted to isomorphisms.

## Focused next steps

- Implement the constructive discrepancy method and certify an explicit balanced prime tower.
- Close the gap between affine and full graph isomorphism at625.
- Turn joint subgroup/quotient types into reusable enumeration states without losing action embeddings.
- Test the enhanced3125 stabilizer strata before making an all-depth rigidity conjecture.
- Seek independently certified endpoint information for the Collatz residual method, rather than optimizing a supersolution tautology.

The current packet closes this overnight session. These are research directions, not claims of completed classifications or scheduled follow-up work.
