# Session239: efficient exact sequence engines and proof-facing refinements

**Status: original planning checkpoint, 10 November 2026.** The completed
results and independent audits are in [README.md](README.md). No theorem was
promoted by the original checkpoint. Base: `4904ebec2923d9a014296607226a5013acb365bd`.

## User target and typed objects

Generate A000568 and related sequences more efficiently and use the functions
toward proofs. A000568 counts isomorphism classes of tournaments on n vertices,
including the empty tournament at n=0. It is neither the set of Hamiltonian
path counts H(T) nor their maximum. Chan–Pak remains OPEN.

## Inheritance and portfolio

- Anchor: exact Davis/Burnside evaluation. THM-514 (unit-tail core kernel),
  THM-516 (divisor-profile quadratic form), THM-4845 (Gaussian form, existing
  infinite-state recurrence and non-P-recursiveness), and existing integer,
  C/GMP and response-state implementations must be compared, not rediscovered.
- Niche: 2-adic precision in normalized counts. THM-305 already proves the
  odd-order valuation. Seek precision-dependent cycle bounds and actual
  congruences, not a new claim for the inherited parity theorem.
- Wildcard: value-refined fixed-permutation strong-component recurrence, which
  preserves H-value information destroyed by scalar class counting. Recover
  THM-4803, THM-3380, THM-4744, and Session238's action-cost bounds first.

## Concept board and decisive probes

1. Burnside cycle types -> factorial-scaled response states. Preserve exact
   counts; test against Fraction recurrence and independent Davis evaluation.
2. Fixed points -> soluble unit tail. Already known; any improvement must pay
   for the core construction, not just advertise the quadratic outer address.
3. Precision -> bounded number of cycles. Compare exact rational valuations
   with a uniform cutoff and attack repeated parts/factorial denominators.
4. Strong components -> ordered multiplicative H-values. Compare refined
   Burnside rows against directly generated small tournaments.
5. Class counts -> existence claims. Scalar counts forget the requested path
   value; retain the histogram and exhibit that loss instead of inferring a
   Chan–Pak result from fast enumeration alone.

## Verification and publication

Code, exact outputs, hostile controls and reproduction commands will be
committed together. Benchmark old and new algorithms on the same machine;
separate all-n proofs, finite checks and timings. Independent audit follows
implementation. Preserve unrelated checkout line-ending changes and stage
only explicit owned paths. Publish nonforce verified checkpoints to main.

## Sources

- [A000568, offset 0 and Davis formula](https://oeis.org/A000568), checked
  10 November 2026.
- Canonical files cited above are existing repository results, not external
  priority claims. Relevant primary-source statements will be cited wherever
  imported.
