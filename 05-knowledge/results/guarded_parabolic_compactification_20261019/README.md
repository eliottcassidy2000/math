# Guarded parabolic compactification

[The proof note](GUARDED_PARABOLIC_COMPACTIFICATION.md) establishes a positive Collatz result with an explicit additional guard: a contracting parity word followed by itself cannot have positive endpoint overshoot. Applied to the K,V,F,H first-crossing bank, this removes the episode cap under the self-matching hypothesis and retains the smaller-companion-cycle alternative.

It also identifies the exact parabolic Möbius map obtained from near-neutral bank words, proves that integer repetition guards force sources out of its nonlinear crossover scale, and gives a necessary-and-sufficient compact-boundary criterion for fixed-offset uniform original-height repair at large sources. Actual source/displacement limits remain distinct from formal infinite bank paths and the K-infinity rational anchor.

## Verification

- [Main exact checker](verify_guarded_compactification.py) and [receipt](verification-receipt.json): 137,320 checks
- [Independent proof review](independent-review/INDEPENDENT_REVIEW.md): PASS, within scope; 1,124,584 independent exact checks
- [Independent checker](independent-review/independent_controls.py), [normal receipt](independent-review/independent-normal.json), and [optimized receipt](independent-review/independent-optimized.json)
- [File manifest and proof revision hashes](MANIFEST.json), [checksums](SHA256SUMS)

Run the main checker normally and with Python's `-O` option; both outputs reproduce its receipt. The independent checker is self-contained and can also be run normally or with `-O`.

Status: scoped proofs plus finite exact controls. No unconditional pointwise bank entry, universal rational coefficient-stopping theorem, or Collatz convergence claim.
