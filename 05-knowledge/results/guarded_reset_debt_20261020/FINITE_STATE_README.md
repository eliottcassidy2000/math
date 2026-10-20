# Finite-state barriers addendum

Date: 2026-10-20. This is an add-only extension of
[guarded reset-debt composition](RESET_DEBT_COMPOSITION.md).

**PROVED:**

- A finite directed graph admits positive state barriers exactly when every
  cycle has multiplier product at least one and every product-one cycle
  has zero debt on all its edges, under positive multipliers and
  nonnegative debts.
- Rational edge data admit an exact rational feasibility test and a fully
  rational constructive barrier certificate.
- A finite graph of fixed nonempty Collatz valuation words, carrying
  physical-height parameters and positive state barriers, cannot give
  universal continued coverage. An explicit infinite family of all-ones
  prefixes defeats every initial graph state within V edges.

Read [the addendum](FINITE_STATE_BARRIERS_ADDENDUM.md),
[the earlier boundary review](FINITE_STATE_INDEPENDENT_REVIEW.md), and
[the independent publication audit](FINITE_STATE_PUBLICATION_AUDIT.md).
No new source correction was needed after the frozen addendum's explicit
L>=1 boundary repair.

The no-go theorem is restricted to finitely many fixed words and the
physical-height barrier model. It does not cover unbounded-depth
parameterized rules, infinite or unbounded state, different signed/matrix
resources, or source entry by an unbounded preliminary mechanism. It does
not prove Collatz convergence or settle the native signed-row debt problem.

## Reproduce

Both independent implementations need only Python 3's standard library.

    python3 -B verify_finite_state_barriers.py > /tmp/state-producer.json
    python3 -B -O verify_finite_state_barriers.py > /tmp/state-producer-O.json
    cmp /tmp/state-producer.json /tmp/state-producer-O.json
    cmp /tmp/state-producer.json finite-state-receipt.json
    python3 -B independent_finite_state_audit.py > /tmp/state-independent.json
    python3 -B -O independent_finite_state_audit.py > /tmp/state-independent-O.json
    cmp /tmp/state-independent.json /tmp/state-independent-O.json
    cmp /tmp/state-independent.json finite-state-independent.normal.json
    sha256sum -c FINITE_STATE_SHA256SUMS

The producer records 83,452 explicit checks. The independent audit records
146,537, using rational vertex enumeration rather than Fourier–Motzkin
elimination and a separate multiplicative shortest-path implementation.
Finite controls do not replace the arbitrary-graph and all-depth proofs.

[FINITE_STATE_SOURCE_INVENTORY.json](FINITE_STATE_SOURCE_INVENTORY.json)
records the frozen producer hashes, the base-note dependency, and runtime
requirements. This addendum does not replace or alter the base checksum
index or any base file.
