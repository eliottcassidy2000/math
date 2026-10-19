# Higher-order bottleneck clocks

Read [the proof](HIGHER_ORDER_BOTTLENECKS.md) for uniform discrete escape-time estimates and their base-independent interpretation.

The pure even-power family has count K_r epsilon^(-(r-1)/r)+O(1). A quartic family sharing the same leading term has a parameter-dependent correction of order epsilon^(-1/4). These estimates make the lower-digit information missing from a leading scaling approximation explicit. They do not establish a prime or Collatz statement.

The proof has passed independent mathematical review. Higher-order continuous-time scaling and quadratic escape asymptotics have substantial cited prior literature; no historical priority is asserted.

## Reproduce

Only Python's standard library is required.

- Run `python verify_higher_bottlenecks.py` here; compare with `verification.json`.
- Run `python -O verify_higher_bottlenecks.py`; the output must be identical.
- Run `python independent-review/independent_checks.py`; compare with `independent-review/independent-verification.json`.
- Run the independent checker with `python -O` as well.

The author's escape counts use directed repeated Decimal products; the independent checker uses exact integer-grid enclosures. Transcendental comparisons in the first checker are explicitly numerical controls. The written proof, not a finite sample, establishes parameter-uniform bounds.

See [the independent review](independent-review/REVIEW.md) for checked constants, endpoint estimates, attribution, and four additional exact controls. Source hashes in its receipt refer to the packaged proof, author verifier, and author receipt. The only adaptation of the independent checker for this package is its parent-directory lookup.
