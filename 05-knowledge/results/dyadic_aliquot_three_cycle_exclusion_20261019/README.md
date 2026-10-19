# All-exponent dyadic aliquot three-cycle exclusion

The [proof](DYADIC_THREE_CYCLE_EXCLUSION.md) resolves A3 and D3 from the [earlier dyadic frontier](../aliquot_dyadic_frontier_20261014e.md): no primitive aliquot three-cycle can have every node equal to 2^a p, with a>=1 and p an odd prime.

The distinct-core theorem excludes even a formal integral solution, using both cyclic orientations, a strict size interval, a congruence, and 2-adic valuations. It has no exponent or source-height cutoff. The earlier repeated-core classification supplies the remaining cases. It does not exclude a general aliquot three-cycle with more complicated odd parts.

## Reproduce

Python 3 and SymPy are required.

- Run `python verify_dyadic_three_cycle.py` here and compare output with `verification.json`.
- Repeat under `python -O`; output must be identical.
- Run `python review/independent_audit.py` and compare with `review/independent_audit.json`.
- Repeat the independent checker under `python -O` as well.

The [independent review](review/INDEPENDENT_REVIEW.md) verifies the unrestricted algebraic proof, its inherited repeated-core cases, and independent exact controls. A nondyadic formal example shows why the powers-of-two hypothesis matters. Finite checks supplement the proof rather than establish its all-exponent quantifier. Historical priority has not been established or claimed.
