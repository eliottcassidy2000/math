# Matched factors and the repetition-poor first-crossing boundary

These independently audited extensions use the [primary expansion-first theorem](../INTERLEAVED_REPETITION_BARRIER.md) and preserve its original proof and review.

- [Matched-block extension](MATCHED_BLOCK_EXTENSION.md): sufficiently long identical chronological factors force either strict endpoint descent or an exact internal return. The arbitrary-bit-start and factor-complexity versions retain the same cycle alternative.
- [Repetition-poor greedy family](REPETITION_POOR_GREEDY_BOUNDARY.md): a fixed mechanical construction has lawful cylinder lifts but no macroscopic pure-expansion power at the stated thresholds. The cutoff is m>900 for its displayed parse and m>915 for arbitrary exact K^4F/V reparsing. Its literal factor-complexity bound still gives the separate matched-factor obstruction.

An internal periodic orbit can have a nonzero preperiod. The exact-cycle alternative is not eliminated; the notes neither construct non-descending sources nor establish unbounded membership in a fixed positive-displacement source set. No general Collatz convergence conclusion follows.

## Reproduce

Only Python 3 and the standard library are required. From this directory:

- `python3 -B verify_boundary_extensions.py`, then `python3 -B -O verify_boundary_extensions.py`: compare with the author normal/optimized receipts (49,321 checks).
- `python3 -B review/independent_extension_controls.py`, then `python3 -B -O review/independent_extension_controls.py`: compare with the independent receipts (178,380 checks).

The [independent audit](review/MATCHED_BLOCK_EXTENSION_REVIEW.md) covers the proofs and corrected parse thresholds. All executable checks stay active under optimization. Finite controls supplement the universal arguments. `SOURCE_PROVENANCE.json` records frozen source hashes and the presentation-only cleanup; `review/SHA256SUMS` checks the published review artifacts.
