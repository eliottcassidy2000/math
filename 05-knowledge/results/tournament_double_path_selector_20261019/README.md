# The forest primitive floor and the eight-tree selector

Two precise consumers of the integral source-free reconstruction basis:

- For every oriented forest of every order, nonconstant primitive coordinates begin at source-free core order seven. Source-stack substitutions with a core of at most six vertices preserve every spanning forest parity; the inherited seven-tree witness makes that threshold sharp.
- The eight-vertex internally rooted double-path D has exactly two orientation-parity functions. Of its 128 labelled orientations, 64 are identically odd and 64 have one common nonconstant defect. These give 32 Rédei and 32 nonconstant oriented isomorphism classes. The defect lives entirely at primitive degree seven and is erased by one more isolated-vertex padding.

Read [the proof](DOUBLE_PATH_SELECTOR.md). The source-stack theorem is not arbitrary minor or substitution invariance. This packet does not classify all good trees. Its Hamiltonian-path appendix explicitly identifies an inherited OCF identity, without claiming a new forbidden-value theorem.

The producer performs 2,428 exact checks, including all 128 orientations on each of two fixed hosts, an independent embedding-count method, complete automorphism/orientation orbit counts, and the sharp seven-core witness. The universal assertions come from the proofs, not from that finite host sample.

## Review and exact replay

The [independent proof review](review/DOUBLE_PATH_PROOF_REVIEW.md) passes the all-host selector, complete orientation-class split, degree-seven localization, sharp all-order forest primitive floor, and credited OCF translation. Its separate verifier replays the inherited hosts and exact counts without importing producer code.

From the repository root:

```bash
cd 05-knowledge/results/tournament_double_path_selector_20261019
python3 -B verify_double_path_selector.py > /tmp/double-path-selector.json
python3 -B -O verify_double_path_selector.py > /tmp/double-path-selector-optimized.json
cmp selector-receipt.json /tmp/double-path-selector.json
cmp selector-receipt-optimized.json /tmp/double-path-selector-optimized.json
python3 -B review/verify_double_path_independent.py > /tmp/double-path-independent.json
python3 -B -O review/verify_double_path_independent.py > /tmp/double-path-independent-optimized.json
cmp review/double_path_independent_receipt.json /tmp/double-path-independent.json
cmp review/double_path_independent_receipt_optimized.json /tmp/double-path-independent-optimized.json
```

Both verifiers use only the Python standard library. SHA-256 pins are included in `SHA256SUMS` and the publication manifest.
