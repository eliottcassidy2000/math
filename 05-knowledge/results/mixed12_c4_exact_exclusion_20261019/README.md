# Frozen mixed-(1,2), c=4 exclusion packet

Date: 2026-10-19. Exact, self-contained, standard-library verification.

## Result

Read `mixed12_c4_exact_exclusion_20261019.md` for the theorem and proof. The selected maximum-span majority-root morphology with one lower slave, two upper slaves, and four common neighbors is impossible in a seven-thin poset with no balanced pair. The order-dual branch is included. Other seven-thin morphologies and the general 1/3–2/3 conjecture are outside this result.

Exact coverage:

- 96,428 canonical cores before degree pruning;
- 79,121 cores satisfying every core degree bound;
- 56,630 cores excluded by the separately scoped analytic witnesses;
- all 22,491 residual cores covered by explicit positive integer certificates;
- 4,253,947 admissible residual majority orders;
- 69,126 stored weighted certificates, all with one or two rows;
- zero linear-programming calls in the complete construction.

The earlier c<=3 packet is linked from the proof at the sibling repository path `../mixed12_witness_migration_20261019/mixed12_witness_migration_20261019.md`. The c=4 replay does not depend on that sibling packet.

## Reproduce

From this directory:

```sh
python -O verify_frozen_c4_packet_20261019.py
```

No third-party packages, network access, external data, or prior sample receipts are required. The verifier rebuilds the core catalogue, checks every analytic certificate, recomputes every retained-generator margin from linear extensions, and verifies coverage of every admissible residual majority order. It also checks all package hashes when `manifest.json` is present. Checks are explicit exceptions and remain active under Python optimization. Tested with Python 3.12.14.

The command takes several minutes and prints progress every 1,000 residual cores. It writes `replayed_c4_controls_20261019.json`. Compare its mathematical counts and hashes with `frozen_c4_replay_receipt_20261019.json`. The optional `--limit N` performs only a clearly labeled prefix control and is not a full proof replay.

## Contents and provenance

- `verify_frozen_c4_packet_20261019.py`: final exact replay.
- `verify_saturated_near_twin_lemma_20261019.py`: bundled standard-library catalogue/extension routines used by the replay; its additional specialized audit is not needed for the main command.
- `c4_combined_analytic_coverage_20261019.json.gz`: every analytic witness, its precise logical scope, and the residual catalogue indices.
- `normalized_residual_certificates_20261019.jsonl.gz`: self-contained weighted certificates for every residual core.
- `one_defect_exchange_reductions_20261019.md`, `max_span_directed_defect_reduction_20261019.md`, and `seven_thin_saturation_increment_lemma_20261019.md`: the structural proofs used in the theorem. Each note retains the scope of its own sufficient lemma; the complete combination is in the main theorem note.
- Construction, independent-replay, repacking, and final replay receipts: exact evidence and source hashes. The construction reused earlier checks, while the independent replay supplied 85 additional explicit coefficient certificates to remove all reliance on omitted prior sample coefficients.
- `c4_whole_core_fixed_certificates_20261019.json.gz`: an optional compact explanatory view of 14,554 already-verified fixed-orientation certificates. It contains 11,580 same-common bridges, 1,231 cross-common bridges, and 1,743 single rows. The remaining 7,937 cores retain the full stored list; this does not assert that a fixed certificate is impossible for them.
- `manifest.json`: SHA-256 and byte size of each packaged input, excluding the manifest itself and newly generated replay output.

The normalized stream was losslessly recompressed from 7,369,534 bytes to 1,533,074 bytes. Its decompressed content is unchanged, with SHA-256 `69f3a5ea1c426bd1345da0430ff18226f9b1d224273dc65b9ab664ea2cf7a371`.

The complete construction, independent integer replay, and final compact-packet replay are distinct recorded passes. The final replay uses only the files in this directory. No historical novelty claim is made.
