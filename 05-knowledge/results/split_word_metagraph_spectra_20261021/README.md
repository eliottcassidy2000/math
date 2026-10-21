# Split-word metagraph spectra

Research date: 2026-10-21. Read the [proof and scope](session23-slack-spectrum/README.md) and [independent audit](session23-slack-spectrum-review/AUDIT.md). Status: PROVED relative to the inherited slack-three classification, with separate FINITE-EXACT literal controls.

The packet gives exact split-word path-orbit alphabets and slack-three spectra, a Fibonacci-product bijection, and a literal obstruction to determining wiring from scalar invariants. Two distinct Cartesian-square constructions are kept separate. Degree-support completeness is restricted to slack three; arbitrary strong-atom inventories and full metagraph adjacency remain open.

## Reproduction

From this packet directory, use output paths outside the archive:

    python3 -B session23-slack-spectrum/check.py --receipt /tmp/session23-spectrum-replay.json
    python3 -B -O session23-slack-spectrum/check.py --receipt /tmp/session23-spectrum-replay-optimized.json
    python3 -B session23-slack-spectrum-review/check_literal_ports.py --receipt /tmp/session23-spectrum-review-replay.json
    python3 -B -O session23-slack-spectrum-review/check_literal_ports.py --receipt /tmp/session23-spectrum-review-replay-optimized.json

The original proof, checkers, receipts and source pins are unchanged. Recorded finite controls do not replace the all-order proof and do not establish a full tournament census. Publication is established by the enclosing verified checkpoint.
