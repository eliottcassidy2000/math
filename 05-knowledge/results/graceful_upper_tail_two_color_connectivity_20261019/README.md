# Arbitrary graceful upper-tail fibers: two-color connectivity

Packet 93, 2026-10-19. The main proof gives a finite-matroid rank-staircase reconfiguration lemma and applies it to all integer-difference upper-tail completion fibers. With r free colors the diameter is at most 2^(r-1) for r>=1; r=0 is a singleton. Tree shapes may change. No graceful-tree existence conjecture or one-color connectivity claim is made.

## Replay

From this directory, using standard-library Python:

    python verify_staircase.py > /tmp/staircase-replay.json
    cmp controls.json /tmp/staircase-replay.json
    python -O verify_staircase.py > /tmp/staircase-replay-O.json
    cmp controls.json /tmp/staircase-replay-O.json
    python check_binary_matroids.py > /tmp/binary-replay.json
    cmp binary_matroid_controls.json /tmp/binary-replay.json
    sha256sum -c SHA256SUMS

The binary-matroid control uses assertions and must run without -O. The graphic verifier uses explicit runtime checks and is also replayed with -O. All proof and control dependencies are co-located; no sibling packet is required for replay.

## Related proved results

- [Leaf-buffer theorem](../graceful_two_color_buffer_20261019/)
- [At most one disconnected color](../graceful_connected_color_contacts_20261019/)
- [At most two disconnected colors](../graceful_two_bad_color_contacts_20261019/)
- [Three bad colors with a rank condition](../graceful_three_bad_rank_contacts_20261019/)
- [Earlier frozen-fiber obstruction](../graceful_fiber_obstructions_20261019.md)

These are contextual links, not runtime dependencies. Earlier special-case bounds can be sharper.
