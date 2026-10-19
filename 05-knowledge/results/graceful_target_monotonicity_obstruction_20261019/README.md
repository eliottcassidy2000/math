# Target-monotonicity limitation and attribution

Packet 97, 2026-10-19. The five-state multistar satisfies the rank staircase but requires temporarily removing an already-correct target edge. Exact endpoint distance is four. This is a limitation of target-preserving strategies, not an exponential lower bound.

The bounded integer-difference controls found no corresponding obstruction through m=8 in either one-fixed-edge fibers or three-free-color fibers. These are finite negative checks, not a universal theorem.

Open refinement: is every **nonempty** arbitrary fixed-color graceful completion fiber connected under at-most-two-color moves? An arbitrary fixed forest need not have any completion, unlike an upper tail. Neither existence for arbitrary fixed forests nor this connectivity refinement is asserted.

## Replay

All dependencies are co-located; Python standard library only. From this directory:

    python check_target_monotonicity.py > /tmp/target97.jsonl
    cmp controls.jsonl /tmp/target97.jsonl
    python check_three_free.py > /tmp/three97.jsonl
    cmp three_free_controls.jsonl /tmp/three97.jsonl
    sha256sum -c SHA256SUMS

`PRIMARY_SOURCE_ATTRIBUTION.md` is a focused primary-source comparison for the earlier rank-staircase theorem. It establishes no historical-priority claim.

Related: [rank-staircase and all upper-tail fibers](../graceful_upper_tail_two_color_connectivity_20261019/), [three-bad-color rank criterion](../graceful_three_bad_rank_contacts_20261019/). These links are contextual, not replay dependencies.
