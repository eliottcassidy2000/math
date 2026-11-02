# Square-sum recovery and sharp order-ten stability

2026-11-02. Scope: the user's clarified target is to extend the small-order
**weighted stability** range. Raw tournament Hamilton-path counts 7 and 21
were already proved absent at all orders; this session does not rebrand that
old theorem as a finite-order result.

- **PROVED, computer-assisted:** [STABILITY_10.md](STABILITY_10.md) proves the
  sharp coefficient 39/64 for every real weighted tournament of order ten,
  with positive directed-cycle reserves and the fair-pair equality locus.
  The finite reduction checks all 12,005,168 unlabelled simple graphs with
  two exact kernels. The small-order range becomes 3..10, not an all-order
  solution; the inherited large-order range begins at 42.
- **PROVED, elementary synthesis:** [NATIVE_CONTEXTS.md](NATIVE_CONTEXTS.md)
  turns the old abstract boundary-pairing separator into an actual induced
  square-sum extension with unchanged integer source labels. It requires
  exactly two new vertices per cycle connector. Sparse extension is distinct
  from initial-interval coverage.
- **Recovery:** [CATALOG.md](CATALOG.md) and `CATALOG.json` index all 936 matching
  tracked Markdown files at the pinned baseline, comprising 665 distinct
  blobs. This is an inventory, not 936 proof audits. The mechanism-level
  review and its limits are recorded in [RECOVERY.md](RECOVERY.md).

Proof review is same-session, with independently implemented exact counting
paths, not an independent-agent audit or Lean verification. The generator
dependency, full reproduction commands, positive/hostile controls and the
prepublication multiplicity repair are explicit in [AUDIT.md](AUDIT.md).
See [BOARD.md](BOARD.md) and [HANDOFF.md](HANDOFF.md) for the remaining targets.

[INCOMING.md](INCOMING.md) records the mid-session integration: the new Q25
guarded tournament response, an independent check of the50396-prime path,
its singleton scope correction, and the triangle surplus for the next order.
Run `python3 seal.py --check` to verify packet and cited-source integrity;
this is a byte check, not a substitute for the mathematical audit.
