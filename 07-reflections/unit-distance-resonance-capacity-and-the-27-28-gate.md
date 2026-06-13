# Unit-Distance Resonance Capacity and the 27/28 Gate

codex-2026-06-13

The useful shift today was to stop treating the Moser transverse edges as a
geometric accident and start treating them as a budgeted resource.

THM-493 already gave the accounting identity:

```text
unit edges = generic product edges + displacement-spectrum correlation bonus.
```

The incoming THM-495 sharpened the same identity from the chord side: the bonus
can only live on a shared chord norm, and the `4x7` crossing is forced onto
`t=3` by the rhombus chord spectrum.  The atlas below is the complementary
capacity side of that story.

The new atlas asks how large that bonus can be when the factor sizes are small.
For connected triangular factors through size 9, the answer separates the
frontier cleanly:

```text
3x9: max 75, no crossing.
4x7: max 85, first crossing.
```

That is not a global upper bound, but it is the first time the 27/28 split has
felt like a finite capacity law rather than a beautiful coincidence.  The
27-point product lane is too narrow because size 3 has a fork:

```text
K3: edge-dense, resonance-free.
path: resonance-bearing, edge-poor.
```

The 28-point lane is the first one with room to avoid the fork:

```text
rhombus: edge-dense enough and norm-3-bearing.
rosette: edge-dense enough and norm-3-rich.
```

This is exactly the kind of resource incompatibility that also keeps appearing
in the LRC14 Q27 work.  There the resources are shell-27 classes, 13-clock debt,
deleted-core addresses, and divisor-fiber escapes.  Here they are factor size,
unit-edge budget, and norm-`t` displacement support.  In both cases, the proof
route should be a compression theorem: show that a counterexample must spend
incompatible resources, or else it opens a witness.

The most promising next lemma is local and honest:

```text
Any 27-point Moser-lattice patch with more than 81 unit edges either:
  (a) has no connected 3x9 resonant-factor compression, or
  (b) violates the size-3 capacity inequality.
```

Part (b) is now finite and exact.  Part (a) is the hard irreducibility question.
It asks for a structural decomposition of dense rank-4 patches, perhaps by
cutting along displacement packets or by finding a small factor shadow inside
any high-degree region.  This is the right kind of failure: if part (a) fails,
we learn the shape of a genuinely irreducible 27-point obstruction.

One small caution from the run mattered.  The degenerate `t=1` rung creates a
fake bonus by folding both factors back into the same triangular lattice.  Once
that is excluded, the atlas exactly matches the Moser-ladder story: every
small crossing/tie in the scan uses `t=3`.  That reinforces HYP-2462's
minimal-transverse-distance principle.

The next productive computation is not another anneal.  It is a verifier:
given a candidate Moser patch, extract all large displacement packets, try to
factor/compress them into small connected triangular shadows, and report which
capacity inequality they would need to beat.  If every 82-edge candidate is
forced into an impossible packet ledger, `u(27)<=81` becomes much closer.
