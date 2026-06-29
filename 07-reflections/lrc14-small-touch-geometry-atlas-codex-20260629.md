# LRC14 Small-Touch Geometry Atlas Reflection

HYP-3478 clarifies what the six small-touch/no-hard rows are doing.

The important geometric fact is simple:

```text
dead-cover projection edges = 0
```

for every one of the six rows.  HYP-3476 already hinted at this, but the atlas
makes it visible component-by-component.  The dead components are not a weakly
connected projection graph.  They are isolated rank-2 pockets.

Every dead pocket has:

```text
one B0 owner
one B1 owner
a unique label pair in its row
an exact mirror partner
at least two complete E/branch gate touches
```

Five rows have one mirror pair of pockets.  `random_covering_001` has two
mirror pairs.  The owner-pair grammar is:

```text
001: (165,179), (81,153)
039: (63,129)
062: (9,81)
074: (15,99)
086: (133,169)
101: (7,175)
```

This shifts the proof picture.  A bigger gate tuple will not create a
projection edge, because no dead-label support is shared.  The right object is
a singleton-current lemma for a mirror pair of isolated pockets, with the
complete/touching E/branch gates and owner-pair residue/span word retained.

The S319 split is also helpful.  Four rows are branch-unit singleton rows:

```text
001, 062, 086, 101
```

Two rows are S319 delta-sidecar rows:

```text
039, 074
```

So the next proof attempt should not start with all six.  Prove the four
branch-unit rows first as the pure singleton-current packet, then add the
cover-delta sidecar needed for `039` and `074`.

A companion best-touch audit narrows the first packet further.  Within the
four branch-unit rows, `062`, `086`, and `101` have clean best current
`(1,1)`, while `001` has four isolated components and best touching current
`(2,1)`.  The most efficient order is therefore: prove the three clean
best-touch rows, then the asymmetric branch-unit row `001`, then the two
delta-sidecar rows `039` and `074`.

The guardrail stays the same as HYP-3474/HYP-3476: raw row names, gate counts,
or color counts are shadows.  The proof object must retain the dead-pocket
interval, mirror partner, owner pair, residue/span word, and complete gate
touch sidecar.
