# LRC14 Random031 Topology Atlas, codex-2026-06-29

This pass tried to stop looking at `random_covering_031` as the one stubborn
row and instead ask what shape it is.

The shape is surprisingly crisp.

## Reframe 1: Mirror-Punctured Annulus

The dead-cover projection has four components, largest component one, and no
edges.  The dead components are not a connected bad graph needing a Menger cut;
they are four isolated punctures, paired by mirror:

```text
(0,3), (1,2)
```

Each island has rank `2` and a singleton B0/B1 label pair.  This explains why
pair-current never fixes the HYP-3472 edge-cut failure: there is literally no
edge in the dead projection to cut.  The right local lemma is island current,
not larger gate tuples.

## Reframe 2: Bypassed Saddle Seam

The HYP-3455 max-delta gates on components `(43,54)` are a true mirror seam.
They carry the full seven-owner boundary

```text
(23,45,93,113,147,169,173)
```

but the q=`14V` phase flow never crosses them:

```text
hard_gate_hits=0
```

The same components also carry a lower-delta mirror bypass:

```text
43 branch1 delta 2
54 branch0 delta 2
```

and phase flow hits that bypass twelve times, six on each side.  So the hard
pair is not a sealed wall.  It is a saddle whose legal bypass already lives on
the same two components.  The remaining proof obligation is to show the
seven-owner seam cannot force global obstruction after the bypass is retained.

## Next Proof Move

Promote random031 to a two-layer local theorem:

```text
mirror-punctured islands
  + bypassed saddle seam
  + seven-owner boundary
  => local island current or named owner/two-adic/SPEC debt
```

This should plug into HYP-3473 as a terminal packet next to HYP-3476
pair-current, HYP-3477 hard-orbit discharge, and HYP-3478 small-touch
geometry.
