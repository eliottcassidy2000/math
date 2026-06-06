# LRC14 Apex Opening Modes S678b

S677 gave the right residual:

```text
primitive apex debt => p_0(V,1/14)>0.
```

S678b improves the route by refusing to treat that as one shapeless measure
claim.  If a speed `a` is divisible by `14`, then near the old unit clocks it is
a shutter of radius `1/(14a)`.  Deleting `a` often recovers a one-sided cone at
some `j/14`; if the cone is wider than the shutter, the proof is local and
finished.

That covers `235` of the `414` primitive apex-debt probes in the S677 coherent
carry atlas.

The remaining rows show why the older `n`-clock language was too narrow.  They
do not need a clock point.  They need endpoint owners.  In the full row, every
positive safe interval has two boundary owners; recording whether those owners
are apex or non-apex gives the useful trichotomy:

```text
apex_free       230 rows
one_apex        182 rows
apex_both_only    2 rows
```

After clock priority, the actual certificate routes are:

```text
clock_shutter          235
apex_free_side_door    106
one_apex_hinge          73
```

The two apex-both-only rows are caught by clock shutter.  So the sampled
primitive apex-debt branch has no mystery bucket.

The conceptual gain is small but load-bearing: owner-private deletion becomes
a boundary question.  Instead of trying to prove that the apex is harmless in
bulk, ask which endpoint of a specific safe interval the apex owns.

This suggests a four-lemma package:

```text
clock shutter:  deleted cone wider than apex shutter
side door:      both endpoint owners non-apex
hinge:          one apex owner and one non-apex owner
aperture:       both apex owners, but the apex period leaves an open slit
```

The first lemma is almost tautological after the local derivative calculation.
The last three are the geometric form of the owner ledger.  They also mesh with
the recent metagraph parity lesson: a collapsed scalar graph loses the line
address; here a collapsed `p_0` scalar loses the interval boundary owner.

For the next proof attempt, I would not widen the search first.  I would prove
the four local lemmas in isolation, then build a classifier over normalized
`Res_27` carry-owner rows that forces at least one mode.  The hard case to fear
is a row where every deleted clock cone is swallowed, every safe interval is
apex-owned on both sides, and those apex apertures all close.  S678 finds no
such row in the coherent atlas; the proof should show that such a row would
have to be a scalar wall, hence not primitive apex debt.
