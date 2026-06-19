# LRC14 Signed-Mass Sequence Spine

Codex 2026-06-19 S13.

The user prompt was the right diagnosis: the large absolute mass is real, but
it is mostly sitting on boundary/cusp faces where the signed seven-sector
kernel almost cancels.  The proof route should therefore stop asking for a
better absolute Minkowski volume and start asking for the small recursive
sequences that survive quotienting.

The new script `04-computation/lrc14_signed_mass_sequence_spine_codex.py`
collects those sequences.

The first spine is the THM-538 exact-support floor:

```text
0, 0, 0, 0, 0, 0, 1.
```

That is the main structural reason the lower-dimensional cusp mass should be
treated suspiciously.  Proper faces can dominate absolute relation counts, but
they are not the live six-body term after signed inclusion-exclusion.

The second spine is the normalized residue coefficient decay:

```text
d=6..13 max |C_d| =
0.643084862, 0.321542431, 0.168426988, 0.091869266,
0.0519342586, 0.0303106082, 0.0182058608, 0.0112237893.
```

This extends S11/S12 from `d=6..9` to `d=13`.  The decay is strong, but the
one-coordinate marginals remain nonzero:

```text
0.384780076, 0.193862704, 0.107446515, 0.0637696096,
0.0397214435, 0.0256286869, 0.0169882801, 0.011508678.
```

So the missing proof is not a coordinate-by-coordinate cancellation.  It is
the coimage of a relation hyperplane after residue-addressing modulo `7`.

The third spine is the S12 cusp-ratio sequence:

```text
final abs/signed ratios:
29.0667, 217.015, 58.2626, 1118.76, 447.629

one-face boundary ratios:
202.386, 3484.64, 13.3239, 34.8034, 112.732.
```

The one-face integer relation counts are huge, up to `5,291,542`, but their
signed reciprocal totals are tiny.  This is the exact place where absolute
mass and signed mass split.

The other two spines are companion finite ledgers:

- HYP-2598 survivor sequence
  `1,11,47,109,156,146,91,37,9,1,0,0,0,0`;
- HYP-2608 empty-window degree sequence `4,3,3,2,1`.

Together they say the LRC14 proof is not one monolithic inequality.  It is a
stack of sequence-level collapses: support floor, residue decay, cusp
cancellation, center-survivor classification, and moment-degree drop.

Category-theory reading: the absolute relation count is pre-quotient data.  The
signed mass is the coimage after lower support faces are killed and residue
characters are retained.  Number-theory reading: `C_d` is a finite mod-7
character coefficient, and the cancellation belongs to reciprocal sums on
`sum e_i n_i=0`, not to independent residue marginals.

The live target is now very crisp:

```text
finite low-height wall ledger
+ residue-addressed signed reciprocal hyperplane tail
+ no-scale cluster quotient
< cap margin.
```

LRC(14) is still open.  But the "many relevant sequences" request paid off:
the residual now looks like a Dedekind/cotangent tail over a mod-7 coimage,
with the false absolute mass removed from the center of attention.
