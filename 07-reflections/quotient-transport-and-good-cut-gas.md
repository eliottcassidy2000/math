# Quotient Transport and the Good-Cut Gas

**Session:** opus-2026-05-29-S15

The good-cut story split into two complementary formal objects today.

First, the interval-union recurrence says the raw tiling cube sees `g` as a
one-dimensional gas on the cut path.  A tiling chooses intervals of length at
least two; the good cuts are their union; connected runs contribute independent
cover weights `c_L`.

Second, the quotient-transport theorem says every quotient bucket has a
conserved ordered half-line budget:

```text
stay + escape = fiber_size * move_count.
```

The new computation puts these together.  Treat the good-cut bucket as the
source coordinate and decompose every outgoing half-line by merged tournament
geometry.  The striking pattern is:

```text
Delta g != 0  =>  merged tournament class changes
```

through every exact scan n=3..6 and every selected mask family in S15.  The
`bad` column in `goodcut_transport_excess_s15.out` is always zero.

This makes `g` feel less like a statistic and more like a quotient-compatible
height: if the height changes, the transport has crossed a real class boundary.

There is also a geometry gradient.  At n=6, ordinary wiggly half-lines with
`|Delta g|=1,2,3` are mostly sea:

```text
|Delta g|=1: sea 73.8%
|Delta g|=2: sea 67.7%
|Delta g|=3: sea 72.7%
```

but the extreme all-down-to-top jump `|Delta g|=5` is pure spine.  The bottom
tiling sees the principal line; the middle buckets see the sea.

That suggests a useful proof shape.  The interval gas governs which heights are
possible and how many tilings live there.  Quotient transport governs how those
heights leak across the merged metagraph.  The spine/ribs/sea decomposition is
not an extra layer on top; it is the escape geometry of the gas.

Engineering reading: `g` should not only be a scalar feature.  A practical
Tournament TDA extractor should include ordered transport features:

- `g` bucket;
- neutral/escape ratio at a few move families;
- `Delta g` histogram;
- spine/ribs/sea or projection-defect profile conditioned on `Delta g`.

Those features are normalized by the Lean balance theorem, so they come with a
built-in row-sum checksum.

**Codex-2026-05-30 update:** THM-352 now formalizes the target-bucket row
version of that checksum in Lean. The off-diagonal sum of
`transportHalf(b,c)` is exactly the escaping half-line term, and the
Boolean-cube xor-mask specialization gives the matrix audit directly.

**Codex-2026-05-30 continuation:** THM-353 pushes the audit into the concrete
staircase model. The good-cut height is now packaged as a finite Lean quotient
`goodCutBucket : StTiling n -> Fin (n+1)`, with bottom equal to all-down and
top equal to strong connectivity of the induced tournament. The single-tile,
all-nonzero, and complement-mask transport rows for this quotient are all
machine-checked instances of the same checksum.

**Codex-2026-05-30 gap update:** THM-355 proves that the finite image of
`goodCutBucket` is exactly `{0} ∪ {2,...,n-1}` for `n >= 3`. So the missing
height `1` is not merely a spectrum fact; it is a zero row and zero column in
the transport matrix. The extra sentinel bucket `n` in `Fin (n+1)` is also a
certified gap.
