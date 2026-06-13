---
id: THM-353
name: staircase-tiling-transport-checksum
status: PROVED
date: 2026-05-30
session: codex-2026-05-30
depends_on:
  - THM-330
  - THM-336
  - THM-351
  - THM-352
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/StaircaseTileModel.lean
  - 04-computation/lean/TournamentH7/TournamentH7/StaircaseBucketTransport.lean
  - 04-computation/lean/TournamentH7/TournamentH7/Verify.lean
results:
  - 05-knowledge/results/lean_staircase_bucket_transport_codex_2026-05-30.out
---

# THM-353: Staircase Tiling Transport Checksum

## Statement

The concrete staircase tiling space

```text
StTiling n := StTile n -> Bool
```

is a finite Boolean cube.  In Lean, `StTile n` is equivalent to the finite
subtype of legal coordinate pairs

```text
{(hi,lo) : Fin n x Fin n // lo.val + 2 <= hi.val}.
```

Therefore the Boolean-cube transport row checksum of THM-352 applies directly
to staircase tilings.

For any finite bucket map

```text
q : StTiling n -> beta
```

and any selected finite family `M` of nonzero staircase masks,

```text
2*internalLineCount_q,M(b) + sum_{c != b} |transportHalf_q,M(b,c)|
  = |q^{-1}(b)| * |M|.
```

This is formalized for three concrete mask families:

- all nonzero staircase masks;
- all single-tile masks;
- the complement mask `allUp`, for `n >= 3`.

It is also formalized for the concrete good-cut quotient

```text
goodCutBucket : StTiling n -> Fin (n+1).
```

The bottom bucket is exactly the all-down tiling condition, and the top bucket
is exactly strong connectivity of the induced staircase tournament.

## Lean Formalization

`TournamentH7.StaircaseTileModel` proves the finite/equality infrastructure:

- `StTile.toGapPair`
- `StTile.ofGapPair`
- `StTile.equivGapPair`
- `StTile.instDecidableEq`
- `StTile.instFintype`
- `StTiling.instFintype`
- `StTiling.instDecidableEq`

`TournamentH7.StaircaseBucketTransport` defines:

- `StTiling.nonzeroMasks`
- `StTiling.singleTileMasks`
- `StTiling.complementMask`
- `StTiling.goodCutBucket`
- `StTiling.topGoodCutBucket`

and proves:

- `StTiling.singleUp_isNonzeroMask`
- `StTiling.allUp_isNonzeroMask_of_three_le`
- `StTiling.unordered_balance_allNonzeroMasks`
- `StTiling.unordered_balance_singleTileMasks`
- `StTiling.unordered_balance_complementMask`
- `StTiling.transport_row_balance_allNonzeroMasks`
- `StTiling.transport_row_balance_singleTileMasks`
- `StTiling.transport_row_balance_complementMask`
- `StTiling.goodCutBucket_eq_zero_iff_all_down`
- `StTiling.goodCutBucket_eq_top_iff_toTournament_SC`
- `StTiling.transport_row_balance_goodCutBucket_allNonzeroMasks`
- `StTiling.transport_row_balance_goodCutBucket_singleTileMasks`
- `StTiling.transport_row_balance_goodCutBucket_complementMask`

`TournamentH7.Verify` audits these theorems.  The transport checks depend only
on Lean foundations (`propext`, `Classical.choice`, `Quot.sound`).

## Proof

First, `StTile n` is sent to its `(hi,lo)` coordinate pair plus the legal gap
proof.  The inverse rebuilds the tile from that subtype element.  These maps
are inverse by cases on the tile/subtype, so `StTile n` inherits finite and
decidable equality instances.

Second, `StTiling n` unfolds to `StTile n -> Bool`, hence it is a finite
Boolean cube.  A mask is nonzero exactly when some tile coordinate is `true`.
Single-tile masks are nonzero by evaluating at their own tile, and the
complement mask is nonzero for `n >= 3` by evaluating at the apex tile.

Third, THM-352 applies to `q : StTiling n -> beta` and the chosen mask family.
The good-cut quotient is a finite target because `goodCutCount <= n - 1`, so
it lives in `Fin (n+1)`.  The existing good-cut and connectivity theorems give
the interpretations of the bottom and top good-cut buckets.

## Use

THM-353 is the concrete checklist for staircase quotient computations.  A
computed transport matrix for good-cut buckets, merged classes, feature
clusters, or any other finite quotient of `StTiling n` must satisfy the row
checksum for the mask family being used.

It also clarifies the formal boundary: THM-351 and THM-352 are fully abstract
Boolean-cube bookkeeping; THM-353 is the bridge proving that actual staircase
tilings are one of those cubes.

The follow-up THM-355 adds the gap layer: the finite good-cut quotient has
image exactly `{0} ∪ {2,...,n-1}` for `n >= 3`, and empty quotient fibers give
zero transport rows and columns.

## Related

- THM-330: strong connectivity iff all base-path cuts are crossed upward.
- THM-336: good-cut bucket spectrum and the missing bucket 1.
- THM-351: Boolean-cube mask bucket balance.
- THM-352: quotient transport row checksum.
- THM-355: quotient gap transport vanishing.
- `05-knowledge/variables/tiling-bucket-balance.md`.
- `07-reflections/staircase-transport-is-boolean-transport.md`.
