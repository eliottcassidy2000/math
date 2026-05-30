---
id: THM-351
name: boolean-cube-mask-bucket-balance
status: PROVED
date: 2026-05-30
session: opus-2026-05-30-S1
depends_on:
  - THM-346
  - THM-348
  - THM-350
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean
  - 04-computation/lean/TournamentH7/TournamentH7/Verify.lean
results:
  - 05-knowledge/results/lean_boolcube_bucket_balance_opus_2026-05-30-S1.out
  - 05-knowledge/results/lean_verify_boolcube_bucket_balance_opus_2026-05-30-S1.out
---

# THM-351: Boolean-Cube Mask Bucket Balance

## Statement

Let `Q = {0,1}^I` be a finite Boolean cube, let

```text
q : Q -> B
```

be any finite bucket map, and let `M` be any finite set of nonzero masks
`u : I -> Bool`. The mask action is coordinatewise xor:

```text
step(u,x) = x xor u.
```

For every bucket `b`, the unordered bucket balance holds:

```text
2*internalLineCount_b(M) + |crossHalf_b(M)| = |q^{-1}(b)| * |M|.
```

Thus the full Lean bridge from the oriented half-line theorem to the tiling
hypercube instance of THM-346 is complete: nonzero Boolean masks are
fixed-point-free involutions.

## Lean Formalization

`TournamentH7.BucketBalance` defines:

- `BucketBalance.BoolCube index := index -> Bool`
- `BucketBalance.xorMask u x := fun i => Bool.xor (x i) (u i)`
- `BucketBalance.IsNonzeroMask u := ∃ i, u i = true`

and proves:

- `BucketBalance.xorMask_involutive`
- `BucketBalance.xorMask_fixedPointFree_of_nonzero`
- `BucketBalance.unordered_balance_boolCube_masks`

The final theorem specializes THM-350:

```text
2 * internalLineCount q xorMask moves b
  + (crossHalf q xorMask moves b).card
= (fiber q b).card * moves.card
```

whenever every selected mask is nonzero.

## Proof

For a fixed mask `u`, coordinatewise xor is involutive because

```text
(x_i xor u_i) xor u_i = x_i
```

for each coordinate `i`.

If `u` is nonzero, choose a coordinate `i` with `u_i=true`. Then
`(x xor u)_i` is the negation of `x_i`, so `x xor u ≠ x`. Hence every
selected nonzero mask acts as a fixed-point-free involution on the Boolean
cube.

THM-350 applies directly to the move system `step(u,x)=x xor u`, yielding the
unordered bucket balance.

## Consequence

The Lean status of THM-346 is now closed at the abstract Boolean-cube level.
For staircase tilings, instantiate `I = StTile n`; then `StTiling n` is
definitionally `StTile n -> Bool`, so the tiling-mask action is exactly the
Boolean-cube mask action.

The follow-up THM-352 refines this into a target-bucket transport-row
checksum: the off-diagonal sum of oriented lines from bucket `b` to all
`c != b` equals the `crossHalf` term, so the same Boolean-cube proof audits
computed quotient transport matrices row by row. THM-353 then instantiates
the Boolean-cube theorem for the concrete staircase tiling cube `StTiling n`.

The remaining work around bucket transport is no longer the balance theorem
itself. It is structural and computational: understanding which quotient
buckets exchange mass, how that mass decomposes into spine/ribs/sea geometry,
and how to turn normalized escape/neutrality profiles into engineering
features.

## Related

- THM-346: tiling quotient bucket balance.
- THM-348: finite bucket half-line balance.
- THM-350: finite unordered bucket balance layer.
- THM-352: quotient transport row checksum.
- THM-353: staircase tiling transport checksum.
- `05-knowledge/variables/tiling-bucket-balance.md`.
- `07-reflections/unordered-bucket-balance-orbits.md`.
- `07-reflections/boolean-cube-balance-as-checksum.md`.
