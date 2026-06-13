---
id: THM-352
name: quotient-transport-row-checksum
status: PROVED
date: 2026-05-30
session: codex-2026-05-30
depends_on:
  - THM-348
  - THM-350
  - THM-351
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean
  - 04-computation/lean/TournamentH7/TournamentH7/Verify.lean
results:
  - 05-knowledge/results/lean_transport_row_bucket_balance_codex_2026-05-30.out
---

# THM-352: Quotient Transport Row Checksum

## Statement

Let `alpha` be finite, let

```text
q : alpha -> beta
```

be a finite bucket map, and let `moves` be a finite set of moves acting by

```text
step : move -> alpha -> alpha.
```

For buckets `b,c`, define the oriented transport entry

```text
transportHalf(b,c) =
  {(x,u) : q(x)=b, u in moves, q(step(u,x))=c}.
```

Then the off-diagonal transport row is exactly the escaping half-line count:

```text
sum_{c != b} |transportHalf(b,c)| = |crossHalf(b)|.
```

Consequently, when the selected moves are fixed-point-free involutions,

```text
2*internalLineCount(b) + sum_{c != b} |transportHalf(b,c)|
  = |q^{-1}(b)| * |moves|.
```

For finite Boolean cubes with nonzero xor masks, this row checksum applies
directly. THM-353 specializes this statement to the concrete staircase tiling
cube `StTiling n`. THM-355 records the gap complement: empty source fibers
give zero transport rows and empty target fibers give zero transport columns.

## Lean Formalization

`TournamentH7.BucketBalance` defines:

- `BucketBalance.transportHalf q step moves b c`

and proves:

- `BucketBalance.transportHalf_self`
- `BucketBalance.sum_transportHalf_card_eq_incidentHalf_card`
- `BucketBalance.sum_transportHalf_card_offdiag_eq_crossHalf_card`
- `BucketBalance.transport_row_balance_of_even_selfHalf`
- `BucketBalance.transport_row_balance_of_involutive_fixedPointFree`
- `BucketBalance.transport_row_balance_boolCube_masks`

`TournamentH7.Verify` audits the off-diagonal transport identity, the
fixed-point-free involutive row checksum, and the Boolean-cube xor-mask row
checksum. These audits depend only on Lean foundations.

## Proof

The sets `transportHalf(b,c)` partition the incident half-lines from bucket
`b` by their target bucket `c`. Removing the diagonal bucket `b` from that
partition leaves exactly those incident half-lines whose target is not `b`,
which is the definition of `crossHalf(b)`.

THM-350 already proves that fixed-point-free involutive move systems have
even internal half-line cardinality and therefore satisfy

```text
2*internalLineCount(b) + |crossHalf(b)| = |q^{-1}(b)| * |moves|.
```

Substituting the off-diagonal transport sum for `|crossHalf(b)|` gives the
row checksum.

THM-351 supplies the Boolean-cube specialization because every nonzero xor
mask is a fixed-point-free involution.

## Use

This theorem turns the bucket-balance law into a concrete matrix audit. For
any computed quotient transport matrix, each row must satisfy the displayed
identity. This is the formal checksum needed for merged-tournament,
even-graph, good-cut, projection-defect, and TDA feature transports.

## Related

- THM-346: tiling quotient bucket balance.
- THM-348: finite bucket half-line balance.
- THM-350: finite unordered bucket balance layer.
- THM-351: Boolean-cube mask bucket balance.
- THM-353: staircase tiling transport checksum.
- THM-355: quotient gap transport vanishing.
- `05-knowledge/variables/tiling-bucket-balance.md`.
- `07-reflections/boolean-cube-balance-as-checksum.md`.
