---
id: THM-355
name: quotient-gap-transport-vanishing
status: PROVED
date: 2026-05-30
session: codex-2026-05-30
depends_on:
  - THM-336
  - THM-348
  - THM-353
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean
  - 04-computation/lean/TournamentH7/TournamentH7/StaircaseBucketTransport.lean
  - 04-computation/lean/TournamentH7/TournamentH7/Verify.lean
results:
  - 05-knowledge/results/lean_bucket_gap_transport_codex_2026-05-30.out
---

# THM-355: Quotient Gap Transport Vanishing

## Statement

Let `q : alpha -> beta` be a finite quotient map.  If a bucket `b` has empty
fiber,

```text
q^{-1}(b) = empty,
```

then every incident half-line row out of `b` is empty.  In particular:

```text
incidentHalf(q,b) = empty
selfHalf(q,b)     = empty
crossHalf(q,b)    = empty
transportHalf(q,b,c) = empty  for every c.
```

Likewise, if a target bucket `c` has empty fiber, then no transport row can
land in `c`:

```text
transportHalf(q,b,c) = empty  for every b.
```

Equivalently, the corresponding matrix entries, source-row sum, and
target-column sum all have cardinality zero.

For the concrete good-cut quotient

```text
goodCutBucket : StTiling n -> Fin (n+1),
```

and `n >= 3`, the finite image is exactly

```text
{0} union {2,3,...,n-1}.
```

Thus good-cut buckets `1` and `n` are genuine finite quotient gaps: their
fibers are empty, their source rows are empty, and their target columns are
empty.

## Lean Formalization

`TournamentH7.BucketBalance` proves the generic finite quotient facts:

- `BucketBalance.fiber_eq_empty_iff`
- `BucketBalance.incidentHalf_eq_empty_of_fiber_eq_empty`
- `BucketBalance.selfHalf_eq_empty_of_fiber_eq_empty`
- `BucketBalance.crossHalf_eq_empty_of_fiber_eq_empty`
- `BucketBalance.transportHalf_eq_empty_of_source_fiber_eq_empty`
- `BucketBalance.transportHalf_eq_empty_of_target_fiber_eq_empty`
- `BucketBalance.transportHalf_card_eq_zero_of_source_fiber_eq_empty`
- `BucketBalance.transportHalf_card_eq_zero_of_target_fiber_eq_empty`
- `BucketBalance.sum_transportHalf_card_eq_zero_of_source_fiber_eq_empty`
- `BucketBalance.sum_transportHalf_card_eq_zero_of_target_fiber_eq_empty`

`TournamentH7.StaircaseBucketTransport` proves the good-cut quotient image and
gap instances:

- `StTiling.goodCutBucket_image_iff`
- `StTiling.goodCutBucket_ne_one`
- `StTiling.goodCutBucket_ne_overTop`
- `StTiling.goodCutBucket_fiber_one_eq_empty`
- `StTiling.goodCutBucket_fiber_overTop_eq_empty`

`TournamentH7.Verify` audits both the generic and staircase-specific theorems.
They depend only on Lean foundations.

## Proof

The generic proof is direct finite-set bookkeeping.  Membership in
`fiber q b` is equivalent to `q x = b`.  If no such `x` exists, then the
product defining `incidentHalf q moves b` is empty, and all filters of that
product are empty.  If a target bucket `c` has empty fiber, then any
half-line landing in `c` would exhibit an element of `q^{-1}(c)`, a
contradiction.

For the good-cut quotient, THM-336 already proves the exact natural-number
spectrum

```text
exists u, u.goodCutCount = r
  iff r = 0 or (2 <= r and r <= n-1).
```

THM-353 packages `goodCutCount` as the finite quotient `goodCutBucket`.  The
finite image theorem is obtained by transferring the THM-336 spectrum through
`Fin.val`.  The bucket `1` is excluded by `goodCutCount_ne_one`; the bucket
`n` is excluded by `goodCutCount <= n-1` when `n >= 1`.

## Use

This theorem gives a precise transport meaning to quotient gaps.  A missing
bucket is not only absent from a spectrum; it is a zero row and a zero column
in every transport matrix built from the quotient.

That is useful in three different directions:

- good-cut computations can include buckets `1` and `n` as explicit sentinel
  zero rows/columns;
- quotient transport audits can distinguish impossible buckets from merely
  low-mass buckets;
- residue/projection studies can ask not only which fibers vanish, but which
  rows and columns become zero after a projection.

## Related

- THM-336: good-cut bucket structure.
- THM-348: finite bucket half-line balance.
- THM-353: staircase tiling transport checksum.
- HYP-1783: quotient gap/residue principle.
- `07-reflections/fiber-gaps-and-residue-boundaries.md`.
