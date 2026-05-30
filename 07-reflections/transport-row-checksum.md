# Transport Row Checksum

Session: codex-2026-05-30

THM-352 is the point where the Boolean-cube bucket balance becomes an
engineering-grade row audit.

THM-351 says the unordered balance is forced by nonzero xor masks:

```text
2*internal + escaping = bucket_size * mask_count.
```

THM-352 names the missing matrix object. For buckets `b,c`,
`transportHalf(b,c)` counts oriented half-lines starting in `b` and landing in
`c`. The target buckets partition all incident half-lines, and deleting the
diagonal target `b` leaves exactly the escaping half-lines:

```text
sum_{c != b} transportHalf(b,c) = escaping_b.
```

So every computed quotient transport matrix must satisfy

```text
2*internal_b + sum_{c != b} transportHalf(b,c)
  = |q^{-1}(b)| * |M|.
```

This is small as a theorem but useful as a boundary marker. The checksum does
not know what the buckets mean. They can be merged tournament classes,
even-graph classes, good-cut buckets, projection-defect buckets, or feature
clusters for TDA. The semantics live in the labels and in the distribution of
off-diagonal mass; the row integrity is pure finite-set bookkeeping.

That gives a clean workflow for the next computational pass:

1. Build a quotient transport matrix.
2. Validate every row with THM-352.
3. Only then study how the off-diagonal mass splits into spine, ribs, sea, or
   whatever geometry the quotient exposes.

The proof also confirms the right abstraction level. We do not need a
tournament-specific theorem to trust row sums. We need tournament-specific
work only to interpret the transport channels after the checksum has passed.

References: THM-350, THM-351, THM-352,
`04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean`.

**Codex-2026-05-30 continuation:** THM-353 now instantiates this row audit
for the actual staircase tiling cube.  `StTile n` has finite coordinate-pair
instances, `StTiling n` is a finite Boolean cube, and the checksum is proved
for all nonzero masks, single-tile masks, the complement mask, and the finite
good-cut quotient `goodCutBucket : StTiling n -> Fin (n+1)`.

**Codex-2026-05-30 gap layer:** THM-355 adds the complementary row/column
fact. If a quotient fiber is empty, the corresponding source row and target
column in `transportHalf` are empty. For `goodCutBucket`, the finite image is
exactly `{0} ∪ {2,...,n-1}` for `n >= 3`, so the sentinel buckets `1` and `n`
are certified zero transport buckets.
