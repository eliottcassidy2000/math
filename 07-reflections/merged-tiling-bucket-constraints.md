# Merged Tiling Bucket Constraints

**Session:** kind-pasteur-2026-05-29-S5
**Computation:** `04-computation/tiling_quotient_bucket_balance_s5.py`
**Results:** `05-knowledge/results/tiling_quotient_bucket_balance_s5.out`
**Canon:** THM-345, THM-346

## The Small Formal Object

THM-345 gives the sharper merged-bucket parity theorem: SC merged buckets
have odd size, NS merged buckets are `2 mod 4`, and every Hamming layer has
ordered transport constraints. The extension in THM-346 is the quotient-level
version: any quotient of the tiling hypercube is a bucket system. Once stated that way,
the basic constraint is almost embarrassingly clean:

```text
2*self_b(M) + incident_cross_b(M) = |bucket_b| * |M|.
```

Here `M` is any set of nonzero masks. A line internal to bucket `b` consumes
two half-lines; a line leaving `b` consumes one. That is the whole theorem.

The reason it matters is that the project has many quotient maps from the same
triangle:

- tilings to merged tournament classes `G_n/Z_2`;
- tilings to even graph classes `E_n`;
- tilings to H-levels, score buckets, SC/NS type, or projection-defect profiles.

Each quotient inherits the same half-line accounting. So any proposed weighted
quotient matrix has row constraints before we ask anything subtler.

## What The Exact Run Adds

The S5 quotient-balance script verified the balance law exactly for `n=3..6`, for both the
merged tournament quotient and the even-graph quotient. It checked `d=1`,
middle layers, the complement-tiling layer `d=m`, and the full waggly graph.
There were zero balance violations and zero complement-parity violations.

Some useful exact totals:

| n | quotient | bucket count | `d=1` self/cross | `d=m` self/cross | all-waggly self/cross |
|---|---|---:|---:|---:|---:|
| 5 | `G_n/Z_2` | 10 | 20 / 172 | 4 / 28 | 264 / 1752 |
| 5 | `E_n` | 7 | 24 / 168 | 6 / 26 | 366 / 1650 |
| 6 | `G_n/Z_2` | 34 | 494 / 4626 | 26 / 486 | 23850 / 499926 |
| 6 | `E_n` | 16 | 556 / 4564 | 98 / 414 | 56686 / 467090 |

The merged tournament cross-lines also decompose into the geometric alignment:

| n | layer | spine SC-SC | ribs SC-NS | sea NS-NS |
|---|---|---:|---:|---:|
| 5 | `d=1` | 124 | 44 | 4 |
| 5 | `d=m` | 16 | 12 | 0 |
| 5 | all | 1108 | 624 | 20 |
| 6 | `d=1` | 276 | 1572 | 2778 |
| 6 | `d=m` | 40 | 156 | 290 |
| 6 | all | 25362 | 188160 | 286404 |

This is a nice sanity check on the project narrative: at `n=5` the spine still
dominates many tiling-level cross-lines, while at `n=6` the sea has already
taken over.

## The Parity Handle

For the complement-tiling involution, `|M|=1`, so every bucket satisfies

```text
incident_cross_b == |bucket_b| mod 2.
```

That is tiny but useful. If a complement-tiling edge table ever violates this
parity, it is almost certainly using tournament complement or class edges
instead of tiling complement lines.

## Non-Equitable Buckets

The aggregate constraints do not mean the quotient is equitable. At `n=6`,
for the merged tournament quotient:

- `d=1`: 18 of 34 buckets have nonconstant internal-neighbor counts.
- middle layer `d=5`: 28 of 34 buckets are non-equitable.
- complement-tiling: 10 of 34 buckets are non-equitable.
- all-waggly: 0 of 34, trivially, because every tiling sees all other tilings.

So the right formal object is not "every tiling in a class behaves the same."
The right object is conserved half-line mass per bucket. That distinction
keeps us from overreading quotient matrices.

## Engineering Reading

For `tournament_tda.py`, THM-345 and THM-346 give two cheap feature families:

- **escape profile:** `incident_cross_b / (|bucket_b|*|M|)`;
- **neutrality profile:** `2*self_b / (|bucket_b|*|M|)`.

These are normalized perturbation-response features. They can be computed for
the tournament quotient, the even-graph quotient, or any cheap coarsening such
as score buckets. The balance equation makes the features self-checking:
escape plus neutrality must be 1.

For systems work, the same equation is a row-sum audit for weighted quotient
matrices. It catches line/edge confusion and complement confusion without
requiring a theorem prover.

## Next Question

Condition this bucket accounting on the principal line:

> Do bucket escape ratios vary monotonically from transitive floor through the
> SC spine toward the regular ceiling, or is the escape pressure primarily a
> sea phenomenon?

This is the spine/ribs/sea version of the projection-defect program. It should
be the next exact computation before attempting `n=7` sampling.
