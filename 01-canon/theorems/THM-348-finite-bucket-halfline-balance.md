---
id: THM-348
name: finite-bucket-halfline-balance
status: PROVED
date: 2026-05-30
session: kind-pasteur-2026-05-30-S1
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean
results:
  - 05-knowledge/results/lean_bucket_balance_kind_pasteur_2026-05-30-S1.out
---

# THM-348: Finite Bucket Half-Line Balance

## Statement

Let `alpha` be a finite set, let

```text
q : alpha -> beta
```

be any bucket map, and let `M` be a finite set of moves with an action

```text
step : M -> alpha -> alpha.
```

For a bucket `b`, define:

- `fiber(b) = {x : alpha | q(x)=b}`.
- `incidentHalf(b) = fiber(b) x M`.
- `selfHalf(b) = {(x,u) in incidentHalf(b) | q(step(u,x)) = b}`.
- `crossHalf(b) = {(x,u) in incidentHalf(b) | q(step(u,x)) != b}`.

Then

```text
|selfHalf(b)| + |crossHalf(b)| = |fiber(b)| * |M|.
```

Also, `crossHalf(b)` is empty iff every selected move from `fiber(b)` remains
inside bucket `b`.

## Proof

The proof is pure finite-set bookkeeping.  The set `incidentHalf(b)` is the
cartesian product `fiber(b) x M`, so it has cardinality `|fiber(b)|*|M|`.
The predicate `q(step(u,x)) = b` partitions this product into two disjoint
pieces: internal half-lines and escaping half-lines.  Taking cardinalities
gives the formula.

## Lean Formalization

`TournamentH7.BucketBalance` proves this theorem without project-specific
axioms:

- `BucketBalance.halfLine_balance`
- `BucketBalance.crossHalf_card_eq_zero_iff`
- `BucketBalance.selfHalf_card_le_total`
- `BucketBalance.crossHalf_card_le_total`
- `BucketBalance.selfHalf_card_eq_total_of_crossHalf_zero`
- `BucketBalance.crossHalf_card_eq_total_of_selfHalf_zero`

The `TournamentH7.Verify` audit reports only Lean foundations
(`propext`, `Classical.choice`, `Quot.sound`).

The S2 extension in the same module also proves the unordered layer
(THM-350): the partner map `(x,u) -> (step(u,x),u)` preserves internal
half-lines for involutive moves, is nontrivial for fixed-point-free moves, and
the unordered balance follows for fixed-point-free involutive move systems via
the finite orbit-parity lemma.

## Relation to THM-346

THM-346 is the unordered-line specialization for tiling hypercubes:

```text
2*self_b(M) + incident_cross_b(M) = |q^{-1}(b)| * |M|.
```

THM-348 formalizes the oriented half-line conservation law underneath that
identity. THM-350 formalizes the unordered algebra and the finite
fixed-point-free involution orbit-cardinality theorem. THM-351 recovers full
THM-346 at the Boolean-cube level by specializing to each nonzero mask:

```text
(x,u) <-> (x xor u, u).
```

That pairing turns internal oriented half-lines into unordered internal lines
counted twice, while cross lines are counted once from a fixed bucket.

## Engineering Use

This theorem is a reusable checksum for quotient construction.  Any finite
feature extractor that buckets states and samples moves can audit its rows by
checking:

```text
internal_half_lines + escaping_half_lines = fiber_size * number_of_moves.
```

The theorem is independent of tournaments, so it applies directly to future
`tournament_tda.py` perturbation features, even-graph quotient features, and
generic finite-state testing harnesses.

## Related

- THM-345: merged bucket parity and Hamming-layer transport constraints.
- THM-346: tiling quotient bucket balance.
- THM-350: finite unordered bucket balance layer.
- THM-351: Boolean-cube mask bucket balance.
- INV-194: merged tiling bucket constraints.
- `05-knowledge/variables/tiling-bucket-balance.md`.
