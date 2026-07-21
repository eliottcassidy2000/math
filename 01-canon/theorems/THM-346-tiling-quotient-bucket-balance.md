---
id: THM-346
name: tiling-bucket-balance
status: PROVED
date: 2026-05-29
session: kind-pasteur-2026-05-29-S5
scripts:
  - 04-computation/tiling_quotient_bucket_balance_s5.py
results:
  - 05-knowledge/results/tiling_quotient_bucket_balance_s5.out
---

# THM-346: Tiling Bucket Balance for Quotients of the Hypercube

## Statement

Let `Q_m = {0,1}^m` be the tiling hypercube. Let

```text
q : Q_m -> B
```

be any quotient map into a finite set of buckets. For a nonempty set `M` of
nonzero masks in `Q_m`, put an undirected line between `x` and `x xor u` for
each `u in M`.

For a bucket `b in B`, define:

- `B_b = q^{-1}(b)`.
- `self_b(M)` = number of unordered `M`-lines with both endpoints in `B_b`.
- `w_bc(M)` = number of unordered `M`-lines with one endpoint in `B_b` and
  the other in `B_c`, for `c != b`.

Then every bucket satisfies the exact balance law

```text
2*self_b(M) + sum_{c != b} w_bc(M) = |B_b| * |M|.
```

Equivalently, each bucket has exactly `|B_b|*|M|` incident half-lines, with
internal lines consuming two half-lines and cross-bucket lines consuming one.

## Important Instances

The theorem applies to all quotient maps used in the tiling project, including:

- `Q_m -> G_n/Z_2`, the merged tournament isomorphism-class quotient.
- `Q_m -> E_n`, the even-graph quotient.
- Any refinement or coarsening of these quotients, including buckets stratified
  by `H`, score sequence, SC/NS type, spine/ribs/sea region, or projection
  defect profile.

For the merged tournament quotient, every cross-bucket line also has a
class-level geometry:

- `SC-SC` line: spine.
- `SC-NS` line: rib.
- `NS-NS` line: sea.

Thus each merged bucket has a refined balance:

```text
2*self_b + spine_incident_b + rib_incident_b + sea_incident_b = |B_b|*|M|,
```

where terms impossible for the type of `b` are zero.

## Complement-Tiling Parity Corollary

If `|M| = 1`, in particular for the complement-tiling mask that flips all
tiles, then

```text
sum_{c != b} w_bc(M) == |B_b|  (mod 2).
```

So every odd-size bucket has an odd number of cross half-lines under that
involution. This is a cheap parity check for any computed complement-tiling
quotient matrix.

## Proof

Fix a bucket `b`. Count ordered incident pairs `(x,u)` with `x in B_b` and
`u in M`. There are exactly `|B_b|*|M|` such pairs.

Each internal unordered line `{x, x xor u}` with both endpoints in `B_b`
contributes two ordered incident pairs, one from each endpoint. Each cross
line with exactly one endpoint in `B_b` contributes one ordered incident pair.
This gives the formula.

The parity corollary follows when `|M|=1` because the formula becomes

```text
2*self_b + incident_cross_b = |B_b|.
```

Reducing modulo 2 gives `incident_cross_b == |B_b| mod 2`.

## Lean Core

`TournamentH7.BucketBalance` formalizes the oriented finite-set core as
THM-348:

```text
|selfHalf_b(M)| + |crossHalf_b(M)| = |q^{-1}(b)| * |M|.
```

This is the exact half-line conservation law behind the displayed unordered
formula.

S2 adds THM-350, the next abstract Lean layer.  It defines the partner map
`(x,u) -> (step u x,u)`, proves internal half-lines are closed under that map
when the selected moves are involutions, proves fixed-point-free moves have no
self-partnered internal half-lines, and proves the unordered balance whenever
`selfHalf.card` is even:

```text
2*internalLineCount_b(M) + |crossHalf_b(M)| = |q^{-1}(b)| * |M|.
```

A later Codex pass closes the abstract orbit-parity step: Lean now proves that
a finite fixed-point-free involution has even cardinality, and THM-350 derives
the unordered balance directly for fixed-point-free involutive move systems.
THM-351 closes the Boolean-cube specialization: xor by any nonzero mask is a
fixed-point-free involution, so the full tiling-cube balance is now formalized
at the abstract Boolean-mask level.

## Computation

`tiling_quotient_bucket_balance_s5.py` verifies the balance implementation
exactly for `n=3..6` for both `G_n/Z_2` and `E_n`, using the following mask
families:

- `d=1` wiggly masks.
- middle Hamming layers.
- `d=m` complement-tiling mask.
- all nonzero waggly masks.

The verification found zero balance violations and zero complement-parity
violations in all tested cases.

The same run also shows that finite layers are usually not equitable inside a
bucket: individual tilings in the same bucket can have different numbers of
same-bucket neighbors, even though the aggregate balance law is exact. This is
why bucket constraints are safer than assuming quotient regularity.

## Engineering Use

The formula gives a linear checksum for quotient construction. Given bucket
sizes and a proposed weighted quotient matrix, every row must satisfy the
half-line balance. This is useful for:

- validating merged-metagraph and even-graph computations,
- detecting stale complement-tiling versus tournament-complement mistakes,
- building normalized "escape" and "neutrality" features for
  `tournament_tda.py`.

## Related

- MISTAKE-031 and MISTAKE-033: tiling complement is not tournament complement.
- THM-345: merged-bucket parity and Hamming-layer transport constraints.
- THM-348: finite bucket half-line balance (Lean core).
- THM-350: finite unordered bucket balance layer.
- THM-351: Boolean-cube mask bucket balance.
- THM-352: quotient transport row checksum.
- THM-353: staircase tiling transport checksum.
- INV-194: merged tiling bucket constraints.
- INV-236: projection-defect profiles across tournament and even-graph quotients.
- `07-reflections/merged-tiling-bucket-constraints.md`.

## 2026-07-21 Frobenius walk-congruence addendum

The whole-layer mechanism of THM-2041 gives an additional exact transport law
which does **not** require the quotient to be equitable. Let the Boolean tiling
cube be the additive group `G=F_2^m`, let `T_u` denote translation by a mask
`u`, and for any mask set `M` put

```text
A_M = sum_(u in M) T_u.
```

The translations commute. Therefore Freshman's dream gives, over `F_p`,

```text
A_M^p = sum_(u in M) T_u^p.
```

For every odd prime `p`, `T_u^p=T_(pu)=T_u` in the exponent-two group, while
in characteristic `2`, `T_u^2=I`. Hence

```text
A_M^p = A_M          in characteristic p, p odd,          (F1)
A_M^2 = |M| I        in characteristic 2.                 (F2)
```

The `(x,y)` entry of `A_M^r` counts ordered `r`-step mask walks from `x` to
`y`. Thus for odd `p`,

```text
# ordered p-step M-walks x->y = 1_M(x xor y) mod p.        (F3)
```

Define the oriented one-step bucket transport by

```text
W_M(b,c)=sum_(x in B_b, y in B_c) 1_M(x xor y).
```

Summing (F3) over `x` and `y` in arbitrary quotient buckets `B_b,B_c` shows
that the aggregate oriented `p`-step bucket transport is congruent modulo `p`
to `W_M(b,c)`. On a diagonal bucket this convention counts both orientations,
so it is twice an unoriented internal-edge count. The congruence remains true
for every quotient in this theorem even when individual vertices in one
bucket have different neighbor profiles. It is therefore a genuine
congruence beyond the parity row checksum, and it respects rather than evades
the non-equitability warning.

For Tournament Analysis, take mask-walk lengths as vertices and orient an edge
toward the length whose transport matrix retains more exact residue data. The
pairwise observable is equality of all bucket-to-bucket entries modulo `p`.
All lengths `1,p,p^2,...` are congruent and hence tied; the gauge reduces them
to the one-step matrix and the declared tie path is `1 -> p -> p^2 -> ...`.
Equations (F1)--(F3) therefore make the tie-oriented tournament transitive.
This quotient preserves aggregate transport and destroys within-bucket local
degrees; the challenged assumption is that quotient regularity is needed for
prime-step congruences.
