# A Formalized Staircase Calculus for Tournament Tiling

**A condensed paper on good cuts, strong connectivity, and bucket transport**

**Status:** research draft
**Date:** 2026-05-30
**Lean development:** `04-computation/lean/TournamentH7`

## Abstract

This paper condenses a recent formalization program around the fixed-base
staircase model of tournaments.  The main object is a binary tiling of the
non-consecutive pairs above the base Hamiltonian path

```text
n-1 -> n-2 -> ... -> 1 -> 0.
```

Each upward tile crosses an interval of base-path cuts.  The set of crossed
cuts, called the good-cut set, gives a height statistic

```text
g(b) = |goodCuts(b)|.
```

The Lean development now proves three connected facts.

First, the good-cut statistic has a complete abstract spectrum:

```text
g(b) is either 0 or lies in {2, ..., n-1};
for n >= 3 every such value is realized.
```

Thus bucket 1 is the only interval-geometric obstruction.

Second, the concrete tiling has been connected back to tournaments.  Lean now
constructs the tournament induced by a staircase tiling, proves it has the
base path, proves concrete good cuts are exactly the crossing-upward cuts used
in the strong-connectivity cut theorem, and obtains

```text
g(b) = n-1  iff  the induced tournament is strongly connected.
```

Third, quotient transport out of the tiling cube has been isolated into a
finite bucket calculus.  For any finite bucket map and finite move set, the
oriented half-lines from a bucket split into internal and escaping half-lines.
If the moves are fixed-point-free involutions, internal half-lines pair into
two-element orbits, yielding the unordered balance

```text
2 * internalLineCount + |crossHalf| = |fiber| * |moves|.
```

The corresponding Lean theorems are axiom-free apart from Lean foundations.
The Boolean-cube mask specialization is now also formalized: xor by a nonzero
mask is a fixed-point-free involution, so finite Boolean cube quotients satisfy
the unordered bucket balance directly.

## 1. Introduction

The project studies finite tournaments through three intertwined lenses:

1. Hamiltonian path counts and forbidden H-values.
2. A fixed-base tiling model for tournaments with a chosen Hamiltonian path.
3. Quotients of the tiling cube by tournament isomorphism and complement.

The flagship Lean theorem in `TournamentH7` is the proof-modulo-cited-axioms
statement that no tournament has exactly seven Hamiltonian paths:

```lean
Tournament.H_ne_seven
```

This paper does not re-prove that theorem.  Instead, it records the
formalized staircase calculus that has grown around the tiling model.  The
calculus is useful because it takes objects that first appeared as pictures,
coordinates, or computations and turns them into reusable Lean statements.

The central theme is that much of the observed structure is not accidental
enumeration.  It comes from two elementary but rigid mechanisms:

1. A tile covers an interval of cuts of length at least two.
2. A reversible perturbation contributes two internal half-lines for every
   unordered internal line.

The first mechanism produces the good-cut spectrum and connects its top
bucket to strong connectivity.  The second mechanism produces quotient bucket
balance and explains why transport across quotient classes has forced row
sums and parity constraints.

The recent formalizations close two conceptual loops:

- The good-cut height is not merely a coordinate artifact.  Its top level is
  exactly the strongly connected stratum of the concrete staircase tournament
  family.
- The unordered bucket balance is not a special trick of tiling quotients.  It
  is an instance of finite half-line conservation plus fixed-point-free
  involution orbit parity.

## 2. The Fixed-Base Staircase Model

Fix `n`.  The staircase model starts with the directed base path

```text
n-1 -> n-2 -> ... -> 1 -> 0.
```

For every non-consecutive pair of vertices `hi, lo` with

```text
lo + 2 <= hi,
```

there is a binary tile.  The default, or down, orientation sends the larger
vertex downward:

```text
hi -> lo.
```

The upward orientation reverses that non-consecutive arc:

```text
lo -> hi.
```

In Lean this coordinate is represented by

```lean
StTile n
StTiling n
```

where a tiling is a Boolean-valued function on tiles.

The concrete arc relation was formalized in
`TournamentH7.StaircaseTileModel`:

```lean
StTiling.arc
StTiling.toTournament
StTiling.toTournament_hasBasePath
```

The induced tournament has:

- no loops;
- exactly one directed arc between every pair of distinct vertices;
- the fixed base path `i+1 -> i`;
- non-consecutive arcs determined by the corresponding tile bit.

This is the bridge that lets coordinate statements about tilings become
statements about tournaments.

## 3. Good Cuts

For `1 <= k < n`, the cut `k` separates the vertices into

```text
lower side: {0, ..., k-1}
upper side: {k, ..., n-1}.
```

A tile `(hi, lo)` crosses the cut `k` exactly when

```text
lo < k <= hi.
```

Thus a tile contributes the interval

```text
{lo+1, lo+2, ..., hi}
```

to the set of cuts.  Since tiles are non-consecutive, every such interval has
length at least two.

For a tiling `b`, define:

```lean
StTiling.IsGoodCut b k
StTiling.goodCuts b
StTiling.goodCutCount b
```

The informal meaning is:

```text
k is good if at least one upward tile crosses k.
```

The count

```text
g(b) = b.goodCutCount
```

is the good-cut height of the tiling.

## 4. The Good-Cut Spectrum

The first formalized layer is pure interval geometry.

### Theorem 4.1: Bucket 1 is impossible

For every staircase tiling,

```lean
StTiling.goodCutCount_ne_one
```

proves

```text
goodCutCount(b) != 1.
```

The proof is short and robust.  If one upward tile exists, it crosses both
`lo+1` and `lo+2`, so at least two good cuts appear.  If no upward tile exists,
there are no good cuts.  There is no way to obtain exactly one.

Lean also records the stronger dichotomy:

```lean
StTiling.goodCutCount_eq_zero_or_two_le
```

which states

```text
goodCutCount(b) = 0  or  2 <= goodCutCount(b).
```

### Theorem 4.2: Bucket 0 is exactly all-down

Lean proves:

```lean
StTiling.goodCutCount_eq_zero_iff_all_down
```

Equivalently,

```text
g(b) = 0
iff
every tile is down.
```

This pins the bottom of the height function to the transitive all-down
staircase.

### Theorem 4.3: Exact abstract spectrum

For `n >= 3`, Lean proves:

```lean
StTiling.goodCutCount_spectrum
```

with the content

```text
(exists b, g(b) = r)
iff
r = 0 or (2 <= r and r <= n-1).
```

The constructive direction is explicit.  Bucket 0 is realized by the all-down
tiling.  For any `2 <= r <= n-1`, choose the single upward tile from `0` to
`r`; it crosses exactly the interval `{1, ..., r}` and realizes `g=r`.

This theorem is important because it separates interval geometry from
tournament structure.  If a later quotient, H-layer, or isomorphism class
misses some height value other than 1, the obstruction is not caused by the
good-cut coordinate itself.

## 5. From Top Good-Cut Bucket to Strong Connectivity

The good-cut statistic was first defined in tiling coordinates.  A separate
theorem, the SC cut theorem, characterizes strong connectivity for base-path
tournaments:

```lean
thm330_SC_iff_all_cuts_crossing
```

It says that a base-path tournament is strongly connected iff every legal cut
has a crossing-upward non-consecutive arc.

The recent formalization supplies the missing dictionary between the concrete
tiling and this abstract tournament theorem.

### Theorem 5.1: A concrete tiling induces a base-path tournament

Lean constructs:

```lean
StTiling.toTournament
```

and proves:

```lean
StTiling.toTournament_hasBasePath
```

This is not just a wrapper.  The proof checks the concrete `Bool` arc relation
case-by-case:

- self-arcs are false;
- consecutive pairs follow the base path;
- non-consecutive pairs are oriented by the tile bit;
- exactly one direction holds between distinct vertices.

### Theorem 5.2: Good cuts are crossing-upward cuts

Lean proves:

```lean
StTiling.isGoodCut_iff_crossesUpward_toTournament
```

In words:

```text
StTiling.IsGoodCut b k
iff
CrossesUpward b.toTournament k.
```

The proof is a direct coordinate translation.  An upward tile `(hi, lo)` that
crosses `k` gives the arc `lo -> hi` across the cut.  Conversely, any
crossing-upward non-consecutive arc in the induced tournament is generated by
the corresponding upward tile.

### Theorem 5.3: Top bucket iff strong connectivity

Combining the good-cut dictionary with the SC cut theorem gives:

```lean
StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected
```

That is,

```text
g(b) = n - 1
iff
IsStronglyConnected b.toTournament.
```

The top good-cut bucket is therefore not just "all cuts are covered" as a
coordinate condition.  It is exactly the strongly connected stratum of the
concrete staircase tournament family.

Lean also packages extremal witnesses:

```lean
StTiling.allUp_toTournament_stronglyConnected
StTiling.allDown_toTournament_not_stronglyConnected
```

For `n >= 3`, all-up is strongly connected.  For `n >= 2`, all-down is not.

## 6. Quotient Bucket Balance

The second formalized layer is a finite transport theorem.  It abstracts away
from tournaments almost entirely.

Let:

```text
alpha  = a finite population,
beta   = a set of buckets,
q      : alpha -> beta,
moves  = a finite set of moves,
step   : move -> alpha -> alpha.
```

For a bucket `b`, define its fiber:

```text
fiber(b) = {x in alpha : q(x) = b}.
```

An oriented half-line from bucket `b` is a pair `(x,u)` with `x` in the fiber
and `u` a selected move.  It is internal if `step(u,x)` remains in the same
bucket, and escaping otherwise.

Lean formalizes these as:

```lean
BucketBalance.fiber
BucketBalance.incidentHalf
BucketBalance.selfHalf
BucketBalance.crossHalf
```

### Theorem 6.1: Oriented half-line conservation

Lean proves:

```lean
BucketBalance.halfLine_balance
```

with statement

```text
|selfHalf| + |crossHalf| = |fiber| * |moves|.
```

This is the finite row-conservation law: every incident half-line either stays
inside the bucket or leaves it.

The theorem is deliberately generic.  It does not know about tournaments,
tilings, isomorphism, complements, or H-values.  It is just finite set
bookkeeping, and its genericity is a feature.

### Theorem 6.2: Internal half-lines pair under an involutive move

For an oriented internal half-line `(x,u)`, define its partner:

```text
pairHalf(step)(x,u) = (step(u,x), u).
```

Lean proves:

```lean
BucketBalance.pairHalf_pairHalf
BucketBalance.pairHalf_mem_selfHalf
BucketBalance.pairHalf_ne_of_fixedPointFree
```

If each selected move is involutive, pairing twice returns to the starting
half-line.  If the move is fixed-point-free, a half-line is never its own
partner.

### Theorem 6.3: Fixed-point-free involutions have even finite support

The recent orbit formalization proves:

```lean
BucketBalance.even_card_of_fixedPointFree_involutiveOn
```

In words:

```text
If a finite set is closed under f,
f(f(x)) = x on the set,
and f(x) != x on the set,
then the set has even cardinality.
```

The proof is by strong induction on the cardinality of the finite set.  If the
set is nonempty, choose `x`, let `y=f(x)`, remove both `x` and `y`, prove the
remaining set is still closed under `f`, and recurse.  This is the Lean
version of the elementary fact that a fixed-point-free involution partitions a
finite set into two-element orbits.

Applying this generic theorem to `selfHalf` gives:

```lean
BucketBalance.selfHalf_card_even_of_involutive_fixedPointFree
```

### Theorem 6.4: Unordered bucket balance

Define:

```lean
BucketBalance.internalLineCount
```

as `selfHalf.card / 2`.  Lean first proves the algebraic version:

```lean
BucketBalance.unordered_balance_of_even_selfHalf
```

If `selfHalf.card` is even, then

```text
2 * internalLineCount + |crossHalf| = |fiber| * |moves|.
```

The new orbit theorem removes the separate evenness assumption for the common
case of fixed-point-free involutive move systems:

```lean
BucketBalance.unordered_balance_of_involutive_fixedPointFree
```

That theorem is the abstract core of the tiling quotient identity.

### Theorem 6.5: Boolean-cube nonzero masks

The final cube-level specialization is now formalized as:

```lean
BucketBalance.unordered_balance_boolCube_masks
```

Lean defines `BoolCube index := index -> Bool` and the mask action
`xorMask u x`.  It proves:

```lean
BucketBalance.xorMask_involutive
BucketBalance.xorMask_fixedPointFree_of_nonzero
```

Thus every finite Boolean cube quotient, together with any finite family of
nonzero xor masks, satisfies the unordered bucket balance.  For staircase
tilings, the index type is the tile set.

## 7. Interaction Between Good-Cut Height and Bucket Transport

The two pieces of formalization meet conceptually.

Good-cut height is a coordinate:

```text
g(b) = number of base-path cuts crossed by upward tiles.
```

Bucket balance is a transport law:

```text
2 * internal + escaping = row size.
```

When a quotient map sends tilings to merged tournament classes, a move in the
tiling cube either stays inside a quotient bucket or escapes it.  The
half-line balance tells us that every bucket row has a fixed total mass.  The
good-cut coordinate then gives a way to stratify that transport by height
change.

Computations recorded in the repository suggest that nonzero changes in
good-cut height are strongly tied to quotient-boundary crossings.  In tested
layers for `n <= 6`, ordered half-lines with nonzero `Delta g` always changed
the merged tournament class.  The `Delta g = 0` stratum is where silent
defects, self-lines, ribs, and sea geometry concentrate.

This paper does not promote that observation to a theorem.  It records the
formal infrastructure that makes such a theorem possible:

- good-cut height is now a Lean object with exact spectrum;
- its top bucket is exactly strong connectivity;
- bucket transport has an audited conservation law;
- unordered internal transport is explained by finite orbit parity.

## 8. Lean Audit Summary

The relevant modules are:

```text
TournamentH7.StaircaseTileModel
TournamentH7.GoodCuts
TournamentH7.StaircaseConnectivity
TournamentH7.BucketBalance
TournamentH7.Verify
```

The key audited theorem names include:

| Theme | Lean theorem |
|---|---|
| No bucket 1 | `StTiling.goodCutCount_ne_one` |
| Bucket 0 | `StTiling.goodCutCount_eq_zero_iff_all_down` |
| Spectrum | `StTiling.goodCutCount_spectrum` |
| Tiling to tournament | `StTiling.toTournament_hasBasePath` |
| Good cut dictionary | `StTiling.isGoodCut_iff_crossesUpward_toTournament` |
| Top bucket equals SC | `StTiling.goodCutCount_eq_top_iff_toTournament_stronglyConnected` |
| Oriented row law | `BucketBalance.halfLine_balance` |
| Partner closure | `BucketBalance.pairHalf_mem_selfHalf` |
| No self-partner | `BucketBalance.pairHalf_ne_of_fixedPointFree` |
| Orbit parity | `BucketBalance.even_card_of_fixedPointFree_involutiveOn` |
| Internal evenness | `BucketBalance.selfHalf_card_even_of_involutive_fixedPointFree` |
| Unordered row law | `BucketBalance.unordered_balance_of_involutive_fixedPointFree` |

The `Verify.lean` audit reports the new good-cut, staircase-connectivity, and
bucket-balance theorems as depending only on Lean foundations:

```text
propext, Classical.choice, Quot.sound.
```

They do not depend on the project-specific axioms used by the H-value
forbidden-number theorems.

## 9. What Is Proved vs. What Remains Open

### Proved in Lean

The following are now formalized:

1. The good-cut set is the union of upward tile intervals.
2. No tiling has exactly one good cut.
3. For `n >= 3`, the exact abstract good-cut spectrum is
   `{0} union {2,...,n-1}`.
4. Concrete staircase tilings induce base-path tournaments.
5. Good cuts are exactly crossing-upward cuts of the induced tournament.
6. The top good-cut bucket is exactly strong connectivity.
7. Finite bucket half-lines satisfy oriented conservation.
8. Internal half-lines pair under fixed-point-free involutive moves.
9. Finite fixed-point-free involutions have even cardinality.
10. Fixed-point-free involutive move systems satisfy unordered bucket balance.
11. Nonzero Boolean xor masks satisfy the unordered bucket balance on finite
    Boolean cube quotients.

### Still open

The main formalization target listed in the previous draft is now closed at
the Boolean-cube level.  What remains is structural rather than bookkeeping:
attach the Boolean-cube theorem to concrete staircase tiling coordinates only
if a semantic wrapper proves useful, and study how quotient transport mass
decomposes across spine/ribs/sea and even-graph projections.

Other mathematical questions remain:

1. Does `goodCutCount` descend structurally to merged tournament classes, not
   merely computationally through the tested range?
2. Can the good-cut height and bucket transport laws predict the geometry of
   spine, ribs, and sea in `G_n/Z_2`?
3. Is there a closed asymptotic description of the connected-run covers in
   the good-cut interval gas?
4. Can bucket-balance row laws become cheap consistency checks for future TDA
   feature extractors?

## 10. Conceptual Summary

The formalized picture can be compressed into two slogans.

First:

```text
Good-cut height is interval geometry until the top bucket,
and the top bucket is strong connectivity.
```

The good-cut statistic begins as a one-dimensional union of intervals on the
cut path.  That alone explains the missing bucket 1 and the exact abstract
spectrum.  But once every cut is covered, the statistic becomes structural:
it is exactly the strong-connectivity predicate for the induced tournament.

Second:

```text
Unordered bucket transport is oriented conservation plus two-element orbits.
```

The row law counts every oriented half-line from a bucket.  Internal
half-lines come in partner pairs when moves are fixed-point-free involutions.
Escaping half-lines are seen one-sided from the source bucket and are not
halved.  The unordered balance is precisely this bookkeeping.

These two slogans explain why the recent Lean work matters.  It does not only
certify isolated lemmas.  It turns a collection of computations and geometric
pictures into a small calculus:

- interval covers for cut height;
- cut crossings for strong connectivity;
- half-line conservation for quotient transport;
- orbit parity for unordered internal mass.

That calculus is now available for future formalization and computation.

## Repository References

- `04-computation/lean/TournamentH7/TournamentH7/GoodCuts.lean`
- `04-computation/lean/TournamentH7/TournamentH7/StaircaseTileModel.lean`
- `04-computation/lean/TournamentH7/TournamentH7/StaircaseConnectivity.lean`
- `04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean`
- `04-computation/lean/TournamentH7/TournamentH7/Verify.lean`
- `01-canon/theorems/THM-336-good-cuts-structure.md`
- `01-canon/theorems/THM-346-tiling-quotient-bucket-balance.md`
- `01-canon/theorems/THM-348-finite-bucket-halfline-balance.md`
- `01-canon/theorems/THM-350-finite-unordered-bucket-balance-layer.md`
- `01-canon/theorems/THM-351-boolean-cube-mask-bucket-balance.md`
- `07-reflections/good-cut-spectrum-complete.md`
- `07-reflections/staircase-top-bucket-is-strong-connectivity.md`
- `07-reflections/unordered-bucket-balance-orbits.md`
- `07-reflections/boolean-cube-balance-as-checksum.md`
