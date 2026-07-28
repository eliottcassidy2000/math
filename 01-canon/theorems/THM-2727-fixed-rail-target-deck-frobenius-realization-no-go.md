---
id: THM-2727
title: "Fixed-rail target-deck Frobenius realization no-go"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED.  On the complete two-slice c2-safe bank, the maximal
  formal fixed-semantic coefficient extension of the THM-2716 target deck
  fixes rail, owner, support, and private edge, sends digit (0,1) to (12,0),
  shifts carry by +6, and applies z -> z^(-1) in F13[C7].  It does not pair
  the forward and reverse physical banks: both same-clock semantic-signature
  intersections are empty, including after a non-inherited root-preserving
  edge swap.  Exactly fourteen nonzero inherited matches survive; all are
  clock-constant units on rail 104.  This is a realization no-go, not a row
  exclusion, endpoint current, semantic-owner transport, or LRC(14) result.
source: root/bs13-fixed-rail-target-deck-2026-07-28
depends_on:
  - THM-2041-frobenius-stability-of-exact-period-projectors
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary
  - THM-2706-relative-segal-macro-cycle-and-minimal-ghost-midpoint-completion
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
related:
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
  - THM-2717-minimal-c2-safe-ghost-transit-rebuild-and-deletion-cancellation-boundary
script: 04-computation/lrc14_fixed_rail_target_deck_frobenius_no_go_thm2727.py
output: 05-knowledge/results/lrc14_fixed_rail_target_deck_frobenius_no_go_thm2727.out
script_sha256: a636819051de5b8b7045e1b5eb1ca694d1408b227e33dad17f17d236a70a91fc
output_sha256: f4cc6e406d8361a856aafab52fdba5fad8578904007da5d80825c5ae53f00504
hash_basis: LF-normalized bytes
---

# THM-2727 -- fixed-rail target-deck Frobenius realization no-go

## Status and scope

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  This theorem tests the cheapest possible physicalization
of the proved THM-2716 arm bitorsor on the complete two-slice `c2`-safe bank
which THM-2717 is separately packaging.
It defines the maximal formal fixed-semantic coefficient extension already
licensed by the delayed `C4` orbit, the private-root law, and Frobenius on
`F13[C7]`, while deliberately leaving rails and semantic supports fixed.  That
extension does not give a forward/reverse ghost-bank bijection.  A small
clock-constant coefficient subbank survives, but no semantic endpoint current
or row consequence follows.

The companion is self-contained relative to the immutable THM-2698 exact
carrier and pins that dependency's LF hash.  It does not import a THM-2717
implementation detail.  Canonical promotion is intentionally deferred until
THM-2717 itself is promoted and the theorem ID is rechecked against the live
remote namespace.

## 1. The lawful domain of the target deck

On the exact THM-2706 delayed four-cycle, the target deck is

```text
J=B^2,                    1/17 <-> 16/17.                 (1)
```

For a non-seam point with delayed half-digit

```text
d=2h+kappa in {0,...,25},                                (2)
```

reflection of the delayed phase sends

```text
d -> 25-d,
(h,kappa) -> (12-h,1-kappa).                             (3)
```

Thus the two forced slices are

```text
(h,kappa)=(0,1) -> (12,0).                               (4)
```

The proved THM-2640 private-root formula is

```text
r=2c+floor(d/13)+1_(edge=left) mod 13.                   (5)
```

Frobenius fixes the `F13` root colour.  Keeping the private edge fixed, the
unique carry relabelling which keeps `(5)` fixed under `(4)` is

```text
c -> c+6 mod 13,                                         (6)
```

because `2*6+1=0 mod 13`.  In particular the displayed forward carry `7`
maps to the displayed reverse carry `0` and both have root `2`.

For a normalized seven-clock profile

```text
Y(z) in F13[z]/(Phi_7),                                  (7)
```

absolute Frobenius is coefficient inversion:

```text
Y(z)^13=Y(z^13)=Y(z^(-1)),                               (8)
```

since `13=-1 mod 7`.  Equations `(3)`, `(6)`, and `(8)` therefore define the
maximal formal fixed-semantic coefficient extension on the already-typed data

```text
J_coeff:(j,edge,c,0,1,Y)
       ->(j,edge,c+6,12,0,Y(z^(-1))).                    (9)
```

The rail index `j`, its source/owner/theta metadata, the present-factor
support, and the private edge are fixed.  Strictly speaking, THM-2716 acts
only on the delayed `C4` incidence; `(8)` is the canonical algebraic sidecar,
and identity is the only inherited action on the remaining labels.  Acting
on rails would be an additional semantic intertwiner, exactly the datum under
test rather than an inherited symmetry.

The companion also tests the largest obvious private-orientation enlargement:
swap left/right edges while keeping the root fixed.  Equation `(5)` then
forces carry changes

```text
left c -> right c,             right c -> left (c-1).     (9a)
```

This edge involution is not supplied by the target deck; it is included as a
hostile against the claim that fixing the edge caused the failure.

Likewise, `(9)` does not define a map on physical addresses `N`: many lifts
have the required carry residue, and canon supplies no distinguished lift.

## 2. Full two-slice coefficient census

The exact probe rebuilds `(h,kappa)=(0,1),(12,0)` on all 162 THM-2717 rails,
using the global safe-lattice content `26`.  The results are

```text
forward nonzero profiles                  3388
reverse nonzero profiles                  3432
one-sided reverse profiles                  44
paired nonzero profiles                   3388

forward primitive units                   3018
reverse primitive units                   3012.           (10)
```

Frobenius preserves the multiplication determinant of every source profile,
as it must.  Realization by the opposite safe row is much rarer:

```text
exact Frobenius equalities                  366
equalities of the zero cyclotomic class     352
nonzero exact equalities                     14
nonzero equalities up to F13* scalar         54
nonzero equalities up to scalar/clock shift 314.          (11)
```

The last two lines are not same-clock functors: they choose a scalar and,
in the last line, a monomial clock origin profile by profile.

Under the extra edge swap `(9a)`, the corresponding counts are

```text
paired nonzero profiles                    3080
nonzero exact Frobenius equalities           39
nonzero scalar equalities                     50
nonzero scalar/clock-shift equalities         302.         (11a)
```

All thirty-nine exact classes are units.  This is a larger coefficient-only
relation, but it depends on an edge action absent from THM-2716.

All fourteen genuine equalities in `(11)` are primitive units on the single
rail

```text
j=104,        (source,owner,theta)=(8,6,12),              (12)
```

and every surviving class is clock-constant in the canonical quotient basis,

```text
(a,0,0,0,0,0),                    a!=0.                  (13)
```

There are seven allowed carries on each of the two private edges.  This is
the strongest inherited positive survivor: Frobenius is exactly realized on a
fourteen-profile constant subbank.  It is not a source-one rail, carries no
nontrivial `C7` clock information, and is disjoint from the forced ghost
semantics.

## 3. The forced coefficient mismatch

Using the global content `26`, the canonical reduced class of the forced
forward point is

```text
Y_f=(7,0,0,0,2,10).                                      (14)
```

Its lawful image is

```text
Fr_13(Y_f)=(7,0,10,2,0,0).                               (15)
```

The opposite digit slice on the **same** semantic rail `j=9` and the actual
reverse forced rail `j=2` give, respectively,

```text
Y_same=(0,0,0,0,8,8),
Y_reverse=(0,0,3,0,0,0).                                 (16)
```

The additional edge-swapped same-rail candidate is

```text
Y_edge=(9,0,0,0,6,6).                                    (16a)
```

Their canonical clock-support sizes are

```text
3, 2, 1.                                                  (17)
```

The support sizes in `(17)` immediately rule out exact equality and equality
up to a nonzero scalar.  Because canonical reduction modulo `Phi_7` can alter
support size after a monomial shift, the companion separately tests all seven
clock translations and all nonzero scalars; none maps `(15)` to either class
in `(16)`.  Thus support size is the cheapest mismatch, and the exhaustive
affine-clock check is the stronger invariant.
The same exhaustive check rejects `(16a)`, so changing the private edge does
not repair the forced coefficient.

## 4. Same-clock physical cell banks

The probe next scans the actual source-one rail/present/private-root cells at
both ghost phases while holding the semantic clock and rail set fixed.  It
finds

```text
clock 4, rails {8,9}:
  phase 1/17                         5850 cells
  phase 16/17                       5848 cells

clock 1, rails {2,3}:
  phase 1/17                         4958 cells
  phase 16/17                       4958 cells.            (18)
```

At clock four the cardinality defect `2` already forbids a bijection.  At
clock one the cardinalities agree, but after `(9)` the multisets of full
semantic signatures

```text
(j,edge,carry,root,normalized C7 class)                  (19)
```

have intersection zero.  The same is true at clock four.  Repeating the test
with the edge involution `(9a)` again gives intersection zero at both clocks.
Thus equal mass is not enough: even after forgetting the physical address,
no fixed-rail signature-preserving pairing exists, with or without private
edge reversal.

The displayed forward and reverse THM-2706 witnesses are weaker still as a
candidate pair: they use clocks `4` and `1` and rails `9` and `2`.  Mapping
one to the other would first require an action on owner clocks and rails,
which `(9)` lawfully refuses to invent.

## 5. Consequence and first false implication

The exact result is

```text
Frobenius preserves every abstract primitive C7 unit,
but generally does not realize that unit in the opposite c2-safe row;

the THM-2716 arm bitorsor has no same-semantic-support realization
inside the current THM-2717 safe transit bank.             (20)
```

The first false implication is therefore

```text
algebraic Frobenius/unit preservation
  ==> opposite-ghost safe-profile realization.           (21)
```

The fourteen constant profiles show that `(21)` can hold on a special
coefficient subbank, while `(17)` and `(18)` prove that it fails exactly
where the macro bridge needs it.

The strongest inherited architecture is a **partial coefficient relation**:
retain the fourteen constant matches and the general abstract Frobenius image
as formal target objects, but do not call either a physical transit functor.
If an independent private-edge involution is later constructed, the
thirty-nine matches in `(11a)` enlarge that formal relation.  To reach a
current one still needs a new rail/owner intertwiner plus a distinguished
address lift; an edge choice, scalar, or clock origin is insufficient.

## 6. Exact reproduction

Run normally and with optimization:

```text
python 04-computation/lrc14_fixed_rail_target_deck_frobenius_no_go_thm2727.py
python -O 04-computation/lrc14_fixed_rail_target_deck_frobenius_no_go_thm2727.py
```

The companion directly rebuilds the two safe digit slices on all 162 rails,
checks private-root support and global-content divisibility, verifies the
Frobenius determinant identity profile by profile, classifies exact/scalar/
affine-clock matches, proves `(12)--(17)`, and scans all four same-clock cell
banks in `(18)`.
