# LRC(14): witness-hypergraph capacity after support compression

**Session status (2026-09-01):** the abstract deletion statements below are
**PROVED** for finite carrier systems; the numerical specialization is
**FINITE-EXACT** and promoted as
[THM-4309](../01-canon/theorems/THM-4309-lrc14-endpoint-595-support-threshold-residual-hypergraph-compression.md).
The recursive singleton-ideal sidecar is **FINITE-EXACT** and promoted as
[THM-4306](../01-canon/theorems/THM-4306-lrc14-index-265-recursive-ideal-two-mask-replacement.md).
LRC(14) remains **OPEN**.

## Inheritance pass

- Closest proved mechanism: THM-4302's inactive exchange preserves a
  9,019-mask carrier while repairing endpoint 596.
- Canonical hostile: THM-4303 leaves 145 pair-tagged bodies at
  `(96,595)`, `(100,595)`, and `(210,595)`.
- Corrected near miss: strict common inactivity has capacity only two, while
  the exact mixed response demand is ten.  It can therefore produce a
  9,027-mask repair, but cannot preserve size.
- Least-used sidecar: mask activity across the *whole preservation band*, not
  only disjointness at a current failed body.

The productive concept board was: pair-tagged response, common inactivity,
row-support size, residual witness edge, transversal dual, and recursive
singleton ideal.  The decisive change was to stop asking whether a mask was
individually deletable and ask which complete witness sets a simultaneous
deletion would erase.

## The finite carrier theorem

Let `R` be a finite carrier, `K` a finite row set, and `O` the body
obligations on those rows.  For `o=(p,B)`, define

```text
W_o={m in R : m is active at p and m is disjoint from B}.
```

Closure says exactly that every `W_o` is nonempty.  After deleting `D subset
R`, the surviving witness family is `W_o\D`; hence

```text
R\D closes K  iff  W_o is not a subset of D for every o in O.       (1)
```

This elementary equivalence has several useful consequences.

1. Safe deletion sets form a hereditary simplicial complex.  Its minimal
   nonfaces are the inclusion-minimal witness sets, not necessarily all raw
   `W_o`.
2. For a proposed pool `Q`, retaining `H subset Q` repairs deletion of `Q`
   exactly when `H` hits every witness edge contained in `Q`.
3. If `D_0` is already safe and `U_o=W_o\D_0`, only residual edges
   `U_o subset Q` constrain additional deletion from `Q`.  If their
   transversal number is `tau`, the exact additional deletion capacity is

   ```text
   |Q|-tau.                                                        (2)
   ```

4. The failure count

   ```text
   F(D)=sum_o 1[W_o subset D]
   ```

   is increasing and supermodular.  The slack for `A,B` counts precisely the
   witness sets split between the two exclusive sides.  Thus simultaneous
   deletion has increasing returns in the *bad* direction; singleton
   redundancy is intrinsically unable to predict a large deletion.

There is a dual set function on row sets.  If `I_m` is the inactivity set of
mask `m`, then

```text
f_C(S)=sum_m 1[S subset I_m]
```

is decreasing and supermodular.  Common-inactive exchange is its extreme
zero-support boundary.  Equations (1)--(2) explain why that boundary can be
far from the true carrier-compression capacity.

## Exact endpoint-595 specialization

Append the exact ten-mask response cover to THM-4302's carrier, producing
`R` of size 9,029.  On the 391-row preservation target, protect the 421 joint
masks and let `D_350` contain every other inherited mask active on at most 350
rows.  Exact support reconstruction gives

```text
|D_350|=5,141, FNV=03921cf597ee9863.
```

Deleting all of it is a useful hostile: the raw 5.59-billion-case replay
exposes 84 bodies on 21 rows.  Their residual hypergraph inside `D_350` has
70 responding masks and 53 response types.  A 37-mask cover and a 37-obligation
integral packing prove `tau=37`.  Retaining that cover and deleting the other
5,104 masks gives

```text
carrier size 3,925, ranks 3,858/67, FNV=6fbd0bffcf0ed78b,
391 rows, 5,594,095,650 row-body tests, zero failures.
```

The exact minimum 3,925 is only within the fixed family
`R\(D_350\H)`.  The threshold-360 hostile leaves 391 bodies exposed on 36
rows before repair.  That failure is informative: raw deletion count is the
wrong optimization target; the relevant quantity is the net capacity
`|Q|-tau(Q)`.

## Orthogonal ideal signal

The separate singleton ideal `H_265` contains 367 rows.  Removing its unique
inactive joint coordinate exposes eight private bodies.  No common-active
rank-eight mask covers all eight, while two masks do; matching packing and
cover certificates prove exact replacement minimum two.  Rebuilding that
deck closes all 367 rows.  Only its typed row consequence is unioned with
THM-4305: the decks themselves remain separate.

This suggests two complementary recursions:

- a **horizontal carrier recursion**, enlarging the row band and optimizing
  `|Q|-tau(Q)` over candidate deletion pools;
- a **vertical signature recursion**, replacing one inactive joint coordinate
  over a whole ideal and adding only the newly certified row set.

They share a response-hypergraph language but not a common carrier.  Treating
that shared syntax as permission to merge decks would lose the proof object.

## Next decisive probes

1. Directly replay the 3,925-mask carrier on the 22 residual endpoint-594
   rows.  Freeze the complete pair-tagged failure ledger before optimizing.
2. For several nested deletion pools `Q`, compute the exact or bounded ratio
   between `|Q|` and `tau(Q)`.  Support threshold should be a proposal
   generator, not an assumed optimum.
3. Compare support size with witness-edge degree, first blocker endpoint, and
   response-type rarity.  A useful score must predict net capacity on held-out
   thresholds and survive the threshold-360 hostile.
4. Recurse on the next largest singleton-signature ideals, but preserve the
   rule that each rebuilt deck is a separate proof node.
5. The main conceptual wall is still physical entry: a finite fixed-pool row
   certificate does not yet produce the arbitrary counterexample owner,
   phase, and arrival map required by LRC(14).
