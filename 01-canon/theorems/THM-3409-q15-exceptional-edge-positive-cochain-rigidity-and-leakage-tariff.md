---
id: THM-3409
title: "Q15 exceptional edge positive-cochain rigidity and leakage tariff"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED
source: root-2608-crouzeix-puzzle-2026-08-15
depends_on:
  - THM-3398-general-finite-mode-sheet-cover-cochain
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
verified_by:
  - 04-computation/lrc_q15_exceptional_positive_cochain_tariff_thm3409.py
  - 05-knowledge/results/lrc_q15_exceptional_positive_cochain_tariff_thm3409.out
---

# THM-3409 -- q15 exceptional edge positive-cochain rigidity and leakage tariff

## 1. Status, inheritance, and scope

This theorem exactly classifies one finite object: the owner edge

```text
q=15,                    U_*=(1,7,8,9,11,13).          (1)
```

The pinned full-clutter and mobile-centre audits identify `(1)` as the unique
rank-six physical edge on owners `1,...,14` which is absent from the capped
zero-cochain locus.  That uniqueness is inherited finite-exact input.  The
new result classifies **every** physical realization of this edge and proves
that the previously displayed cochain norms are exact minima.

The inheritance pass is:

- closest proved mechanism: [THM-3398](THM-3398-general-finite-mode-sheet-cover-cochain.md),
  which identifies physical covers with compatible selected modes and an
  exact complete affine cochain;
- canonical hostile: `(1)`, whose old single witness had nonzero cochain but
  did not prove a lower tariff;
- corrected near miss: the unnumbered mobile-centre audit called the old
  values witness norms, not minima;
- least-used sidecar: the signed pair cochain itself, compressed either on all
  pairs or on an optimal spanning tree.

This is not a uniform theorem over all physical edges, does not alter a live
LRC row, and gives no LRC(14) decrement.

## 2. Packet and tariff definitions

For `u in U_*` and `t in R/Z`, write

```text
D_u(t)={ell in Z/15Z: ||u(t+ell/15)||<1/14}.          (2)
```

A packet at `t` chooses one THM-3398 selected-block mode for every owner such
that its open mode interval contains `t` and the six selected blocks cover
`Z/15Z`.  If `x_i` is the containing lift of the chosen mode centre for owner
`u_i`, put, in the fixed increasing owner order `(1)`,

```text
p_ij=30 u_i u_j (x_i-x_j),             0<=i<j<6.     (3)
```

The lift is unique: every mode radius is at most `1/(14u_i)`, strictly less
than half the centre-lattice step `1/u_i`.  The values `(3)` are integers,
lie in the THM-3398 pair fibres, and obey

```text
u_k p_ij-u_j p_ik+u_i p_jk=0           (i<j<k).       (4)
```

For this fixed edge define the complete-pair leakage tariffs

```text
lambda_1(U_*)   = min_packets sum_(i<j)|p_ij|,
lambda_infty(U_*)=min_packets max_(i<j)|p_ij|,
lambda_2^2(U_*) = min_packets sum_(i<j)p_ij^2.         (5)
```

Because `(4)` lets any spanning tree reconstruct the complete cochain, also
define

```text
tau_1(U_*)    = min_(packet,tree T) sum_(ij in T)|p_ij|,
tau_infty(U_*)= min_(packet,tree T) max_(ij in T)|p_ij|. (6)
```

The tree tariff is the cheaper faithful sidecar; the complete tariff prices
the fully expanded pair ledger.

## 3. Exact classification

The physical cover locus of `(1)` in one time period is the disjoint union of
exactly `30` strict-open intervals.  Every interval has length

```text
1/3696,                                                   (7)
```

so the total physical cover measure is

```text
30/3696=5/616.                                           (8)
```

No endpoint is a cover.

On each interval there is exactly one packet.  Its six selected blocks equal
the six maximal actual danger sets `(2)` and form a partition of all fifteen
sheets.  Thus there are exactly `30` packets, not merely `30` sampled
witnesses.

Common sheet translation acts freely and splits them into two orbits of size
fifteen.  Canonical representatives are

```text
B_+=((0,1),(5,7),(8,10,12),(4,9,14),(2,6,13),(3,11)),

p_+=(2,2,3,3,4, -2,3,-1,2, 6,2,6, -6,-3, 5),          (9)
```

and

```text
B_-=((0,1),(9,11),(4,6,8),(2,7,12),(3,10,14),(5,13)),

p_-=-p_+.                                                (10)
```

The coordinates of `p_+` and `p_-` use lexicographic pairs

```text
01,02,03,04,05,12,13,14,15,23,24,25,34,35,45.          (11)
```

Each orbit has time measure `5/1232`.  Physical time reversal acts by
`t -> -t`, `ell -> -ell`, and `p -> -p` on every one of the thirty packets.
After canonicalizing the reversed blocks by the common sheet translation
`ell -> ell+1`, the displayed representatives `(9)` and `(10)` are related
by `ell -> 1-ell`.  This is an order-sensitive two-state residue invisible
to the zero/nonzero question alone.

## 4. Exact leakage tariffs

Every one of the thirty packets has

```text
sum |p_ij|=50,          max |p_ij|=6,
sum p_ij^2=206.                                         (12)
```

Consequently

```text
(lambda_1,lambda_infty,lambda_2^2)=(50,6,206).          (13)
```

All fifteen pair entries are nonzero.  In particular no realization of
`U_*` has zero complete cochain, independently recovering its exclusion from
the mobile common-centre locus.

There are `6^(6-2)=1296` labelled spanning trees.  Exact enumeration of their
restrictions gives

```text
(tau_1,tau_infty)=(10,3).                               (14)
```

Fifteen trees attain the first value and `104` attain the second.  One tree
attaining both is the star at owner seven.  In lexicographic edge orientation
its edges and signed values are

```text
((0,1),(1,2),(1,3),(1,4),(1,5)),
(2,-2,3,-1,2).                                         (15)
```

Thus a faithful tree compiler reduces complete absolute cost `50` to `10`
and the maximum edge scale from `6` to `3`, but it cannot erase the drift.

## 5. Why the finite proof is exhaustive

For each owner and sheet, predicate `(2)` changes only at

```text
t=k/u-ell/15 +/- 1/(14u).                              (16)
```

These are finitely many rational events modulo one.  Their common exact scale
for `(1)` is `15135120`; after deduplication there are `1260` boundary points
and `1260` complementary open cells.  A strict danger pattern is constant on
each open cell, so one midpoint represents the whole cell.  Boundary points
must be checked separately because `(2)` is strict.

Membership of `t` in a selected-mode interval is the conjunction that every
sheet in its selected block satisfies `(2)`.  Hence it changes only at the
same event set.  At each of the `2520` representatives the companion:

1. reconstructs every literal danger set directly from `(2)`;
2. enumerates every selected mode and independently computes its unique
   containing centre lift;
3. requires direct selected-block firing to agree with mode-interval
   membership;
4. tests all products of the six active mode banks for a sheet cover;
5. checks pair fibres, every triangle identity `(4)`, generalized-CRT
   realization, maximality, partition, and strict endpoints.

This makes the event-cell census a proof of the fixed finite statement, not a
random or bounded numerical probe.  The companion performs `579600` direct
mode-interval comparisons and obtains precisely the packets `(9)--(10)` and
their translations.

## 6. Information ledger and nonclaims

The source-to-target connection is:

| field | content |
|---|---|
| source | the exact physical cover locus of `U_*` |
| target | two signed cochain orbits plus full/tree tariffs |
| map | event cell -> maximal selected blocks -> unique centre lifts -> `(3)` |
| preserved | physical realization, sheet partition, strictness, translation orbit, time reversal, and complete compatibility |
| destroyed by norms | signs, pair locations, block order, and which of the two orbits occurs |
| required sidecar | retain `p_+/-p_+` or at least a signed reconstructing tree when orientation matters |

The scalar norms in `(13)--(14)` are tariffs, not classifiers: both time-
reversed orbits have identical norms.  Conversely, the orbit sign is not an
LRC certificate; one still needs a legal row embedding, owner chronology,
wall weights, and the physical safe-mass ledger.  No tournament is introduced:
the native carrier is a six-block partition with a weighted complete
cochain, and forcing its signed entries into pair orientations would discard
zero/scale and triangle data.

## 7. Reproduction

Run

```text
python3 04-computation/lrc_q15_exceptional_positive_cochain_tariff_thm3409.py
python3 -O 04-computation/lrc_q15_exceptional_positive_cochain_tariff_thm3409.py
```

The standard-library companion contains no floating-point literal and no
optimization-sensitive `assert`.  It imports the pinned THM-3398 primitives,
while the full-clutter and mobile-boundary artifacts are pinned as inherited
uniqueness inputs.

Artifact hashes after LF normalization:

```text
script   1d8528e187fa2804a61e3632dac4bb546473dd7c1e93749f6236f9ccb30f8243
output   ef1a0fa0f17411046f77c3c90fbe71009c6cddb14a91ae1df647e1a8ff2e1b20
semantic 80ecee351d7648a6a4c00f879a649d8fd4dbe87119929e5a5ded72d693c9859d
```

An independent no-import `Fraction` event census reconstructed every event,
packet, reversal pair, cochain, tariff, and hash, and normal/optimized replays
both byte-match the stored output.
