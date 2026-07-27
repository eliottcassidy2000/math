---
id: THM-2609
title: "External target-section itinerary saturation and root-state no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On THM-2600's canonical typed row, let Q_(s,ell) be the union, over its
  two positive constant-six middle rails, of absolute target sections whose
  globally primitive septimal Bockstein is a unit.  The 84 sets have sizes
  13 on 74 cells, 12 on 8 cells, and 11 on 2 cells; every difference set is
  all of F13.  Hence any external correction -7a can be placed between two
  unit sections in one repeated clock cell, while q=0 supplies the six
  intermediate clocks.  Each positive exact slice contains a base-13
  cylinder of depth 20, so the eight chosen events concatenate on one orbit
  in a cylinder of measure 13^-160.  This exhausts positivity, section
  differences, and finite temporal itinerary freedom, but does not create
  the missing state map: q is an external target-shift label and every
  physical arrival/future digit is h=6.  The common-event diagonal q=h is
  only (6,6), with zero correction.  Moreover, the 10 deficient section
  fibres have trivial translation stabilizer and cannot be the free C13
  ancestry torsor erased by positive time.  No transition kernel, semantic
  endpoint, row exclusion, or LRC(14) conclusion is proved.
source: carry-transition-cell-2026-07-28-target-section-itinerary
depends_on:
  - THM-2600-constant-six-middle-rail-common-x-atlas-and-uniform-bockstein-section
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
related:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2599-rootwise-opposite-shift-paired-slice-law
  - THM-2604-unshifted-future-root-accessibility-and-selector-cross-mixing-boundary
  - THM-2607-constant-six-rail-boundary-holonomy-invoice
  - THM-2610-chronological-paired-slice-marked-triangle-graft-and-action-axis-boundary
script: 04-computation/lrc14_external_target_itinerary_thm2609.py
output: 05-knowledge/results/lrc14_external_target_itinerary_thm2609.out
script_sha256: af46d7dd913d15d00c56dce45cc54d02e70b3bc2a815420eb3c6a19f5a595f74
output_sha256: 683c536e65ad2e7ab938085c740aee088436daf6570d93b07de4832a7785a49f
hash_basis: LF-normalized bytes
---

# THM-2609 -- external sections saturate; the event state does not

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2600 supplies a positive unit-Bockstein vertex in every one of its 84
base cells.  Keeping all thirteen absolute target sections reveals much more:
the available sections already have every additive difference.  Base-13
cylinder concatenation then realizes every prescribed correction itinerary
whose event labels lie in their corresponding available section sets, on one
physical orbit.

That is the strongest possible conclusion from these coordinates.  The
section label `q` is an external intervention in the target factor, while the
root digit `h` is a physical state of the same `x`.  The two coordinates do
not become one because their values happen to lie in the same field.  This
theorem proves both the saturation and the exact diagonal obstruction.

## 1. The unit-section sets

Use THM-2600's canonical row, response clock, and base grid

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5),

R=13^6=4,826,809,
T=297,836,897,838,480.                                  (1)
```

For each base cell

```text
(s,ell) in F_13^* x F_7,                                (2)
```

THM-2600 retains every positive middle rail

```text
(v,t)=(6,12) or (6,0).                                  (3)
```

For one such rail and one absolute target section `q`, its primitive owner
polynomial is the element `Y(z)` of

```text
R_7=F_13[z]/(Phi_7(z))                                  (4)
```

defined in THM-2600 (13), after the single global content

```text
g0=4,244,240                                             (5)
```

has been removed.  Define

```text
Q_(s,ell)
 ={q in F_13 : at least one positive rail (3)
                 has Y_q(z) a unit of R_7}.              (6)
```

This is a union over physical rails before any endpoint pair is selected.
It is not an assertion that one rail works for every `q`.

The exact reconstruction gives

```text
|Q_(s,ell)|=13 on 74 cells,
             12 on  8 cells,
             11 on  2 cells.                            (7)
```

In particular

```text
min_(s,ell)|Q_(s,ell)|=11,

0,6 in Q_(s,ell) for all 84 cells.                      (8)
```

The full seven-pattern census is frozen by the companion, together with a
deterministic endpoint-witness digest.

## 2. Every external correction occurs in every cell

For a nonempty set `A` in the prime cyclic group `F_13`, Cauchy--Davenport
gives

```text
|A-A| >= min(13,2|A|-1).                                (9)
```

Applying (9) to (7) gives the exact consequence

```text
Q_(s,ell)-Q_(s,ell)=F_13

for every (s,ell) in F_13^* x F_7.                      (10)
```

The companion also checks all 84 difference sets directly.  Thus for every
nonzero marker `a`, displacement `s`, and starting clock `ell0`, there are
unit sections `q0,q7` in the repeated cell `Q_(s,ell0)` such that

```text
q7-q0=-7a mod 13.                                       (11)
```

Equation (11) is an **external-label correction**.  It is numerically the
invoice appearing in THM-2542 and THM-2607, but no rail/deck or target/root
identification is inferred from that numerical agreement.

## 3. One-orbit realization of every prescribed correction itinerary

There is no finite-word recurrence obstruction left.  Set

```text
ell_i=ell0+i mod 7,

(q_i)_(i=0)^7=(q0,0,0,0,0,0,0,q7).                    (12)
```

By (8), every intermediate label belongs to `Q_(s,ell_i)`; by construction,
the two endpoint labels belong to the repeated starting cell.  This gives

```text
12 displacements * 7 starting clocks * 12 nonzero markers
 =1,008 prescribed eight-event unit-section itineraries. (13)
```

For each of the eight memberships in (12), choose one witnessing positive
rail from the union in (6), and then one positive `(ell5,r)` entry of that
unit aggregate slice.

These abstract labels can be realized on one physical orbit with an explicit
uniform clock gap.  Here is the elementary grid lemma used for that step.

**Grid-cylinder lemma.**  If a nonempty open interval has endpoints on the
`1/D` grid, then it has length at least `1/D`.  Whenever `13^d>2D`, it contains
the closure of a depth-`d` base-13 cylinder.

Indeed, take the first depth-`d` grid point strictly to the right of the left
endpoint.  It lies less than one cylinder length away, and the following grid
point is still strictly left of the right endpoint.

Every positive entry in THM-2600 (8) is a nonnegative finite sum of literal
Boolean intersections.  All present factors, rail factors, and deep probes
are resolved on the `1/T` grid, and pulling the delayed word back by `R`
refines it to the `1/(RT)` grid.  A positive entry therefore has a positive
Boolean component containing an open `1/(RT)`-grid atom.  Exactly

```text
D=RT=1,437,601,819,018,855,810,320,

13^19 < 2D < 13^20.                                    (14)
```

So depth `20` is the first depth certified by the generic grid lemma.  Choose
one depth-20 cylinder inside a positive component of each of the eight slices
in (12), and concatenate their 20 base-13 digits.  For the expanding map

```text
tau(x)=13x mod 1                                        (15)
```

the resulting depth-160 cylinder `C` satisfies

```text
tau^(20i)(C) lies in the chosen (s,ell_i,q_i) slice,
                                                   i=0,...,7,

mu(C)=13^-160.                                          (16)
```

This is an honest common-orbit itinerary.  The positive Boolean component at
each time belongs to a slice whose **aggregate** Bockstein is a unit.  Unitness
does not descend to a hidden sheetwise charge, and (16) does not claim that it
does.

## 4. The exact common-event obstruction

The physical arrival and delayed-word digit in every THM-2600 slice is

```text
v=h=6.                                                   (17)
```

The target section `q`, however, enters by replacing the present target shift
with `s5=-q`.  It is an externally chosen position, not the value of the
future-root digit.  Consequently the state support visible in one base cell is

```text
S_(s,ell)={(q,6):q in Q_(s,ell)}.                        (18)
```

Its common-event diagonal is exactly

```text
S_(s,ell) intersect {(q,h):q=h}={(6,6)}                 (19)
```

for all 84 cells.  Two diagonal endpoint states therefore have only the
correction

```text
6-6=0.                                                   (20)
```

Thus (10)--(16) cannot pay a nonzero marker after imposing the missing
same-event `q=h` state identification.  This is not a theorem that no such
intertwiner exists.  It proves that positivity, difference-set saturation,
and arbitrary finite chronology do not supply it.

### 4.1 The missing ancestry torsor

The obstruction also has an exact covering-space form.  For every `L>=1`,

```text
tau^L(x+r/13)=tau^L(x),                 r in F_13.        (21)
```

Translation by `r/13` acts freely on each physical inverse fibre.  Therefore
any lift that claims to restore the source `C_13` ancestry **equivariantly**
must adjoin a complete free 13-state branch coordinate; after choosing one
section, the other twelve are its deck translates.

The external section sets do not form such a bundle over all 84 cells.  Their
translation stabilizers have the exact census

```text
|Stab_F13(Q_(s,ell))|=13 on 74 cells,
                       1 on 10 cells.                    (22)
```

For the ten deficient fibres, this also follows abstractly: every nonzero
translation generates the prime cyclic group, so an invariant nonempty subset
would have to be all of `F_13`.  On the 74 full fibres, the free action in
`q` is still only an action on external labels.  Equation (18) keeps the
physical state fixed at `h=6`, and no map identifies `q`-translation with the
physical deck translation in (21).

Thus the itinerary in (16) is a path through chosen sections, not a lift of
the ancestry covering.  It chooses one cylinder branch at each time but does
not supply the common-endpoint 13-state fibre needed to compare branches.

THM-2602 explains the same boundary operationally.  The cylinder in (16) is
an intersection of independently chosen one-time gates.  Those gates belong
to the commutative diagonal vertex calculus.  They do not define a positive
ordered kernel

```text
K_ell(q,q')                                              (23)
```

on one common ancestry carrier.  Changing or deleting one time gate leaves
the other event labels well defined; no adjacent target index is transported.

## 5. Why later physical roots do not identify the states

THM-2610 carries one fixed marked source coefficient through every nonzero
future shift colour on a later same-root THM-2599 chamber.  It is the sharp
coefficient-preserving version of the THM-2599/2604 chronological graft.
Even there, the construction preserves each action at its own event:

```text
external target section q at one time,
physical future root k at another time.                 (24)
```

Even if `q=k` numerically, (22) is a temporal cospan, not a common-event state
map.  The future root remains independent of the relation/section residue.
THM-2607's `C_91` mapping-cylinder normal form reaches the same boundary from
the deck side: multiplication by `13` kills the `F_13` fibre that carried the
rail correction before chronology can identify it physically.

The needed sidecar is now precise: a lawful common-event correspondence from
the absolute target intervention `q` to the physical root state `h`, or an
ordered positive kernel retaining both indices and the semantic owner/repair
data.  More section support, more endpoint differences, or longer finite
itineraries cannot substitute for that sidecar.

## 6. Exact scope and evidence

The gain and loss ledger is:

| item | exact status |
|---|---|
| canonical THM-2600 row | retained |
| all 84 displacement/clock cells | retained |
| unit absolute target sections | exact census (7) |
| every external difference | proved in every cell |
| one physical eight-event orbit | depth-160 cylinder, measure `13^-160` |
| unit hidden sheet | not selected |
| same-event `q=h` | only `(6,6)` |
| `C_13` ancestry torsor | absent on 10 deficient fibres; untyped on 74 full fibres |
| adjacent positive transition | not produced |
| semantic owner/repair endpoint | not produced |
| row/LRC consequence | none |

Run

```text
python 04-computation/lrc14_external_target_itinerary_thm2609.py
python -O 04-computation/lrc14_external_target_itinerary_thm2609.py
```

The companion directly rebuilds the 162 middle rails, all `61,248` positive
fine entries, the single global primitive content, and all `1,740` unit
slices.  It forms the 84 unions (6), checks the size and pattern census, checks
all differences, translation stabilizers, and all 1,008 label itineraries,
freezes the endpoint-witness digest, verifies the cylinder invoice, and checks
the common-event diagonal.
All theorem checks raise under optimized Python.

Normal and optimized executions byte-match the stored transcript after LF
normalization.

Two independent immutable hostile audits replayed both executions from the
candidate commit, confirmed the declared LF hashes, rederived the complete
84-cell section census and difference-set argument, checked all 1,008
itineraries and the depth-160 cylinder invoice, and separately audited the
aggregate-versus-sheet unit distinction, the external-`q`/constant-`h=6`
typing, and the deficient-fibre stabilizer obstruction.  The second audit also
repaired the universal itinerary wording to its exact admissible-label
quantifier.  No theorem defect remains.

No row is removed and LRC(14) remains open.

QED.
