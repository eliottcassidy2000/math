# Three support faces: a proper Helly defect without a ternary origin tax

**Status: NON-CANONICAL SYNTHESIS AROUND PROVED, VERIFIED-EXACT,
INDEPENDENTLY HOSTILE-AUDITED THM-3286.**  The canonical result is
[THM-3286](../01-canon/theorems/THM-3286-three-face-availability-helly-defect-and-binary-origin-width.md).
Its primary executable evidence is
[`fc3_three_support_face_availability_hypergraph_scout_20260803.py`](../04-computation/fc3_three_support_face_availability_hypergraph_scout_20260803.py)
with its frozen
[`output`](../05-knowledge/results/fc3_three_support_face_availability_hypergraph_scout_20260803.out),
and the separate hostile reconstruction is
[`fc3_three_support_face_availability_hypergraph_independent_audit_20260803.py`](../04-computation/fc3_three_support_face_availability_hypergraph_independent_audit_20260803.py)
with its frozen
[`output`](../05-knowledge/results/fc3_three_support_face_availability_hypergraph_independent_audit_20260803.out).

## Inheritance and the changed object

The closest proved mechanism is
[THM-3275](../01-canon/theorems/THM-3275-unrestricted-twenty-two-row-face-blind-selector-obstruction.md):
on the support-`(1,2)` and support-`(1,3)` bank-`I2` faces, exactly the raw
states `(3,4,5)` and `(1,3,4,5)` have disjoint all-row availability sets.
[THM-3278](../01-canon/theorems/THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary.md)
then identifies their restriction to the twelve-row core with its unique
bipartition.  Its conclusion is still explicitly two-face.

The canonical hostile pair is therefore the adjacent two-state seam above.
The corrected near miss is more revealing: singleton `(5)` was a positive
control in THM-3275 because the first two faces share row `12` there.  After
adding the natural third support face `(2,3)`, `(5)` remains compatible on
every pair but becomes incompatible on the triple.  Thus the right object is
not another pair graph or tournament.  It is the Čech nerve of row-witness
sets over raw physical states, as suggested by the underused nerve sidecar in
[THM-2292](../01-canon/theorems/THM-2292-common-catalytic-section-and-helly-calibration-nerve.md)
and the `Audit sections, not only fibres` method card.  THM-2292's convex
calibration theorem is not being imported: these row sets are finite and
discrete, and the exact three-way failure is computed directly.

## Typed map and loss ledger

For each named face `f` and nonreset physical state `n`, send

```text
(f,n) |--> A_f(n) subset {1,...,22},
```

where `r` lies in `A_f(n)` precisely when response row `r` has a strict
one-pole ascent along a move decreasing multiset distance to that face's
reset.  Intersections of the `A_f(n)` preserve exactly the existence of one
row serving all listed face copies of the raw state.

This map forgets response magnitude, the chosen target neighbour, chronology,
and every Gaussian-moment or positivity implication.  The tested restoration
sidecar is a static face label.  The cheapest decisive test is exhaustive:
derive each pole bank and reset from the pinned product-Gamma formula,
construct the physical universes by two independent methods, rebuild every
response vector, and enumerate every intersection.

## Exact signal

The third face is support-`(2,3)` in bank `I2`.  Its pole multiplicity profile,
reset and bank size are

```text
((1,2),(2,4),(3,1),(4,3),(5,3),(6,1),(7,1),(8,1)),
(2,3,4,4,5,6,7,8),
3839 states (3838 nonreset decisions).
```

The profile differs from both promoted faces, so this is not a literal replay
of either state bank by renaming.  Direct reconstruction gives at least one
available row at every nonreset state on all three faces.  The three pair
loci are

```text
H(12,13)={(3,4,5),(1,3,4,5)}   on 238 joint decisions,
H(12,23)=empty                  on 94 joint decisions,
H(13,23)=empty                  on 1726 joint decisions.
```

There are `94` triple joint decisions.  Exactly three have empty total row
intersection:

```text
(5), (3,4,5), (1,3,4,5).
```

At `(5)` the three row sets are

```text
A_12={2,8,9,11,12,14,16,18,22},
A_13={3,4,5,6,7,10,12,13,17,19,20,21},
A_23={3,4,5,6,7,8,10,13,16,17,18,19,20,21,22}.
```

Their pairwise intersections are respectively

```text
{12}, {8,16,18,22},
{3,4,5,6,7,10,13,17,19,20,21},
```

but the triple intersection is empty.  This is the unique proper three-way
Helly defect.  The other two empty triples already contain the THM-3275 empty
pair and are therefore inherited rather than higher-order defects.

## The origin cost does not jump to three values

Grouping every raw state by its active faces, the exact histogram

```text
(active faces, minimum local labels):
(2,1):1776, (3,1):91, (3,2):3
```

shows that no state needs three local branches.  Globally, two fixed labels
are necessary and sufficient.  In face order `(F12,F13,F23)`, all successful
binary assignments are

```text
(0,1,0), (0,1,1), (1,0,0), (1,0,1).
```

Thus `F12` and `F13` must remain opposite, exactly as THM-3275 requires, while
`F23` may share either class.  The minimum static origin width remains one
bit, but its raw-vector dependency locus grows from two states to

```text
{(5),(3,4,5),(1,3,4,5)}.
```

This separates two invariants that pairwise analysis conflates: the nerve
has a genuine missing two-simplex, while the witness-colouring width is still
binary.

The hostile reconstruction also exposes a boundary hidden by the bare width:
THM-3278's canonical row pair `(16,17)` has no global three-face extension
under any of the four optimal binary face colourings.  One bit suffices only
when each labelled block may choose a row state by state from all twenty-two
lawful rows.  Binary origin width therefore does not imply a fixed binary
row alphabet.

## Proved theorem statement

On the three named bank-`I2` support faces `(1,2)`, `(1,3)`, and `(2,3)`, with
the twenty-two lawful response rows inherited from THM-3249, every nonreset
physical state has a strict reset-directed row.  The pair obstruction loci
are exactly those displayed above; the triple obstruction locus consists of
the same two inherited states plus singleton `(5)`, which is the unique proper
three-way Helly defect.  A fixed face-origin alphabet has minimum cardinality
two, with exactly the four binary assignments displayed above, and every
lawful selector must consult that bit on exactly the three obstruction
vectors.  This is a finite response-bank theorem only; it proves no
`FC(3)`, `SFC(3)`, Gaussian Moment Conjecture, or positive-functional claim.
This statement is now canonical as THM-3286; this reflection records the
object choice, loss ledger, and frontier interpretation rather than serving
as its proof surface.

## Frontier after the signal

The immediate positive lesson is that completing the support triangle does
not force a ternary controller.  The immediate hostile lesson is that
pairwise compatibility is not a gluing theorem: even THM-3275's singleton
positive control can become a higher-order obstruction.  The failure of the
fixed pair `(16,17)` suggests two separate next tests: search for the smallest
fixed row alphabet that supports a global three-face policy, and test whether
history-aware selectors can reduce that alphabet without losing chronology.
Other signed banks remain a second axis.  Any attempted
moment conclusion must separately restore response magnitude, positivity,
and a target-preserving chronological transition; the present nerve retains
none of them.
