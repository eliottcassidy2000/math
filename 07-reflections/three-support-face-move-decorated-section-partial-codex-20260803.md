# Three support faces after restoring the physical successor

**Status: NON-CANONICAL FINITE-EXACT PARTIAL SCOUT.**  This reflection records
a bounded extension of
[THM-3286](../01-canon/theorems/THM-3286-three-face-availability-helly-defect-and-binary-origin-width.md)
on its inherited `94` triple-decision states.  Exact evidence is the
[move-decorated scout](../04-computation/fc3_three_support_face_move_decorated_section_partial_scout_20260803.py)
and its [stored output](../05-knowledge/results/fc3_three_support_face_move_decorated_section_partial_scout_20260803.out).
No theorem ID is claimed, and no `FC(3)`, `SFC(3)`, GMC, positivity, or
full-domain controller conclusion follows.

## Inheritance pass and concept board

The closest proved mechanism is THM-3286's static row-availability nerve.
Its canonical hostile states are `(3,4,5)` and `(1,3,4,5)`, inherited from
THM-3275, together with the proper three-way row defect `(5)`.  The corrected
near miss is THM-3286's binary sufficiency: a common row is enough only after
the physical successor has been forgotten.  The least-used relevant sidecar
is precisely that successor.

The live concept board was:

1. the projected row nerve `A_f(n)`;
2. the move-decorated fibre `B_f(n)`;
3. reset-direction compatibility;
4. fixed supplied-origin colourings; and
5. the finite internal transition graph and its reset-distance potential.

The applicable method cards were **Controlled forgetting and unlabeled
quotients require a sidecar** and **Audit sections, not only fibres**.  They
led directly to restoring the successor before asking for a common section.

## Typed map and loss ledger

For each face `f` and inherited triple state `n`, define

```text
B_f(n)={(r,n'):
  n -> n' is an actual nonempty one-pole physical move,
  the move lowers l1 distance to the reset of f by exactly one, and
  response row r strictly increases from n to n'}.
```

The projection

```text
(r,n') |--> r
```

recovers `A_f(n)` exactly on all `94` states.  Its preserved predicate is
existence of some reset-directed strict ascent for row `r`.  It destroys
which successor realizes the ascent and therefore whether two faces can use
the same physical move.  The restored sidecar is the literal successor
multiset.  The cheapest hostile test is exact intersection of the decorated
sets, not comparison of their row and successor projections separately.

The scout rebuilds the three pole banks and resets from the pinned THM-3238
product-Gamma source.  It computes the response coordinates on the exact
`943`-state closure consisting of the `94` sources and all their facewise
successors.  Projection reproduces THM-3286's frozen triple-atlas digest.

## Decorated defects

Restricting every pair census to the same `94` triple states, the decorated
empty-intersection counts are

| face pair | empty decorated fibres | empty row projections | incompatible-successor fibres |
|---|---:|---:|---:|
| `F12,F13` | 13 | 2 | 2 |
| `F12,F23` | 5 | 0 | 3 |
| `F13,F23` | 0 | 0 | 0 |

Thus `F13` and `F23` retain a common decorated move everywhere, while both
other pair relations acquire new defects.  The triple decorated intersection
is empty at `16` states rather than THM-3286's three projected states.  Exactly
two are proper decorated Helly defects, with every decorated pair still
nonempty:

```text
(5), (3,4).
```

The full loci and exact digests are frozen in the output rather than promoted
to canon.

## The supplied-origin width survives, but its orientation becomes rigid

Locally, `78` states need one origin value and the `16` decorated triple
defects need two.  No state needs three.  Globally, the minimum fixed alphabet
also remains binary, but decoration eliminates the colourings in which `F23`
shares the `F12` class.  In face order `(F12,F13,F23)`, the complete optimal
list is

```text
(0,1,1), (1,0,0).
```

Hence there is one colouring up to complement:

```text
F12 | (F13,F23).
```

There is no uncoloured global section, but a binary-coloured static section
does survive.  This is sharper than either slogan “the move decoration kills
the section” or “THM-3286's four colourings persist.”

## Finite iteration only

For the canonical colouring `(0,1,1)`, the scout checks every common
decorated edge, not only one chosen policy.  Every edge internal to the
`94`-state graph lowers the reset distance of each face in its colour block.
The complete allowed internal graph is acyclic.  Its longest internal path is
`5` edges in the singleton `F12` block and `6` edges in the shared
`F13,F23` block.  The lexicographically least sections exit the inherited
graph in at most `6` moves in either block.

Exit is deliberately treated as termination of this bounded test.  The scout
does not continue into the full face banks, compose across a changed active
face set, infer a chronological controller, or count adjacency powers as
physical time.

## Failure anatomy and controls

The first new failure, ordered by state size and then face-pair order, occurs
at `n=(4,5)` for `(F12,F13)`.  Both projected sets contain row `5`, and the
faces have common reset-directed successors `(1,4,5)` and `(3,4,5)`.  But row
`5` ascends on `F12` toward those successors while on `F13` it ascends only
toward `(4,5,6)`.  Therefore the separate common-row and common-successor
tests pass while the common pair test fails.  This isolates the first failed
implication:

```text
common projected row + common physical successor
does not imply a common (row,physical successor).
```

The minimal pure reset-incompatibility hostile is `(2,3,4,5)` for
`(F12,F23)`: projected rows `5` and `12` are common, but the two faces have no
common directed successor.

At control state `(5)`, the three faces have common possible successors
`(3,5)` and `(4,5)` but no common decorated pair.  The inherited pair
`F12,F13` still has the single witness `(12,(1,5))`; the other two decorated
pair intersections remain nonempty.  Thus `(5)` remains a proper Helly
defect after decoration.  The inherited row conflicts `(3,4,5)` and
`(1,3,4,5)` remain empty for `F12,F13` before and after decoration.

There is one important restriction boundary.  The canonical row alphabet
`{16,17}` supports the surviving binary colouring on these `94` states.
This does not repair THM-3286's proved hostile: `{16,17}` has no extension on
all `8,394` face/state decisions.  The full-domain failure is pinned rather
than redundantly recomputed here.

## Frontier and stopping reason

The exact positive signal is a move-compatible binary section with a strict
finite potential.  The stopping reason is equally exact: the state universe
is not closed under every selected successor, and exiting it changes the
typing needed for iteration.  A next scout would have to enlarge the carrier
to all relevant face/state obligations and retain active-face changes and
history; it must not infer those data from this terminating restriction.

The THM-3300/3301 invariant Gaussian-moment lane is analytically orthogonal
and is not a dependency or consequence here.  Any future bridge would need an
explicit map from these physical response moves to its moment predicate.

Reproduce with

```text
python3 04-computation/fc3_three_support_face_move_decorated_section_partial_scout_20260803.py
python3 -O 04-computation/fc3_three_support_face_move_decorated_section_partial_scout_20260803.py
```

and compare both transcripts with the stored output.
