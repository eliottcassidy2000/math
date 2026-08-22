# The shared colour lane terminates at a Boolean reset fork

**Status: NON-CANONICAL FINITE-EXACT PARTIAL SCOUT.**  This note continues
the bounded move-decorated experiment in
[the 94-state scout](three-support-face-move-decorated-section-partial-codex-20260803.md).
The exact evidence is the
[all-choice closure script](../04-computation/fc3_three_support_face_move_decorated_all_choice_closure_scout_20260803.py)
and its
[stored output](../05-knowledge/results/fc3_three_support_face_move_decorated_all_choice_closure_scout_20260803.out).
No theorem ID is claimed.  Nothing here proves `FC(3)`, `SFC(3)`, GMC,
positivity, or an unbounded chronological controller.

## Inheritance pass and the changed question

The closest proved mechanism is
[THM-3286](../01-canon/theorems/THM-3286-three-face-availability-helly-defect-and-binary-origin-width.md):
all three named faces have a static binary row-availability colouring.  The
canonical hostile states remain `(3,4,5)` and `(1,3,4,5)`, together with the
proper Helly defect `(5)`.  The corrected near miss is the preceding scout's
treatment of an edge leaving the common 94-state carrier as termination.
That was a legitimate stopping rule for that scout, but not evidence of a
completed chronology.  The least-used sidecars were the exact reset of each
face and the active face set after a move.

The live concept board was therefore:

1. the surviving static colouring `F12 | (F13,F23)`;
2. literal `(row,successor)` fibres;
3. typed lane nodes rather than untyped raw states;
4. the sum of active reset distances;
5. the fixed row alphabet `{16,17}`; and
6. the first closed frontier on which the two reset directions separate.

The operative method cards were **Controlled forgetting and unlabeled
quotients require a sidecar** and **Audit sections, not only fibres**.  The
new question is not whether each row fibre is populated, but whether the
intersection of the two *directed successor cones* remains populated after
iteration.

## The smallest typed all-choice closure

A node is `(A,n)`, where `A` is the currently active subset of one original
colour block and `n` lies in every face decision universe named by `A`.  The
188 starting nodes are

```text
({F12},n) and ({F13,F23},n),     n in the 94-state triple carrier.
```

For a common decorated edge `(r,n')`, a face is deleted from `A` only when
`n'` is its exact pinned reset.  Every retained face must still type `n'` as a
decision.  A face is never silently acquired.  Closing under *every* such
edge gives the least fixed point of this rule.

This construction deliberately treats the two supplied origins as separate
lanes.  It is not a shared-state interleaving in which an `F12` move also
mutates the `F13,F23` state.  That missing coupling is outside the current
carrier.

All `8,397` face/state responses were recomputed on the complete inherited
face universes before taking the closure.  The result is

```text
typed nodes                         854
underlying untyped physical states 760
F12-typed nodes                      94
F13,F23-typed nodes                 760
unrestricted row-labelled edges 30,853
terminal row-labelled edges          69
blocked typed nodes                    4
```

The distinction `854 != 760` is load-bearing: the initial 94 physical states
occur once in each coloured lane.

## A strict potential and the complete maximal-path census

For a typed node define

```text
Phi(A,n) = sum_(f in A) ||n-reset_f||_1.
```

Every lawful edge lowers each active face's reset distance by one, hence

```text
Phi(A',n') = Phi(A,n)-|A|.                              (1)
```

This proves acyclicity without a graph search.  The only active-face changes
in the closure are `69` row-labelled `F12 -> empty` terminal edges.  In
particular, the shared lane never reaches either face reset and never changes
to a singleton lane.

For the 94 starts, every maximal `F12` path reaches its reset in `1..6`
moves, with histogram

```text
1:7, 2:20, 3:30, 4:25, 5:11, 6:1.
```

Every maximal shared-lane path instead reaches a blocker in `3..9` moves,
with histogram

```text
3:4, 4:15, 5:27, 6:25, 7:16, 8:6, 9:1.
```

There are no shared-lane terminal edges.  Thus the four blockers are not
bad choices among otherwise terminating paths: by finiteness and (1), every
maximal all-choice path in that lane ends at one of them.

## The Boolean reset fork

Put

```text
C=(3,4,5,6,7,8).
```

The exact closure blockers are

```text
C, C union {1}, C union {2}, C union {1,2}.             (2)
```

These are a Boolean square in the two low-root coordinates.  At the smallest
state `C`, the two facewise reset-directed successor sets are

```text
F13: add 1 or add another 3;
F23: add 2 or add another 4.                            (3)
```

They are disjoint.  Nevertheless both faces have the same ten projected
available rows there,

```text
{2,5,8,9,11,12,14,16,18,22}.                           (4)
```

The other three states in (2) likewise have nonempty projected row
intersection, always containing row `16`, and empty common-successor
intersection.  The first failed implication is therefore sharper than in
the 94-state scout:

```text
common row, even canonical row 16
does not imply a common reset-directed physical move.  (5)
```

Across all active contexts of this fixed colouring there are seven
unrestricted decorated blockers: the four in (2), plus three states on the
adjacent double-`4` layer.  The 94-state-generated closure reaches exactly
the four single-`4` states in (2).

The obvious two-pole macro at `C` was also checked.  Combining one desired
`F13` addition from `{1,3}` with one desired `F23` addition from `{2,4}`
leaves only the `(+1,+2)` and `(+1,+4)` targets typed on both faces; the two
macros using `+3` exceed the `F23` multiplicity of root `3`.  Both typed
targets have the ten common strict rows in (4), including row `16`, but each
preserves the distance vector `(2,2)` and lands at another unrestricted
blocker.  Moreover, neither macro has a common first one-pole step.  It is
therefore a response-positive, potential-neutral diagonal across the fork,
not a chronological bridge.

## What happens to `{16,17}`

Restricting every edge to rows `16` and `17` produces exactly the same 854
typed nodes and exactly the same four blockers.  It retains `3,126`
row-labelled edges and seven terminal `F12` edges.  A minimal path to the
smallest blocker is

```text
(3,4,5)
  -- row 17, add 6 --> (3,4,5,6)
  -- row 16, add 7 --> (3,4,5,6,7)
  -- row 16, add 8 --> (3,4,5,6,7,8)=C.                (6)
```

At `C`, unrestricted rows fail too.  Thus the canonical row pair does not
cause the first chronological obstruction on this closure; it reaches the
same universal reset fork as all twenty-two rows.

This does not weaken THM-3286's full-domain hostile.  Collapsing equal-label
face obligations gives `6,668` active coloured contexts across the inherited
face universes.  Seven have no unrestricted decorated move, while 38 have no
`{16,17}` move.  The present closure avoids the 31 row-specific failures and
reaches four of the seven universal ones.  THM-3286's statement concerns all
`8,394` separate face/state decisions and proves that no fixed canonical-pair
extension is globally lawful there.  It remains a pinned control, not a
claim contradicted or recomputed away here.

## A scheduling signal, not yet a bridge

There is one unexpectedly clean resource-ordering signal.  For each of the
94 common starts, compare the deterministic `F12` time to reset with the
deterministic shared-lane time to the potential-four fork.  Their difference

```text
pair-fork time - F12-reset time
```

is never negative, with histogram

```text
0:7, 1:24, 2:16, 3:16, 4:31.                           (7)
```

Equality occurs exactly on the seven nonempty subsets of `{3,4,5}`.  This
suggests a token-recycling reframe: the `F12` origin is always free no later
than the moment at which `F13` and `F23` must split.  But reusing it would
require a lawful origin-rebinding or state-copy map from the completed `F12`
lane to one branch of the shared state.  Neither operation exists in the
current typed system.  Equation (7) is therefore a precise scheduling
sidecar and a cheap next test, not a dynamic two-origin controller.

## Failure anatomy and next decisive test

The static binary colouring survives as an availability section, but its
fixed shared colour class does not continue to both resets.  The mechanism
is not a response-row shortage, cycle, or unlucky policy.  It is a terminal
fork of the reset direction cones after a strictly decreasing potential has
removed all shared directions.

The next bounded experiment should add exactly one new coordinate:

```text
origin occupancy + a typed rebind/copy permission at the four states (2).
```

It should test whether the pointwise slack in (7) can be converted into a
commuting square that preserves actual successor ancestry.  A positive test
must exhibit the source track, target track, copied or transported state,
preserved response predicate, and reset potential on both branches.  A
negative test should identify the first of those coordinates that cannot be
transported.  Merely recolouring `F13` and `F23` after the fork would be a
syntax-level repair and is not enough.

Reproduce with

```text
python3 04-computation/fc3_three_support_face_move_decorated_all_choice_closure_scout_20260803.py
python3 -O 04-computation/fc3_three_support_face_move_decorated_all_choice_closure_scout_20260803.py
```

and compare both transcripts with the stored output.
