---
source: codex-2026-06-03-S589
status: synthesis plus exact small-tournament audit
tags: [LRC, rigidity, fixed-point, observer-coupled, tournaments, pincer, HYP-2133]
---

# Rigidity cascade: fixed points, defect fibers, and the labels that permeate

The user asked where rigidity appears, and specifically for local rigidity
around a fixed point versus global rigidity that cascades through symmetrical
objects like tournaments or LRC.  The answer that emerged is:

```text
rigidity is not symmetry;
rigidity is labelled stabilizer propagation.
```

High symmetry can help, but it can also be a decoy.  A quotient can have a
beautiful automorphism story and still forget the one label that decides the
predicate.  For LRC, that forgotten label is usually the observer/basepoint
fiber.

## The sharp local fixed point: source root

The cleanest tournament fixed point is a source.

If a rooted tournament has the root as a source, then:

```text
the source is fixed by every automorphism,
deleting it is canonical,
adding it back is unique,
and the parent class is recovered exactly.
```

The S589 computation checks this through `n=6`:

```text
n=2 source roots = 1  = U(1)
n=3 source roots = 1  = U(2)
n=4 source roots = 2  = U(3)
n=5 source roots = 4  = U(4)
n=6 source roots = 12 = U(5)
```

This is the tournament shadow of THM-381: LRC loneliness is observer-is-source
in the marked observer tournament.  The source is not just a vertex with high
score; it is a fixed point whose local label survives isomorphism and deletion.

That is the first rigidity rule:

```text
local source fixed point -> exact global source-deletion cascade.
```

## The tempting false rigidity: unmarked shadows

Now compare arbitrary roots.  At `n=6`, there are `296` rooted tournament
classes.  If we look at weaker shadows:

```text
unrooted class only        -> 56 values
root score only            -> 6 values
delete-root parent only    -> 12 values
score sequence only        -> 22 values
side split profile only    -> 36 values
```

All of those have large fibers.  The side split profile is the most revealing:
it records the root score and the two internal subtournaments on the out-side
and in-side, but forgets cross edges between the sides.  At `n=6`, it has max
fiber `64`.

That is exactly the LRC failure mode.  You can know the observer's two sides and
still not know the incident/cross coupling that makes the observer safe or
unsafe.  The missing cross edges are the tournament version of threshold words,
endpoint owners, and augmentation-index labels.

The global deletion deck is not a rescue.  S589 finds deck collisions already:

```text
n=3: 2 classes, 1 deck
n=4: 4 classes, 3 decks
n=5: 12 classes, 11 decks
n=6: 56 classes, 52 decks
```

So the mantra should be:

```text
unmarked global shadows do not certify marked local predicates.
```

They may cache, stratify, or guide a search.  They do not decide LRC safety
unless the observer-coupled fiber is retained.

## Relation rigidity: algebra becomes geometry

The fold material says the same thing in arithmetic language.

Balanced relations are observer-blind.  A 4-term relation

```text
a + b = c + d
```

survives translation and lives in the difference geometry.  It is real
structure, but it does not anchor the observer.

Unbalanced relations are observer-coupled.  A fold

```text
a + b = c
```

references the origin; after translating all speeds, it changes.  This is why
3-term folds carry the hard LRC signal and 4-term energy can be a harmless
shadow.

The Lean `relation_inherits` lemma gives the local rigidity mechanism:

```text
speed relation -> position relation at every time t.
```

HYP-2122 adds the clock version:

```text
D = a+b, t=m/D, and D|v  -> v*t is an integer.
```

So a speed divisible by `D` shields the whole `D`-pinch family.  This is a true
local fixed-point phenomenon: the speed lands at the observer at every clock in
that family.  The global cascade is the denominator-gate ledger:

```text
low D killed / D=n survives / D=n killed / no coherent D-scaffold.
```

AP and `V*` are rigid boundary calibrations because their denominator gates
strip the lower clocks while leaving the `D=n` clock.  Unit-shift AP is loose
because it kills `D=n`; far-shift AP is loose because it lacks the scaffold.

## Endpoint rigidity: pins and circuits

Endpoint-cover work gives another local-to-global pattern.

A small endpoint owner is not just a participant.  It becomes an exact pin:

```text
w(k n +/- 1) = j u.
```

If the congruence fails, the component peels.  If it holds but the component is
too long for the half-radius, the component peels.  That is local rigidity
around an endpoint fixed point.

HYP-2108 turns many such local statements into one global circuit:

```text
one integer v must put all component midpoints near compatible arc centres.
```

The scalar `P(S)` is a global resonance pressure.  Local endpoint pins
permeate around the circuit; if any component refuses the pin, the cover fails
and a positive interval escapes.

This is the same pattern again:

```text
local pin -> circuit propagation -> witness/positive escape/core.
```

## Pincer rigidity: the two-front version

HYP-2123 reframes pinches as pincers.  A pinch alone is contact:

```text
D=a+b, t=m/D.
```

A pincer is contact plus the opposing blocker front:

```text
witness front:  pair clocks, unblocked jaws
blocker front:  shields, anchors, endpoint owners, Phi/CRT gates
meeting:        safe witness or positive escape
failure:        labelled middle core
```

The label is crucial.  A pair-sum clock without observer/basepoint, side,
residue, and endpoint-owner labels is not a proof object.  It is a shadow of a
proof object.

This is where the three-state automaton belongs.  `L/M/R` is a tiny model of
rigidity propagation: an obstruction cannot pass from left to right without
entering middle.  Terminal middle states should be attacked as cores, not
blurred into binary edges too early.

## The fixed-frame failure was already this

The old "almost fixed frame" thread now reads as a rigidity warning.  At
`n=3,4`, independent pair flips behave like a fixed frame.  At `n>=5`, the
frame cannot stay fixed: pair channels couple through the shared frame.

That is the same threshold as the perspective defect and split-profile
collisions.  Small systems let a root/perspective/frame masquerade as a
product.  Larger systems expose the missing coupling fiber.

This also explains the user's old perspective count:

```text
n=3 perspectives = 4 = structures on 4
n=4 perspectives = 12 = structures on 5
n=5 perspectives = 48 < 56 = structures on 6
```

The first two equalities are not a natural bijection.  They are zero-defect
coincidences before coupling appears.

## How to use this in the LRC proof

The proof engine should ask, at every quotient:

```text
what is the fixed point?
what predicate does it preserve?
what labels are natural under isomorphism?
what fiber did we forget?
does closure reach a witness, positive escape, or labelled core?
```

The clocks that matter are the ones made rigid by labels:

```text
D=a+b pinch clocks,
D|v shield clocks,
D=n delta clock,
endpoint centre-lattice clocks,
small-owner pin clocks,
source observer states,
terminal L/M/R middle states.
```

The clocks and shadows that are mostly cache layers:

```text
raw unmarked A000568 class,
score sequence,
root score,
delete-root parent alone,
side split profile without cross edges,
balanced 4-term energy alone,
global deletion deck alone.
```

The next proof target is a no-core theorem:

```text
after source semantics,
denominator-gate compression,
endpoint pin/peel,
pincer closure,
and L/M/R middle routing,
no labelled n=14 rigidity core remains.
```

That would be a much sharper statement than "the quotient is rigid".  It says
the only local fixed points that could trap a counterexample either propagate
to a witness or dissolve into positive measure.

## Artifacts

`04-computation/tournament_rigidity_cascade_s589.py`
`05-knowledge/results/tournament_rigidity_cascade_s589.out`
`05-knowledge/hypotheses/HYP-2133-lrc-rigidity-cascade.md`
