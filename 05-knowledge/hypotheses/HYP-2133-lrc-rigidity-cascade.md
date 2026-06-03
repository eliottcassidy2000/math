---
id: HYP-2133
status: OPEN proof architecture; S589 gives exact small-tournament evidence for local source rigidity versus projection defect fibers
source: codex-2026-06-03-S589
related: [HYP-2134, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2128, HYP-2127, HYP-2126, HYP-2125, HYP-2124, HYP-2123, HYP-2122, HYP-2121, HYP-2120, HYP-2110, HYP-2109, HYP-2108, HYP-1977, THM-401, THM-400, THM-385, THM-381]
---

# HYP-2133: LRC rigidity is labelled fixed-point propagation, not symmetry size

## Claim

The proof-relevant rigidity in LRC and tournament quotients is:

```text
a marked local fixed point
+ labels that are natural under isomorphism
+ closure rules that propagate those labels
= a rigidity cascade.
```

Large automorphism groups, unmarked tournament classes, deletion decks, score
sequences, and balanced additive energy are not rigidity by themselves.  They
are shadows.  They become useful only when they preserve the marked predicate:

```text
observer is a source,
pair denominator is shielded,
endpoint owner is pinned,
relation lattice inherits to positions,
pincer jaws meet or export a labelled core.
```

In this frame, local rigidity is a fixed-point event around a basepoint.  Global
rigidity is the functorial cascade of that event through isomorphic copies,
deletion/addition maps, relation inheritance, endpoint circuits, and proof
automata.  Projection defect is exactly the failure of such a cascade after a
quotient forgets labels.

## Evidence From S589

`04-computation/tournament_rigidity_cascade_s589.py` enumerates unlabelled
tournaments through `n=6` and compares several quotient lenses on rooted
tournament classes.

The source root is exact:

```text
n=2: source_rooted=1  U(1)=1   delete_collisions=0  sources_fixed=True
n=3: source_rooted=1  U(2)=1   delete_collisions=0  sources_fixed=True
n=4: source_rooted=2  U(3)=2   delete_collisions=0  sources_fixed=True
n=5: source_rooted=4  U(4)=4   delete_collisions=0  sources_fixed=True
n=6: source_rooted=12 U(5)=12  delete_collisions=0  sources_fixed=True
```

So source deletion/addition is a genuine rigidity cascade: every parent has one
source-rooted child, deleting the source recovers the parent, and the source is
fixed by every automorphism.

Generic rooted views have defect fibers.  At `n=6`, the domain has `296` rooted
classes, but common quotients collapse them:

```text
full_rooted_class          unique=296  collisions=0
unrooted_plus_root_score   unique=196  collisions=100
underlying_unrooted_class  unique=56   collisions=240
split_profile_no_cross     unique=36   collisions=260  max_fiber=64
score_sequence_only        unique=22   collisions=274
delete_root_parent_only    unique=12   collisions=284
root_score_only            unique=6    collisions=290
```

The split profile records the root score and the two side subtournaments but
forgets cross edges between the sides.  Its large collisions are the tournament
analogue of LRC observer-coupling defect: local-looking side data is not enough
unless the incident/cross fiber is retained.

Global deletion decks are also not a magic rigidity source in this small range:

```text
n=3: 2 classes, 1 deck
n=4: 4 classes, 3 decks
n=5: 12 classes, 11 decks
n=6: 56 classes, 52 decks
```

Thus even global unmarked shadows can collide.  The rigid object for LRC is not
the unmarked deck; it is the labelled source/fold/endpoint cascade.

Tournament Analysis over quotient lenses is transitive:

```text
full_rooted_class
-> unrooted_plus_root_score
-> underlying_unrooted_class
-> split_profile_no_cross
-> score_sequence_only
-> delete_root_parent_only
-> root_score_only
```

This is a collision-order, not a theorem of usefulness.  It says exactly which
lenses lose how much observer-coupled payload.

## LRC Dictionary

```text
local fixed point      source observer, least-folded basepoint, endpoint owner
local algebra          relation_inherits, fold_triangle, pinch_pair_sum
clock rigidity         D=a+b with t=m/D; D|v shields the whole D-clock family
endpoint rigidity      small owner pin/peel, circuit midpoint resonance P(S)
pincer rigidity        witness front + blocker front + labelled escape ledger
automaton rigidity     L/M/R state cannot jump L<->R without terminal M
global cascade         source deletion, denominator-gate compression, endpoint peel
projection defect      unmarked/root-only/side-only fibers mix safe and unsafe states
```

The relation-lattice side is important.  THM-400 says balanced relations are
translation-invariant and observer-blind; unbalanced relations are
observer-coupled.  The Lean `relation_inherits` lemma says speed relations
become position relations at every time.  Therefore the true rigid algebra is
not "many relations"; it is the augmentation-graded relation lattice together
with the observer/basepoint labels.

## Rigidity Operator

Define a labelled state

```text
X = (quotient object, marked basepoint, predicate labels, forgotten-fiber log).
```

The rigidity closure `R(X)` repeatedly applies natural rules:

```text
source root          -> delete source / add unique source child
unbalanced relation  -> inherit position relation at all t
pair denominator D   -> enumerate pinch clocks m/D
shield D|v           -> kill the D-clock family for that observer
endpoint owner pin   -> off-lattice or half-radius peel
pincer jaws          -> return witness, positive escape, or labelled middle core
L/M/R automaton      -> route terminal middle circuits
isomorphism          -> transport every retained label, never just the shadow
```

A state is rigid for the LRC predicate if closure reaches either:

```text
accepting source/safe-box witness,
positive-measure escape,
or a finite labelled core whose remaining labels are explicit proof obligations.
```

A state is flexible if some quotient fiber contains mixed predicate values after
the retained labels are forgotten.  HYP-2121's observer-coupling defect and
S589's split-profile collisions are examples of such flexibility.

## Proof Route For n=14

The n=14 residual should be attacked as a rigidity-core problem.

1. Start with observer/basepoint source semantics (THM-381/385), not unmarked
   A000568 classes.
2. Apply denominator-gate compression (HYP-2122): low `D<n` killed,
   `D=n` survives, `D=n` killed, or no coherent denominator scaffold.
3. Route the multiple branch through endpoint-cover circuit positivity
   (HYP-2108) and small-owner pin/peel descent (HYP-2110).
4. Use pincer calculus (HYP-2123) to return a witness, positive escape, or
   minimal labelled middle core.
5. If a core remains, test whether it is a real object or a quotient artifact:
   add incident threshold words, endpoint-owner labels, and augmentation index
   before trusting any mixed bucket.

The conjectural terminal statement is:

```text
No labelled rigidity core survives all source, denominator, endpoint, pincer,
and middle-automaton closures for n=14.
```

## Assumption Challenge

Do not assume tournament vertices are runners.  For this hypothesis the useful
vertices are quotient lenses and proof obligations:

```text
source roots, rooted classes, split profiles, pair denominators, shield speeds,
endpoint owner clauses, relation-lattice grades, pincer middle states, proof gates.
```

The chosen quotient preserves the LRC predicate only when it retains the labels
that make a local fixed point propagate.  It destroys exact phase duration and
unlabelled isomorphism detail deliberately; if that destruction creates mixed
safe/unsafe fibers, the quotient must be lifted.

## See

`04-computation/tournament_rigidity_cascade_s589.py`,
`05-knowledge/results/tournament_rigidity_cascade_s589.out`,
`07-reflections/lrc-rigidity-cascade-fixed-points-s589.md`,
HYP-2134, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2128, HYP-2127, HYP-2126, HYP-2125, HYP-2124, HYP-2123, HYP-2122, HYP-2121, HYP-2120,
HYP-2110, HYP-2109, HYP-2108, THM-401, THM-400, THM-385, THM-381.
