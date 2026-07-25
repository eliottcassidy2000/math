---
id: THM-2197
title: "Scalar chord coverage has a Boolean deficiency quotient, not an intrinsic mask tournament"
status: >
  PROVED + VERIFIED-EXACT. On a guard-safe thirteen-root fibre, adjoining
  unit masks acts exactly by intersection on their safe-sheet deficiency
  sets. Under arbitrary singleton-test continuations this Boolean
  semilattice is the minimal continuation quotient recognizing coverage.
  No rooted-dihedral-equivariant tournament orientation exists on the five
  mask vertices: an actual locally covering scalar mask configuration has a
  rooted reflection which swaps two pairs of masks, whereas no tournament
  can admit such an
  involution on those vertices. Covering and noncovering actual strict scalar
  mask configurations on one fibre can also have the same five-length
  profile. Chirality or coefficient labels can break the tournament symmetry,
  but do not replace the
  safe-sheet deficiency sidecar; deletion and phase propagation further
  require ownership/incidence, cyclic event order, and winding data. This
  does not eliminate any of THM-2192's remaining 244 residue profiles and
  is not a proof of LRC(14).
source: codex-2026-07-24-scalar-chord-tournament-no-go
depends_on:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
related:
  - THM-2183-order-join-is-an-exact-tournament-metric-product
  - THM-2195-transitive-quotients-exactly-control-universal-substitution-products
script: 04-computation/lrc14_scalar_chord_deficiency_tournament_thm2197.py
output: 05-knowledge/results/lrc14_scalar_chord_deficiency_tournament_thm2197.out
script_sha256: 7022c8b26d870282f5a89ab807f69abd0437bf669f7ea48561358526aa0ad3c9
output_sha256: 7d441307751bb4f919677bc13e10b620784f5060cd5a19befc876d0238db4f37
hash_basis: working-tree bytes (LF)
---

# THM-2197 -- scalar chord coverage has a Boolean deficiency quotient

Fix a guard-safe thirteen-root fibre from THM-2192.  Normalize its four
consecutive guard-unsafe sheets to a rooted block

```text
F subset R=F_13,              |F|=4,
B=R\F,                        |B|=9.                 (1)
```

A partial unit carrier is a finite indexed family of root masks

```text
cal A=(A_i),                  A_i subset R,           (2)
```

where a scalar unit mask is a singleton or a chord.  Put

```text
U(cal A)=B intersection union_i A_i,
Z(cal A)=B\U(cal A).                                 (3)
```

Thus `Z(cal A)` is the set of guard-safe sheets still lacking a unit owner.
The fibre coverage predicate is

```text
Cov(cal A) iff Z(cal A)=emptyset.                    (4)
```

## 1. The exact operation is a Boolean semilattice

If `cal A disjoint_union cal C` means that the two indexed mask families
are adjoined without identifying labels, then

```text
Z(cal A disjoint_union cal C)
   =Z(cal A) intersection Z(cal C).                  (5)
```

This is immediate from De Morgan's law, but it identifies the correct
operation: unit-mask coverage is recognized by the meet semilattice

```text
(2^B,intersection),                                  (6)
```

with the empty set as the accepting state.  It is an idempotent union
quotient, not an `l1` tournament-substitution product.

There is a sharp continuation statement.  Permit a continuation to adjoin
any finite family of safe-sheet singleton masks.  Define

```text
cal A equivalent cal A'

iff for every singleton continuation cal C,
    Cov(cal A disjoint_union cal C)
      iff Cov(cal A' disjoint_union cal C).           (7)
```

Then

```text
cal A equivalent cal A' iff U(cal A)=U(cal A')
                         iff Z(cal A)=Z(cal A').      (8)
```

The reverse implication follows from (5).  For the forward implication,
suppose for example that

```text
x in U(cal A)\U(cal A').                              (9)
```

Adjoin the eight singleton masks

```text
cal C_x=({y}:y in B\{x}).                            (10)
```

Then `cal A disjoint_union cal C_x` covers `B`, while
`cal A' disjoint_union cal C_x` still misses `x`.  The other direction of
the symmetric difference is identical.  Hence all `2^9=512` deficiency
states are pairwise distinguishable by singleton tests.

This minimality is deliberately scoped to arbitrary singleton-test
continuations.  It does not assert that all 512 states occur inside one
five-mask LRC branch, nor that eight extra singletons are a legal
continuation of that branch.  It says that any compositional summary which
is expected to remain valid under local sheet tests must retain the full
deficiency set.

## 2. A genuine scalar reflection forbids an intrinsic tournament

The obstruction already occurs on an actual root fibre, not merely in an
abstract matching.  Take

```text
H=1,                         y=3/17,
(v_1,v_2,v_3)=(1,2,3),
R_y={x in R/Z:13x=y},
x_k=(y+k)/13,                k in F_13.              (11)
```

All three divided-deep masks are inactive because
`||v_j y||>1/14`.  Thus `y in G intersection C_H` in THM-2192's notation.
The guard-unsafe and guard-safe sheet labels are

```text
F={11,12,0,1},               B={2,3,...,10}.         (12)
```

For the five positive thirteen-unit coefficients

```text
q=(40,12,25,92,105),                                  (13)
```

the strict danger masks `A_q={k:||q x_k||<1/14}` are

```text
{6}, {2,3}, {4,5}, {9,10}, {7,8}.                    (14)
```

They partition `B`, so this is a carrier of THM-2192 type `M`.  Each chord
has unoriented residue length one, and give the singleton the same invisible
length label one.

Now apply the rooted reflection

```text
r(k)=12-k mod 13.                                    (15)
```

It preserves `F` and `B`, fixes the singleton `{6}`, and swaps the chord
pairs

```text
{2,3} <-> {9,10},
{4,5} <-> {7,8}.                                     (16)
```

Suppose a rule assigned a tournament to the same five mask vertices using
only the rooted static chord carrier and commuted with every rooted
dihedral isomorphism.  The reflection `r` would be an automorphism of the
assigned tournament.  But if the arc between `{2,3}` and `{9,10}` points
from the first mask to the second, applying `r` makes it point from the
second to the first, a contradiction.  The opposite orientation gives the
same contradiction.  Therefore:

> **No rooted-dihedral-equivariant tournament orientation exists on the
> five vertices of the static scalar mask carrier.**                 (17)

More generally, a tournament has no nontrivial involutive automorphism:
an involution moving `a` to `b` would send the unique arc between them to
its reverse.  Equivalently, the automorphism group of every finite
tournament has odd order.  The reflection witness is therefore reusable
beyond the displayed matching.

This is the smallest possible symmetry mechanism: it uses one monomer and
two reflected pairs of chords in a nine-sheet cover.  A chirality choice
for the root cycle, or an external ordering of the coefficient labels,
breaks this particular involution.  Such a choice is a sidecar, not a
property of the unoriented rooted carrier.  Retaining the integer
coefficients also breaks the displayed symmetry, but then the winding data
has explicitly not been quotiented away.

## 3. The five-length profile does not recognize coverage

Keep `H=1`, `y=3/17`, and the singleton coefficient `40`.  Replace the four
chord coefficients in (13) by

```text
(12,25,105,27).                                      (18)
```

Their masks are

```text
{2,3}, {4,5}, {7,8}, {8,9}.                          (19)
```

Every visible chord again has length one and the singleton is still `{6}`
with length label one.  Thus (14) and (19) have the identical residue-length
profile

```text
(1,1,1,1,1).                                         (20)
```

The first carrier has deficiency `emptyset`.  The second double-owns sheet
`8` and misses sheets `10`, so

```text
Z(second carrier)={10}.                              (21)
```

Both lists come from exact strict scalar danger inequalities on the same
root fibre with positive distinct thirteen-unit coefficients.  Consequently
the hostile control is stronger than a freely placed abstract chord
diagram: the `244/252` length-profile quotient in THM-2192 is a necessary
static filter, but it cannot decide root-sheet coverage.

## 4. Vertices, observable, ties, target, and lost data

The tournament audit is:

```text
candidate vertices: the five unit masks;
pairwise observable:
    rooted two-mask chord position, intersection, or cyclic order;
ties:
    rooted reflections can exchange a pair, and disjoint masks have no
    intrinsic winner;
target:
    every one of the nine guard-safe sheets has at least one owner;
preserved by a mask tournament:
    at most chosen pairwise comparisons among the five masks;
lost:
    unary mask-to-sheet incidence, uncovered-sheet identity, multiplicity,
    chirality, phase-event order, and coefficient winding;
static minimal continuation sidecar:
    Z subset B, equivalently the covered-sheet union U;
faithful carrier:
    the rooted, coloured mask--sheet incidence graph.                 (22)
```

One can orient every cross edge of that coloured bipartite carrier by

```text
mask -> sheet iff the mask owns the sheet.            (23)
```

Then coverage says exactly that every safe-sheet vertex receives an arc
from the mask colour class.  Orientations within either colour class are
irrelevant and generally noncanonical; the colour partition and rooted
sheet cycle remain load-bearing sidecars.  Calling an arbitrary completion
of (23) a tournament adds no information.

## 5. Why the static quotient is not the propagation state

Equation (5) is exact for adjoining masks.  Multiplication-by-thirteen
propagation changes a mask by deleting an old endpoint and inserting a new
one.  A deficiency set alone cannot update after deletion: it does not say
whether the departing endpoint had another owner.  For example, on two safe
sheets `x,y`, the labelled families

```text
({x},{x,y})       and       ({x},{y})                 (24)
```

have the same covered union.  After deleting the first labelled mask, the
first family still covers `x` and the second does not.  An updateable phase
state must therefore retain at least the safe-sheet ownership counts, and a
coefficient-wise update requires the labelled incidence sets.

Nor do static incidences determine which update occurs next.  The all-depth
problem additionally needs:

```text
the anchored cyclic order of boundary events,
the mask label losing or gaining each endpoint,
the active/inactive bits of the two shallower deep combs,
and the signed integer winding coefficients.         (25)
```

This explains the scope of THM-2183 and THM-2195 here.  Their exact
tournament product law concerns unlabeled arc-reversal distance.
Scalar coverage is a domination/union predicate, whose exact product is
(5).  A transitive quotient can make reversal distance additive only after
choosing an external order; a nontransitive quotient permits block
transport, even beyond the quotient-automorphism gauge.  Neither case
recovers the ownership current in (25).

The theorem therefore rules out a static tournament shortcut and identifies
the correct finite state to carry into the next computation.  It does not
exclude any additional residue profile, construct the phase transition
word, control the unique deepest blocker of THM-2192, or prove LRC(14).
QED.
