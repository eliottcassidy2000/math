# Borel Anti-Diagonalization And Embedded Maximality, S673

The useful correction is this:

```text
diagonalization creates an outside object;
Borel anti-diagonalization asks whether an invariant outside selector can exist.
```

Those are not the same move.  Cantor-style diagonalization is observer-coupled:
given a presentation of a list, flip along the presentation and build something
off the list.  Friedman's Borel diagonalization theme is harsher.  If the
selector is Borel and invariant under natural reorderings of the presentation,
then some presentation catches it; the chosen value becomes a coordinate of the
very object it was supposed to escape.

That is exactly the repo's recent pattern, just in a cleaner language.

## The Four-Level Selector Stack

The stack now looks like:

```text
raw diagonal:
  choose an outside witness

constructive named-order selector:
  give an explicit finite procedure for choosing outside

invariant/Borel anti-selector obstruction:
  prove no uniform invariant choice can stay outside every object

embedded recursion:
  attach the address and prove every allowed extension captures, drops rank,
  or pays tax
```

Constructive mathematics sits between raw diagonalization and embedded
recursion.  It is not just "be explicit."  It asks what information the witness
procedure is allowed to inspect, whether the chosen witness is invariantly
defined, and what happens when the ambient object grows.

## The Finite Toy

The S673 script uses a finite universe `U` and a named subset `A`.  An outside
selector wants `y in U\A`.  If a symmetry group `G` acts on `U`, an invariant
selector at `A` can exist only when some outside point is fixed by every
automorphism stabilizing `A`.

That gives a tiny address-tax measurement.

For `n=6`:

```text
ordered/trivial:   63/63 selectable, tax 0
path reflection:   56/63 selectable, tax 1
cyclic rotations:  54/63 selectable, tax 1
dihedral cycle:    36/63 selectable, tax 2
full symmetric:     6/63 selectable, tax 5
```

This is the whole story in miniature.  If the quotient forgets too much, there
is no invariant way to choose the "outside" point.  If enough anchors are named,
the selector exists, but outer extension immediately captures it:

```text
A=() -> select 0; extension names 0
A=(0) -> select 1; extension names 1
...
```

So an outside witness is not yet a proof.  A proof is an address plus a
transition theorem.

## Embedded Maximality Translation

S667 said:

```text
maximal(object, ambient embedding, allowed extensions)
```

S673 says the same thing in selector language:

```text
outside(object, selector, ambient embedding, allowed extensions)
```

A local maximum can be destroyed by adding a point in the next cut.  A local
outside witness can be destroyed by naming it in the next stage.  Both are the
same failure: the object was maximal/outside only relative to an ambient that
was too small or too unnamed.

The repair is not "choose harder."  The repair is:

```text
name the missing address and prove the extension rank cannot keep escaping.
```

## Tangible Incompleteness Reading

The finite toy is tiny and decidable.  But it shows where the proof strength
lives.  Each finite closure sequence is transparent; the uniform theorem over
all finite stages, all allowed selectors, or all bad branches is the hard
object.  That is the same shape as Paris-Harrington:

```text
finite coloring condition: simple
relative-largeness boundary: tiny
uniform witness function: beyond PA-provably recursive growth
```

Tangible incompleteness is not mystical here.  It is what happens when a
transparent finite rule is asked to work uniformly over every extension.

## LRC14 Translation

For LRC14:

```text
raw diagonal witness:
  find a lonely time in a sampled row

constructive selector:
  choose a witness from a clock menu or carry fiber

invariant obstruction:
  the witness must survive quotienting by residue/carry/owner symmetries

embedded recursion:
  every coherent +27 lift either stays in AP/Vstar/2AP, pays strict tax,
  or drops bad-child rank
```

HYP-2241's owner-private deletion bit is no longer just a clever classifier.
It is an address tax.  The proof target is to show that after paying that tax,
there is no invariant way for a bad fiber to keep pushing its witness into the
unnamed tail.

That is a much more theorem-shaped version of the current LRC14 route.

## A000568 Translation

For A000568, the deleted-card collision story now has the same form.  The
unpaired deck nearly chooses the class, but two regular strong collision
buckets remain.  The missing choice is the `L/U` side of the deleted owner.

So:

```text
raw deck = diagonal shadow
half-filter deck = address tax
endpoint n->n+1 theorem = recursive capture law
```

The next enumerator should test `half-filter trace + child-count profile`.
This would be the tournament analogue of adding extension rank to the selector.

## Why This Feels Foundational

The repo has been circling the same object from many directions:

- sum/product/fraction/recursion;
- owner/carry/deletion;
- ultrafilter side choices;
- Paris-Harrington bad branches;
- embedded maximality;
- outer extension usability;
- tournament endpoint enumeration.

S673 gives a compact grammar for the common mechanism:

```text
An invariant proof cannot merely point outside.
It must pay the address tax and survive recursive capture.
```

That feels like the right bridge from Borel anti-diagonalization to the concrete
LRC/A000568/unit-distance work.
