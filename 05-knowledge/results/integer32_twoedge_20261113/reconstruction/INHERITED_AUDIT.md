# Recovery audit of hereditary slots, strong collisions, and one-slot wedges

**Status: independent proof review of inherited statements; no producer edits
or promotions.** Base commit `68d0d86bef480b6a36efc3d9864fd01ff3244bf5`.
The reviewed statements in PROOF.md were RESERVED at this audit; current
status is recorded in ../promotions.json. Exact local source hashes are
in `source-pins.json`.

## Inherited claims and their exact types

**Order-22 hereditary-slot ambiguity — PASS.** The module-run reversal
bijection works on every induced labelled subset, since intersection with a
module remains a module. Each numbered slot is broken into maximal runs in
that module, and every run is reversed without moving its vertex set, slot,
or outside neighbours. This preserves the unmarked slot count. No transitivity
is needed for this unmarked assertion. The inherited vertex-transitive,
non-self-converse 21-vertex core U, followed by a singleton, gives distinct
merged classes while preserving all induced Q. This audit does not enumerate
its Hamilton paths or assert order minimality.

The stronger equivalence between the entire hereditary Q table and labelled
cyclic triples uses Boussaïri–Ille–Lopez–Thomassé's interval-inversion theorem.
Its primary [Theorem 2](https://perso.ens-lyon.fr/stephan.thomasse/liste/bilt.pdf)
was reopened: equal cyclic-triple sets are equivalent to a sequence of
interval reversals. The squarefree logarithm in HEREDITARY_SLOTS.md correctly
converts induced H to the complete induced slot table; the odd-cycle
packing identity is inherited from Grinberg–Stanley. These are cited or
inherited mechanisms, not new results here.

**Strong order-24 first-moment ambiguity — PASS.** In a fixed exterior layout,
the j ordered nonempty module runs can be filled in exactly F_j(U) ways.
Automorphisms of a vertex-transitive U act on these fillings while fixing
run indices, so each fixed run has each vertex as its first, and as its last,
point equally often. Reversal changes a named endpoint into the other endpoint;
uniformity, rather than a false pointwise fixed-start claim, repairs that
change. Internal unordered edge uses are preserved by reversing the same run.
This proves every full-tournament per-slot edge, start, and end marginal at
every slot count. It does not prove every mixed edge/endpoint moment or the
same marginals on every induced subset.

For `C3[U,TT1,TT2]`, the three blocks are the maximal proper modules. A module
meeting two blocks forces all of the third, then all remaining vertices;
hence a proper cross-block module cannot occur. Distinct sizes force an
isomorphism to fix all quotient vertices. Reversing U alone cannot be an
isomorphism, and reversing the quotient cycle cannot be repaired by a
size-preserving permutation. Both tournaments are strong. Replacing the large
block by k ordered copies extends the same argument to 2^k distinct merged
classes at order 21k+3. The labelled six-vertex control in our checker uses a
self-converse C3 module to test the counting bijection; it is not offered as
another nonisomorphic strong collision.

**One-slot positive wedges recover up to converse — PASS relative to the
inherited directed-support decoder.** Legal single-vertex relocation is
connected by the explicit proof in `hamilton_path_relocation_connectivity_20261020/THEOREM.txt`.
Distinct relocations have a directed common deletion word. At n>=5 they
retain at least n−4 common pairs. For n=4 the only two gaps that could split
both common-word pairs are adjacent interior gaps; they cannot both admit
the moved vertex because its comparison to their common vertex would reverse.
For n=3 an interior and endpoint gap similarly conflict, leaving a common
edge for distinct legal moves. Thus the chain of edges in each Hamilton path
glues into a connected graph on all positive edges. Propagating the relative
orientation at its positive wedges recovers the directed used support up to
one global converse.

Recovering unused arcs is an essential separate step. The inherited
`session70-guards/REPORT.md`, G2–G5, uses Busch–Jacobson–Reid's nontraceable-arc
structure theorem to show that missing strong-block arcs are determined by
one-vertex-deleted reachability; across SCCs it restores the forward ordered
joins. The primary article remains a cited dependency, with its precise
scope recorded in that inherited audited source. We reopened its scanned
primary restatement but obtained no usable text extraction here, so do not
claim a fresh independent primary-source verification of BJR. Our independent
implementation reconstructs every ordinary representative through order six,
including its unused arcs. Those checks corroborate, rather than replace,
the cited all-order deduction.

**Two optional slots and the quadratic jet — PASS, inherited.** Rédei's path
on the complementary induced tournament makes every locally directed triple
observable in two slots. A negative two-slot answer means a local sink/source
at the common vertex, unlike a negative one-slot answer. The inherited
fixed M−1 spanning-tree queries and information bound are for labelled
orientations modulo converse and yes/no queries. They are not an unlabelled
class query bound, nor a bound for queries returning unrestricted integers.

## Controls, repairs, and scope

Normal and optimized runs of the new independent checker agree. It imports
no old producer code; only immutable prior canonical representatives and H
metadata are reused. Literal permutations give 1,332 Hamilton paths across
all 74 ordinary representatives at orders 3–6. Every one-slot support is
connected, propagated, and decoded including dark arcs; independent ordered
subset/path products give exactly the local two-slot wedge field.

An initial unpublished checker failed the one-slot round trip because its
SCC sorting key iterated the same list that Python temporarily empties while
sorting. Taking an immutable tuple snapshot before sorting repaired the
implementation. This was not a counterexample to the inherited decoder and
no mathematical claim changed. Both completed receipts include the repaired
routine and every gate remains active under `python -O`.

The crucial surviving distinctions are: labelled versus relabelled data;
full-tournament marginals versus induced marginals; separate marginals versus
joint marks; one-slot positive wedges versus two-slot local beta; and vertex
path-use observers versus the native ordinary-class dart matrix C. The new
local realization/repair theorem applies only to its declared beta carrier.
