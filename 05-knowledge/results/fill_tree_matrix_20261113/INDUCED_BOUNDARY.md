# Five arithmetic tree edges that no pure birth construction can realize

**PROVED reduction plus FINITE-EXACT complete certificate.** These
exclusions concern ordinary Hamilton counts, with no selected-root,
readiness, defect or fixed-order requirement. They do not exclude
repairs that change old arcs or delete old vertices.

**Inheritance:** the November2
[induced-obstruction theorem](../metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md)
already proves the prime-source exclusions13→27,29→57 and31→63, together
with four larger branches and the failure of scalar-arrow transitivity.
Those exclusions are independently recovered here, not claimed as new.
The present union recursion also handles composite source counts and
records the full target-through73 table, including25→51 and35→69.

## The complete finite reduction

For an attained odd h let U(h) be the union of the induced Hamilton-count
sets of all finite tournaments with Hamilton count h. Include the empty
induced tournament, whose count is1. Let V(h) be the analogous union
restricted to nontrivial strong tournaments. Then

    U(1)={1},
    U(h)=V(h) union UNION_(ab=h, a,b>1) U(a)U(b),

where XY={xy:x in X,y in Y}. This is an exact recursion, evaluated in
increasing h. Write a tournament as its ordered SCC sum. Hamilton counts
multiply. Intersecting an induced subset with its components gives the
same product rule, with empty intersections of count1. Singleton
components contribute only1 and can be removed without changing this
set. Conversely an ordinal sum realizes each product in U(a)U(b).

For h<=73, Busch's bound H(C)^3>=5^(|C|-1) implies |C|<=8:
73^3=389017<390625=5^8. The inherited complete strong census through
order8 and this recursion therefore compute U(h) at every tournament
order for all37 positive odd targets through73.

The bank uses converse-merged representatives through order7 and
ordinary representatives at order8. Converse reverses each path and
preserves the induced-count set, so it is safe for this union. The same
convention would require weights when counting ordinary classes.

## The first failed arithmetic branch

The complete result is

    U(27)={1,3,5,9,15,17,27}.

No H27 tournament contains an induced H13 tournament. Thus no sequence
of vertex births takes any H13 source to an H27 target, regardless of
the number of births, intermediate counts, initial representative or
transitive padding. The binary edge13→2·13+1=27 is an impossible pure
growth operation even though both endpoint counts occur.

Structurally, H13 has one strong atom of order5. H27 has ten ordered
nontrivial SCC words: five strong atoms (two of order6, three of order7),
four ordered words consisting of C3 and one of the two strong H9 atoms,
and the word C3→C3→C3. Prime13 cannot assemble across multiple components.
Components of count3 or9 cannot contribute13 by insertion monotonicity.
The strong H27 atoms reduce under converse to three representatives;
their48 five-vertex subsets (6+21+21) contain no H13 atom.

The inherited exact multiplicities are N_n(13)=n−4 for n>=5, whereas

    N_n(27)=n(n−5)(n−4)/6,   n>=4.

The latter follows from

    F_27(z)=(2z^6+3z^7)/(1−z)^2
             +4z^8/(1−z)^3+z^9/(1−z)^4.

A cubically growing target fibre can therefore have no pure-birth edge
from a linearly growing source fibre. Multiplicity growth does not
supply the missing embedding.

## Complete small tree-edge classification and repaired routes

Within targets h<=73, the parent edges of TREE_BLUEPRINT.md fail exactly
at

    13→27, 25→51, 29→57, 31→63, 35→69.

Every other allowed parent/child pair in this range has an induced
witness, retained as literal rows and a subset mask or as a recursive
ordinal-product witness. Pairwise witnesses do not establish a coherent
choice of one representative working for all descendants.

Four failures have a useful physical replacement:

    9→27, 17→51, 19→57, 23→69.

In each case adjoin an ordered C3 component. This multiplies H by3 and
works for every source representative. For63 the analogous factorizations
3·21 and9·7 involve omitted counts; its strong carry atom must be supplied
separately. This is exactly where the inherited two-prime-lane reduction
retains63 as a primitive carry. It suggests combining binary and
multiplicative edges with genuine core-changing repairs, not requiring
one uniform arithmetic operation everywhere.

These five misses refute a pure-birth lift of the full arithmetic tree.
The obstruction is stronger than a unit-ready restriction: it applies
to ordinary unrooted counts. A different arithmetic tree or changes to
old arcs remain outside its scope.

## Reproduction and independent audit

Run `python -B 05-knowledge/results/fill_tree_matrix_20261113/containment_edges.py`
and repeat with `-O`. The universe is148 retained strong representatives,
23,512 induced subsets and the full SCC-product closure at targets<=73.
Each subset is evaluated in forward and backward endpoint coordinates.
The complete class bank is inherited from the reviewed
`actual_hamilton_multiplicities_20261021/session24-multiplicities/` packet;
source hashes are in the receipt.

The independent `matrix/containment_audit.py` enumerates87,388 individual
simple directed paths by depth-first search and recomputes the union
through nondecreasing factor multisets. It agrees on every subset, all37
target unions and all five failed edges. This checks calculations and
closure without regenerating the inherited isomorphism census. No new
unrestricted A000568 value is claimed.
