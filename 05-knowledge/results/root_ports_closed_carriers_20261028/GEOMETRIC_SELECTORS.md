# An eight-port geometric obstruction before six-colour forcing

**Status: PROVED elementary; FINITE-EXACT controls.** Same-agent proof review.

Let G be any finite graph realized by distinct points of the Euclidean plane
with each edge of length1. Nonedges may also have distance1. Designate eight
vertices as ports, pair them into four disjoint pairs, and prescribe a proper
six-colour assignment on the ports which is constant on each pair.

**Theorem.** If G has at most four other vertices, the prescribed assignment
extends to a proper six-colouring of G.

This applies to all105 pairings and therefore to the56 four-direction
matchings in the inherited [oriented Fano bridge](../free_bits_seven_fiftysix_20261028/README.md).
It works over arbitrary real coordinates, including extensions beyond the
earlier cyclotomic fields. It concerns this matching-constant observer;
assignments using five or six distinct boundary colours are outside its scope.

## Proof

Contract equal-coloured boundary vertices into g<=4 groups; keep every
interior vertex separate. Properness makes this quotient loop-free. Add all
edges between its boundary groups. A six-colouring of this augmented quotient
would use distinct colours on those groups, and a permutation of the six
colour names would match the prescribed boundary colours exactly.

The quotient has at most8 vertices. Any graph on at most8 vertices that
requires at least7 colours contains K7. For the only nontrivial case,8 vertices:
two disjoint nonedges give a six-colouring by merging their two pairs; a
triangle of nonedges gives one by merging its three vertices. If neither
occurs, all nonedges form a star (or are empty), and deleting the star centre
leaves K7. The claim for fewer vertices is immediate or follows by adding
isolates to reach8.

Such a K7 in our quotient would contain at least3 interior vertices. Four
mutually adjacent interior vertices are impossible in the unit-distance
plane. Thus it would contain exactly3 interior vertices forming a unit
equilateral triangle and all4 boundary groups. Each group is then exactly
one original pair of distinct ports. To be adjacent to all3 triangle vertices,
each pair must contain a port adjacent to at least2 of them.

Outside a fixed unit triangle there are exactly3 points at distance1 from
at least2 of its vertices: the reflections of the third vertex across each
side. Two unit circles centred at the side endpoints have exactly those two
intersections, one of which is the third triangle vertex. Four disjoint port
pairs would require4 distinct exterior double-contact points, a contradiction.
The quotient is six-colourable, giving the desired extension.

## Scope compared with the existing stronger field-specific theorem

[THM-4614, three-colourable contexts in Q(zeta55)](../../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md)
already handles arbitrarily large three-colourable exteriors of its fixed
eight-port core and **all** proper six-colour boundary assignments. It uses
the exact nine-point double-contact locus in that field. We do not replace
that theorem with a weaker size-limited version. The present trade is explicit:
arbitrary plane ports and arbitrary fields, but only matching-constant
assignments and at most4 private interior points.

The abstract graph K3 joined to8 independent ports defeats a prescribed
four-colour pairing with just3 interior vertices. It is not a unit-distance
realization. This hostile isolates the geometric step rather than confusing
an abstract quotient obstruction with a realizable selector.

## Exact controls and the next search

[check_selectors.py](check_selectors.py) takes all495 choices of4 interior
vertices among12 triangular-lattice points, all105 matchings of the remaining
ports, and all15 partitions of the4 pair-colour groups. Of779625 candidate
patterns,765514 fail proper boundary colouring; all14111 proper patterns
extend. Among the latter,6372 use one of the56 bridge matchings. All arithmetic
is exact in the basis1,exp(pi i/3), with squared distance a^2+ab+b^2.
The receipt also checks the abstract hostile and the three exterior contacts.

For this observer, a nontrivial orientation selector needs at least5 private
interior vertices. A precise next question is whether any such five-interior
selector exists for any eight real ports. Its augmented quotient can now have
9 vertices; the K7-only argument no longer covers all seven-chromatic graphs
(K4 joined to C5 has9 vertices and needs7 colours). That graph is only a
combinatorial candidate: its boundary contractions must be expanded to eight
distinct points and every unit-distance equation authenticated. This is an
OPEN geometric search, not an exhibited unit-distance forcing gadget.
