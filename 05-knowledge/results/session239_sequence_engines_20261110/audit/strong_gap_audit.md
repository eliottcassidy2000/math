# Audit: sharp strong sevenfold symmetry gap

Status: PROVED, with the eight finite seven-circulant inputs independently
verified by complete permutation enumeration. Existing audit findings remain
unchanged.

## Statement audited

If a finite strong tournament T has an automorphism of order7 and more than7
vertices, then H(T)>=2163. Equality is attained by the nine-vertex tournament
C3[Q7,1,1], where Q7 is a seven-circulant with H=175. This is a sharp bound;
no exhaustive classification of all equality cases is asserted here.

## Proof audit

An automorphism of order7 has singleton and seven-element orbits, with at
least one seven-element orbit. Let s be the number of moving orbits.

If s>=2, THM-5094,
`01-canon/theorems/THM-5094-cyclic-orbit-profile-cost-and-sharp-sevenfold-path-floor.md`,
gives H(T)>=175^s>=30625>2163. This alone closes the branch. The stronger
H(T)>=175^s+14 also follows: equality in THM-5094 forces an ordinal sum of the
orbit blocks, impossible for a strong tournament with at least two blocks.
Both H(T) and175^s are odd multiples of7, so the strictly positive difference
is at least14. Oddness uses Redei parity, and7-divisibility uses the free
action on ordered Hamilton paths. Do not attribute the extra14 to the product
bound alone without these additional steps.

If s=1, denote the moving orbit by R. Every other vertex is fixed, so its
orientation towards R is constant. Thus R is a module. Contract R to a vertex
r. The quotient is strong: projecting any directed path in T gives a directed
walk in the quotient, and repeated vertices can be removed. The quotient has
at least two vertices, and in fact at least three because a two-vertex
tournament is not strong.

Every vertex r of a nontrivial strong tournament belongs to a directed
triangle. Indeed let A be its out-neighborhood and B its in-neighborhood;
both are nonempty. If no triangle r->a->b->r exists, every edge between A and B
points B->A. No vertex of A can then reach r or B, contradicting strong
connectivity. Choose the two fixed vertices a,b of such a triangle. The
induced tournament on R union {a,b} is C3[R,1,1].

Hamilton-path count cannot decrease on adding a vertex v to a tournament.
For each old path, insert v at the start if v beats its first vertex, or at
the end if its last vertex beats v. Otherwise the path starts at a vertex
beating v and ends at one beaten by v, so at some consecutive pair x,y it
has x->v->y and v can be inserted there. Choose the first available slot by
any fixed rule. Deleting v recovers the old path, proving injectivity.
Iterating yields H(T)>=H(C3[R,1,1]).

The action on R is regular cyclic, so after choosing a cyclic gauge its
orientation is specified by one member of each pair {1,6},{2,5},{3,4} of
nonzero differences modulo7. These are all8 possible labelled circulants.
The independent exhaustive computation returns six values2163 and two2205
for their triangle inflations. Therefore the minimum is2163. Each inflation
is strong, because its cyclic quotient lets every vertex reach every other
vertex through the other two blocks. It inherits the order7 rotation, so the
minimum witness has the required properties.

## Independent finite certificate

`strong_gap_audit.py` imports none of the researched code. It constructs
all8 circulants directly from their modular difference sets, constructs each
nine-vertex triangle inflation as an arc set, and checks every one of the9!
ordered vertex permutations for a directed Hamilton path. It independently
counts ordered nonempty bipartition path covers and verifies

    H(C3[R,1,1])=3 H(R)+F2(R).

F2 is1638 for every case; internal H is175 in six cases and189 in two. Both
normal and optimized Python executions pass and have identical output.

Reproduction:

    python3 05-knowledge/results/session239_sequence_engines_20261110/audit/strong_gap_audit.py
    python3 -O 05-knowledge/results/session239_sequence_engines_20261110/audit/strong_gap_audit.py

## Inheritance and scope

THM-5094 provides the orbit-profile product and equality classification.
The seven-circulant minimum175 is inherited, not newly discovered.
THM-4987,
`01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md`,
records weaker module-interaction costs under a full-product orbit-module
hypothesis. Its normalized105 bound for a triangle through P7 does not supply
this exact inflation count or this arbitrary strong seven-symmetric bound.
A targeted canon search found no previous mathematical2163 statement; this
is not a comprehensive literature priority check. The new theorem does not
classify all tournaments or settle unrestricted Chan-Pak completeness.
