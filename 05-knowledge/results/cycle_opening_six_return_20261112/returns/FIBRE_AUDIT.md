# Independent audit of the macro fibres and canonical enumeration

**VERIFIED / PROOF AUDITED.** The full-count lane independently reviewed
[FIBRE_GRAPH.md](FIBRE_GRAPH.md), including the canonical-peeling upgrade,
automorphism statement, Johnson graph, adjacent-exchange refinement,
distance polynomials, and exact selected-family enumeration. The proof
passes in its stated strong-source scope. It is not a classification of
all H-neutral moves or an evaluation of all of A000568.

## Canonical reconstruction and isomorphism classes

Suppose every old indegree is at least one and the selected predecessor
p has indegree one. After adding z beaten only by p, exactly p and z have
indegree one. Every other old vertex receives an arc from z. In particular,
the unique predecessor of p now has indegree at least two, while z's
predecessor p still has indegree one. Thus the newest vertex is uniquely
characterized by its indegree and its predecessor's indegree.

This applies to every L birth and to both births of T, because the first
T birth leaves the old port at indegree one. A strong input remains strong,
so the required minimum indegree persists. Repeating the intrinsic deletion
recovers all 8m birth positions down to the prescribed original order n.
Every unrooted isomorphism commutes with that deletion and must preserve
all recovered positions. The first newborn's old predecessor recovers the
original marked source root. Stopping at n is important; no assertion that
the original source itself is irreducible under peeling is needed.

Each macro's private arc between its fifth and sixth newborns determines
whether the macro was A or B. Consequently every macro word is recovered
from the unrooted output. Different words are nonisomorphic, and, for a
fixed base order n and positive depth m, different rooted strong source
isomorphism classes also remain distinct. This is stronger than mere
injectivity of labeled construction words.

The minimum-indegree guard is necessary for this proof. For the weak
transitive triangle (6,4,0), rooted at its middle vertex, one birth produces
(6,12,0,5). Its three vertices of indegree one form a directed triangle;
each has an indegree-one predecessor, so the newest vertex is not uniquely
recognized by this test. The published ordinary-isomorphism promotion
correctly restricts to strong sources.

## Automorphisms and enumeration

At positive depth every output automorphism fixes all newborn positions
and restricts to an automorphism of the base fixing its marked root.
Conversely a rooted base automorphism extends by fixing every newborn:
each new comparison depends only on its port, which is the original marked
root or an already fixed newborn. Restriction and extension are inverse,
proving Aut(output) is isomorphic to Aut(T,r).

Thus a fixed rooted strong source yields exactly 2^m ordinary classes at
order n+8m. Distinct counts are indexed by a+b=m-1 because93 and97 are
coprime and q0 is positive. Their multiplicities are exactly
`2 binomial(m-1,a)`, not just counts of possible addresses. The two final
letters give distinct continuation states and distinct ordinary objects
at the same H. Summing over m gives `x^n/(1-2x^8)` for this one selected
source family, including its single unrooted base at m=0.

The triangle-root stabilizer is trivial. At positive depth its outputs
are asymmetric, so each ordinary class supplies N! labelings. Hence the
stated `2^m N!` labeled count and per-H binomial distribution are correct.
The positive-depth qualification is essential: at depth zero the unrooted
triangle itself has three automorphisms.

Across rooted source classes, the injection is asserted at fixed base
order n and positive depth. One must not sum such formulae over arbitrary
base orders and arbitrary depths without a nonoverlap argument: an output
could itself be selected as a later base. Nor may the depth-zero term count
distinct root orbits as distinct unrooted sources. The current statement
retains these boundaries.

## Exact neutral graph and distances

The complete H value depends on the multiset of the first m-1 letters;
the final letter chooses the output type. Therefore the only neutral
one-letter move is the final flip. The only neutral two-letter move swaps
one A and one B in the prefix. A prefix change together with a final change
cannot repair the changed exponent pair.

Representing A positions by an a-element subset S and retaining the final
letter gives exactly `J(a+b,a) Cartesian-product K2`. A replacement of one
chosen position changes one missing element, so the distance is
`|S minus S'| + indicator(final letters differ)`. This proves the stated
vertex count, degree ab+1, edge count, and diameter min(a,b)+1. It also
proves the distance polynomial

    (1+z) sum_j binomial(a,j) binomial(b,j) z^j.

The boundary cases a=0 or b=0, including m=1, reduce to K2 and satisfy all
formulas. All edges retain the newest physical root and its sole predecessor.
Only the prefix exchanges preserve the entire closed L/T state; a final
flip preserves H but changes that state.

When prefix exchanges must be adjacent, increasing A positions s_i obey
`0≤s_i-i≤b`. Their differences encode an a-by-b Ferrers diagram (reverse
the row order for the usual decreasing-row convention). An adjacent
exchange changes one cell. This proves the distributive-lattice cover
interpretation and distance `sum_i |s_i-t_i|`, plus the independent final
flip. The diameter is ab+1. From an extreme word its distance polynomial
is `(1+z)` times the rectangle's area-generating Gaussian binomial.
This formula is correctly stated from an extreme vertex, rather than
asserted for arbitrary vertices of the nonregular refinement.

The local diamond proof also has the needed continuation: an exchanged
pair lies inside the prefix and therefore has a following macro. Its
letter gives ABA=BAA or ABB=BAB, rejoining the complete L/T state before
any subsequent continuation.

For six macros and a=2,b=3, the verified parameters are20 vertices,
degree7, unrestricted diameter3, and adjacent diameter7. Identifying
2-subsets of a5-set with edges of K5 preserves intersection adjacency;
disjointness is the usual Petersen presentation. This does not identify
the five source cores with K5 vertices.

## Physical realization and verification boundary

A macro-letter change alters exactly four private internal arcs. Later
birth comparisons depend on the same physical last root regardless of the
letter, so they create no additional differences. Consequently the stated
4-times-Hamming arc distance is exact under the construction's canonical
birth alignment. It is not a claim about minimum arc-edit distance under
arbitrary vertex relabeling. Additional neutral moves outside the selected
four-/eight-arc rules are not excluded.

The finite receipt accurately checks every word through depth seven by
exact matrices, graph construction, private-arc recovery, and canonical
peeling before and after relabeling. Direct arc differences are checked
for all within-H-fibre pairs through depth five; that qualifier was repaired
during review. The all-order proofs use canonical deletion, the independent
L/T laws, and the subset representation; no finite bound is extrapolated.

## Additional review: the full L/K binary sublanguage

The subsequent [BINARY_LANGUAGE](BINARY_LANGUAGE.md) also passes independent
review. General unique decoding and the automorphism mechanism are inherited
from [adaptive ANCESTRY, A1–A3](../../adaptive_metagraph_reformulation_20261103/ANCESTRY.md);
the cheap two-root recognizer is a specialization, not a new general
ancestry theorem.

For K, the graph birth is the same as L but the old root remains marked.
A path starting at that root must immediately use the newborn; its start
count and deletion count are both d. Hence

    K(x,b,d)=(2x+2b,0,d),       T=chronological KL.

These are physical marking identities, not just matrix equalities.
Every internal L/K choice is recovered from the next chronological port;
the output mark recovers the final choice. Forgetting the final mark
therefore forgets exactly one binary letter. This proves rooted count2^c
and ordinary count2^(c−1) at positive cost c, with the distinct final root
orbits fixed by the inherited automorphism argument. The ordinary series
is x^n(1−x)/(1−2x); the rooted series is x^n/(1−2x). The multiset full-H
histogram may be divided by two because every graph has exactly those
two final markings, regardless of additional numerical state collisions.

L/T words expand uniquely into L and KL tiles, equivalently binary words
with no consecutive K and final letter L, apart from the empty word.
Thus their ordinary coefficient at cost c is F_(c+1), with generating
function x^n/(1−x−x²). This Fibonacci count is neither the larger classical
unique-Hamilton-cycle count nor the inherited unrestricted adaptive
triangle-descendant count. A gap1/gap2 description uses the initial
marked root as time0; it is not literal vertex-label subtraction for an
arbitrary n-vertex seed.

For C3, one birth gives H5. At two births LL has state (4,2,5), while KL
has (4,1,5); both have H9 and their different defect responses agree with
the two existing five-core representatives. This is a literal construction
link, not a cardinality identification.

Finally [ANCESTRY_BOUNDARY](../compiler/ANCESTRY_BOUNDARY.md) records why
this enlarged binary language still cannot itself finish the ready return
at97 or293: every positive-length strict birth history has two indegree-one
vertices, while the inherited low-readiness classification forces every
strong1/2-ready representative at those targets to have a unique one.
The parent-changing module replacement leaves the restricted language
while remaining an ordinary adaptive construction.
