# A shared apex forces new cycles: an all-subsets gluing inequality

The [older Boolean Eulerian bowtie work](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/eulerian-boolean-kernel-overnight-hexagon-sep05.md)
already distinguishes a support-level move from its multiplicity. This
note asks the corresponding directed question: can two apparently
independent cycles share a vertex without forcing new counts? The
answer gives a different proof for the first forbidden value.

Use the [Rédei–Berge cycle-packing formula](https://arxiv.org/abs/2307.05569):
the Hamiltonian path count H(T) is the sum of 2 to the number of cycles
over all vertex-disjoint packings of nontrivial odd directed cycles.
The empty packing contributes 1. In particular restricting to induced
subtournaments gives an immediate positivity bound.

## From two arms to a whole product, with arbitrary cross orientations

Suppose A_i are induced subtournaments whose vertex sets are pairwise
disjoint except for one common v. In each arm choose a directed triangle
v,u_i,w_i,v through v. Put

    a_i=H(A_i−v),
    b_i=H(A_i−{v,u_i,w_i}),
    e_i=H(A_i)−a_i,

using H(empty)=1. Then for any orientation of the cross arcs,

    H(T) ≥ 2 Π_i(a_i+b_i) − Π_i a_i
           + Σ_i(e_i−2b_i) Π_(j≠i) a_j
           + 2 Σ_(i<j) b_i b_j Π_(l≠i,j) a_l.                 (1)

Each correction term on the right is nonnegative. This bound is about
literal mixed cycles, not an independence or modular assumption.

**Proof.** Packings entirely within the arms may use v either nowhere
or in exactly one arm. Their weight is Πa_i+Σe_iΠ_(j≠i)a_j.
The chosen triangle in arm i and all possible remaining packings show
e_i≥2b_i. For *every* subset I of at least two arms, the induced
tournament on v and the two chosen vertices from each arm is strongly
connected via v. It therefore has a spanning directed cycle, of odd
length 2|I|+1. Here is the elementary cycle argument: in a longest
cycle any outside vertex must point uniformly to the whole cycle or
uniformly away, else one can insert it. Strongness then supplies
an edge from the vertices beaten by the cycle to those beating it;
insert those two vertices consecutively to contradict maximality.
Fix one such odd cycle and complete it with arbitrary packings on
the unused vertices separately in every arm, avoiding v. Its
contribution is at least 2Π_(i∈I)b_i Π_(i∉I)a_i. Distinct subsets
give different selected cycles; within each counted packing this is
the only cycle meeting different arms.

Each pair of arms additionally forces a **mixed triangle** inside its
chosen five-vertex union. If there is none through v, then w_j→u_i
and w_i→u_j must hold, to avoid the cycles v→u_i→w_j→v and its
counterpart. If u_i→u_j, the cycle u_i→u_j→w_j→u_i is mixed;
if u_j→u_i, the reverse-label counterpart is mixed. It is distinct
from the selected spanning 5-cycle. Use the same residual lower
bound on the free packings (leave any additional unused chosen
vertices as singletons). Different pairs are again disjoint cases.
Summing all terms and expanding the product gives (1).

For k bare triangular arms, (1) reads

    H(T) ≥ 2^(k+1)−1+k(k−1).

For k=2 this is sharp at 9, regardless of the unspecified cross
arcs. The forced 5-cycle explains why two triangles meeting in one
vertex cannot behave as an additive model with H=1+2+2 or even
H=1+2+2+2. A directed triangle is not an independent unit of
Hamiltonian count. For exactly two arms the inequality is sharp
even with arbitrary achieved multiplicative factors a_i=b_i:
attach extra factors below the minimal five-vertex core in ordinal
order, assigning them to the two arms.

## An elementary packing proof for the first hole

If H(T)<9, two triangles cannot be disjoint (their packings already
have weight 9), and cannot meet in exactly one vertex by (1).
Three distinct triangle supports pairwise sharing at least two
vertices either share a fixed pair or lie on four vertices. In the
first case, the common directed edge and the three other vertices
form the cyclic join of two singletons and a three-vertex tournament.
Every nonempty subset of those three vertices of odd size gives an
odd directed cycle through the common edge: use a Hamiltonian path
in that subset. This gives four odd cycles and H≥9. In the other
case four vertices have at most two directed triangles: the number
of transitive triples is the sum of binomial(outdegree,2), at least
two since the four outdegrees sum to six.

So there are at most two triangles; if there are two, the union of
their vertices has size four. No longer odd cycle is possible.
Indeed its induced vertex set is strong, and every vertex in a
strong tournament lies on a directed triangle: otherwise its
in-neighbors all beat its out-neighbors, contradicting strongness.
A strong subset of size at least five could not have all its
vertices in the union of the at most two available triangles.
All odd-cycle packings therefore give H=1,3,5 only. This re-proves
the exclusion of 7 without an order floor or tournament census.
It does not give the exclusion of 21 by the same argument: the
latter is already known by the [repo's strong-order proof](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).

Contrast 63 directly. In the one-apex tournament on the transitive
spine 1001100 (1 means the apex beats that vertex), all cycles use
the apex. Counting by their first/last spine vertices gives 8
triangles, 17 five-cycles, and 6 seven-cycles. Thus H=1+2·31=63
on a single strong tournament, not by factoring through the
nonexistent 7. The factorization of the count need not be a
decomposition into independent cycles or tournament components.
