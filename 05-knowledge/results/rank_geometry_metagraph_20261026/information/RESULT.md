# What tournament information survives an unlabeling or a Hodge projection

2026-10-26. **PROVED:** the elementary counting, normalization and
information identities below. **FINITE-EXACT:** exhaustive independent
checks through order five. **OPEN:** a smaller unrestricted A000568
recurrence and an all-order classification by a small inventory of families.
These are repairs and connections of existing mechanisms, not a claim of
priority for canonical augmentation, Hodge decomposition or Shannon's rule.

## Inheritance: recover the question, repair the carrier

The March [fractal tournament codec](../../../../07-reflections/fractal-tournament-codec.md)
correctly isolates a zero-savings model for deletion at a fixed labelled
vertex: its incident fair bits are independent of the remaining labelled
tournament, hence also of that remainder's class. Its score-selected
experiments have a different conditioning law. Neither law implies that
scores reconstruct almost all class information. The March
[unlabeling principle](../../../../07-reflections/the-unlabeling-principle.md)
also needs to distinguish the logarithm of the number of classes from the
entropy induced by sampling labelled objects or path charts.

The June [intransitivity reflection](../../../../07-reflections/tournaments-are-intransitivity-among-n-things.md)
supplies the useful object: a strict complete binary relation. Its unit
flow/potential identification is false. Correcting that identification
also corrects [THM-895](../../../../01-canon/theorems/THM-895-helmholtz-identity-perron-census.md).
Closest current mechanism: signed, **labelled** triangle curl plus one
marked score or a Hamiltonian-path cut reconstructs all arcs
([all-order proof](../../metagraph_reconstruction_20261023.md)). Canonical
hostile: the two order-five classes below agree on scores and H. Corrected
near miss: scalar energy is neither the curl tensor nor the relation.
Least-used sidecar: the sampling measure and the normalization of arc flow.

The live board is: exact source; actual boundary relation; class fibres;
signed curl; primitive/clone quotient; permitted operations. The connection
to the other lanes is about what these maps forget, not a reduction between
Collatz, planar colouring and tournament enumeration.

## 1. Exact class masses in a Hamiltonian-path chart

For n>=1 fix the path 0->1->...->n-1. Its free nonpath arcs form a binary
chart with m=binom(n-1,2) bits (m=0 when n=1). Let C be an ordinary
tournament isomorphism class, H(C) its number of directed Hamiltonian
vertex paths, and a(C)=|Aut(C)|. The automorphism action on such paths is
free: fixing the ordered list fixes every vertex. Thus

    h(C) = H(C)/a(C)

is an integer, and precisely h(C) different chart words represent C.
An orbit of paths produces one word by putting the path in positional
order; conversely equal words give an automorphism between the paths.
Rédei's theorem supplies at least one path per class. Consequently

    sum_C h(C) = 2^m,
    A000568(n) = sum_(chart words w) 1/h(C(w)).

This is an exact enumeration identity, not a compressed recurrence: one
must still know the true fibres. A converse-merged class has fibre h if
self-converse and 2h otherwise. Forgetting that distinction loses the
unmerged A000568 count.

If W is a uniform random chart word, then

    Pr(C(W)=C) = h(C)/2^m,
    H_Shannon(C(W)) = m - sum_C (h(C)/2^m) log2 h(C).

The conditional entropy given C is log2 h(C). Uniform labelled sampling
instead assigns mass n!/(a(C) 2^binom(n,2)). Uniform sampling of classes
has entropy log2 A000568(n). These are three different measures. At n4
the four chart fibres are 1,1,1,5, already refuting the identification of
uniform chart sampling with uniform classes. The chart-to-class map is
many-to-one even though its port data can be retained exactly.

The [augmentation producer](../enumeration/README.md) makes the fibre
inventory constructive at every order without Burnside averaging. It
still represents the classes; the information identities here do not
remove that cost or imply black-core connectivity.

## 2. Scalar predictability versus structural information

There are at most binom(2n-1,n)<4^n sorted lists of n scores in
{0,...,n-1}. Some of those lists are not tournament scores, so this is a
safe upper bound. Since every class has at most n! labelled realizations,

    A000568(n) >= 2^binom(n,2)/n!.

For a uniformly sampled class C and its sorted scores S,

    H_Shannon(C | S)
      >= binom(n,2) - log2(n!) - log2 binom(2n-1,n)
       = n^2/2 - O(n log n).

Thus scores leave quadratically many bits of structural uncertainty even
when they predict a particular scalar statistic well. Adding H and the
automorphism order still leaves this asymptotic obstruction: each has
at most n! possible positive integer values. The exact integer inequality

    2^binom(32,2) > (32!)^3 binom(63,32)

already proves that (sorted scores,H,|Aut|) cannot identify every order-32
class. This is a bound against **reconstructing every class from those
statistics**, not against computing the aggregate class count by a short
formula.

A much smaller explicit collision occurs at n5. In lexicographic upper
pair order, bit1 means i->j for i<j. The minimum-under-all-permutations
codes8 and10 are nonisomorphic, but both have sorted scores
(1,1,2,3,3), H=9, and |Aut|=1. Therefore they also agree on c3=3,
gradient energy and residual energy. The full labelled curl tensor and
path-port incidence can distinguish what these scalars forget.

## 3. The exact Hodge normalization and the binary relation

Set a_ij=1 if i->j, a_ji=-a_ij, and d_i=sum_j a_ij=2s_i-(n-1).
Put x=sum_i d_i^2. The least-squares problem minimizes

    sum_(i<j) (a_ij - (phi_i-phi_j))^2,
    subject to sum_i phi_i=0.

Its normal equations are (nI-J)phi=d, hence phi=d/n. Therefore

    E_gradient = x/n,
    E_residual = binom(n,2)-x/n.

The residual has zero divergence and is orthogonal to the gradient. With
every triangle filled, the clique complex has no harmonic1 component.
It would be wrong to infer that the one-dimensional complete graph has
no cycle space.

The former THM-895 used phi=d/(2n) and x/(4n). A single directed K2 has
d=(1,-1). Its true potential (1/2,-1/2) exactly fits the unit arc. The old
potential only fits half that arc and leaves a residual of nonzero
divergence. Its old test verified the energy of an assigned potential,
not the least-squares property. The repaired checker tests the normal
equations and orthogonality independently with integer residuals.

For an oriented triple define c_ijk=a_ij+a_jk+a_ki. A transitive triple
has c=±1; a cyclic triple has c=±3. Consequently its cyclic indicator is
(c^2-1)/8. The classical score identity gives

    24 c3(T) = n(n^2-1)-3x,
    E_residual = (n-1)(n-2)/6 + 8c3(T)/n.

In particular an acyclic transitive tournament still has residual energy
(n-1)(n-2)/6. The transitive relation is represented by the **signs** of an injective
potential's differences; its unit flow is generally not those differences.
There is no contradiction: the sign map destroys magnitude. This is the
precise correction to the historical phrase "order plus intransitivity".
The excess over the transitive baseline does measure cyclic triangles,
but that scalar does not determine their incidence or the tournament.

For the general decomposition see Jiang, Lim, Yao and Ye,
[Statistical ranking and combinatorial Hodge theory](https://arxiv.org/abs/0811.1067)
(2011). The complete-graph specialization above is derived directly;
no literature priority or new open-problem reduction is asserted.

## 4. Maps, sidecars and cheap decisive tests

| Source -> target | Preserved | Lost / required sidecar | Decisive test |
| --- | --- | --- | --- |
| Chart word -> tournament class | Relation up to relabelling | Marked path; h weights for uniform chart sampling | n4 fibres1,1,1,5 |
| Class -> sorted scores,H,Aut | Those scalar values | Labelled curl and actual path ports | n5 codes8,10 |
| Unit flow -> least-squares potential | Divergence | Divergence-free residual | K2 rejects half-normalization |
| Injective potential -> orientation | Pairwise order | All positive gap magnitudes | TT3 has nonzero unit-flow curl |
| Relation -> primitive clone quotient | Quotient relation | Clone sizes and actual merge rule | Recompute after every extension |

The analogy to the source-preserving rank is limited but useful: a finite
observer may fail to preserve the next transition even when it predicts
coarse behaviour. The [rank lane](../rank/README.md) proves its own exact
obstruction; the entropy bound here proves no impossibility for the
single fixed integer7. Likewise a boundary compatibility permutation is
not a planar distance realization; that obstruction requires the actual
Gram relations in the [unit-distance lane](../unit_distance/RESULT.md).

## 5. A literal Frobenius boundary for the new planar reflection

The [geometry lane](../unit_distance/RESULT.md) finds a planar reflection
of six labelled hexagon vertices which swaps two declared pairs and fixes
the third pair pointwise. Its action on the four residual frames has no
fixed point. Could that same permutation preserve a tournament relation
on these six ports? It cannot.

In general let g be an involution on n vertices with r transposed pairs,
and define D_g(T) as the number of unordered pairs whose oriented relation
changes under g. Then

    min_T D_g(T) = r.

With the full antisymmetric ±1 matrices, the corresponding squared
Frobenius norm is 8D_g(T). This fixes the normalization of the word defect.

Each transposed pair's own arc is necessarily reversed, giving the lower
bound. Every other unordered pair is either fixed pointwise or belongs to
a two-element pair orbit without sign reversal after one full cycle.
Orient those orbits consistently to attain the bound. This is the signed
pair-orbit defect mechanism recovered in the
[previous packet](../../root27_machine_bridges_20261026/README.md), now
applied to the actual geometric reflection. Here r=2, so the minimum is2,
verified independently on all 32,768 labelled six-vertex tournaments.

The independent reviewer strengthened this to the complete polynomial.
If f=n-2r and s=(binom(n,2)-r-binom(f,2))/2, the s moving unordered-pair
orbits each have length two, and

    sum_T t^D_g(T) = 2^(r+binom(f,2)+s) t^r (1+t^2)^s.

Each forced pair contributes 2t, each fixed-point pair contributes2, and
each moving orbit contributes2+2t^2. These choices partition the actual
arc coordinates, so the product is literal. For the geometric reflection
n6,r2 this gives `512 t^2 (1+t^2)^6`, including exactly512 minimizers.
The checker verifies every coefficient, not only the minimum.

The map retains the same six vertices and the same permutation, and asks
whether a *strict complete binary relation* can also be invariant. It
does not turn a colour compatibility table into such a relation. A metric
symmetry, a fixed-point-free action on frames, and an automorphism of a
tournament are different predicates. The exact defect measures their
failure to coincide; it does not furnish a unit-distance colouring
obstruction or a source7 rank. Thus this six/four connection has an exact
map and a quantified obstruction without identifying the underlying
open problems.

## 6. Reproduction and scope

    python 05-knowledge/results/rank_geometry_metagraph_20261026/information/check.py
    python -O 05-knowledge/results/rank_geometry_metagraph_20261026/information/check.py
    python 04-computation/helmholtz_perron_toothpick_macmini_S115.py

The first checker independently enumerates all 1,100 labelled tournaments
of orders0..5 and all their fixed-path charts; canonicalization tries every
permutation. It checks integral normal equations, exact rational energies,
triangle counts, free path actions, class and chart masses, the K2 hostile,
the n5 collision and the order32 integer capacity bound. Order0 uses the
empty path and empty class convention. Regular tournaments (d=0) are
positive controls on which both old and new energy normalizations happen
to agree. The legacy script retains its floating-point Perron census;
its correlations are numerical evidence, separate from these exact gates.
