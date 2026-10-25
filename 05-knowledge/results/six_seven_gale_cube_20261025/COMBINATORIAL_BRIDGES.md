# Six and seven roots: free deck quotients, Fano codes, and a Paley fiber product

**Status: PROVED reductions and FINITE-EXACT incidence censuses, over C.**
The K3/Enriques and Campedelli constructions are classical. This packet adds
the maximal-free-deck bounds, the six-subset partition census, and the
complementary-kernel application. Its Fano code, 16/12/2 census, and seven-subset
calculations are concurrent confirmations of the incoming packets credited
below. It makes no historical novelty claim for the surfaces or binary codes.

## Inheritance and the live concept board

The closest proved mechanism is Theorem T1 in
`05-knowledge/results/session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md`,
independently checked in that packet's `positive-audit/TORIFICATION_AUDIT.md`.
For n distinct lines in P2 with no four concurrent, it constructs the normal
root cover X_n and its minimal resolution Y_n. The deck group is

    G_n = F2^n / <(1,...,1)>,    g_i = [e_i].

The branch inertia at a crossing indexed by S is the span of {g_i:i in S}.
The unique relation among all n generators is their total sum. A triple gives
2^(n-4) A1 nodes. In particular Y_6 is K3; Y_7 has q=0, p_g=7, K^2=16.

Before integration, live main acquired the audited
[Session60 seven-line theorem](../session60_six_seven_geometry_20261025/seven-lines/RESULT.md)
and the parallel [Fano-code packet](../six_seven_geometry_20261025/fano_codes/RESULT.md).
They recover the same Hamming/simplex quotient, freeness test, 30 kernels,
16/12/2 hourglass census and seven-subset counts. Session60 additionally
classifies the bad minimal models and the two free isomorphism classes;
the parallel packet proves a universal partial-Steiner avoidance bound and
the regular action on eight central conics. Those advances are inherited
here. Our independent finite checks remain useful controls, not priority claims.

Our anchor is the six-to-seven change in free deck quotients; the niche is the
simplex/Hamming pair of codes; the wildcard is the Paley orientation. The
five live concepts are branch incidence, point inertia, canonical characters,
Fano coding, and directed triangle systems. Every proposed passage below
specifies which of these it preserves and which needs an additional sidecar.

The canonical hostile is a Fano-dependent actual triple: matching a visible
three-line pattern can create fixed points instead of removing them. The
corrected near miss is the existing warning that the seven-cell hourglass
support is a mitre, not a Fano plane. Our Fano plane is an additional binary
labelling of the *seven branch lines*, not a relabelling of that additive
support. The least-used sidecar here is the full set of concurrent triples;
the number of nodes alone will fail below.

## 1. The exact freeness test

**Proposition (PROVED).** A subgroup H of G_n acts freely on Y_n if and only
if it intersects every point-inertia subgroup trivially.

If these intersections vanish, H acts freely on X_n, and also on its
equivariant resolution: a fixed point upstairs would map to one downstairs.
Conversely a nonidentity stabilizer at a smooth point already fixes that
point. At a node it acts on the exceptional P1 of the minimal resolution.
An automorphism of P1 of finite order has a fixed point over C, including
when its action on P1 is the identity. Thus a nontrivial intersection also
prevents freeness on Y_n. This uses the A1 and characteristic-zero hypotheses
of T1; it is not a rule for arbitrary singularities or characteristic two.

Every pair of distinct branch lines meets. Consequently a free H must avoid
each g_i and each g_i+g_j, whether the corresponding intersection is double
or triple. Equivalently the images of all g_i in G_n/H are nonzero and
pairwise distinct. Actual triple crossings supply the remaining conditions.

## 2. Six roots: exactly the admissible three-plus-three partitions

Represent a deck element by a subset of six indices, modulo complementation.
Choose a representative of size at most three. Sizes one and two always have
fixed points. A size-three representative is free precisely when neither
it nor its complementary triple is concurrent. Thus free involutions are
exactly the admissible unordered partitions of the six labels into two triples.

Moreover **every free deck subgroup has order at most two**. Because six is
even, subset-size parity descends to G_6. Each possible free nonidentity
element has odd parity. The sum of two distinct such elements is nonzero
and even, hence has a representative of size two and has fixed points.

A free involution on Y_6 gives an etale degree-two quotient with chi(O)=1,
q=0 and numerically trivial canonical class. Its canonical class is
nontrivial two-torsion: triviality would give p_g=1, contradicting chi=1
and q=0. This is an Enriques surface. The K3-to-Enriques construction is
classical, for example [Beauville's primary paper](https://arxiv.org/abs/0902.3721).
The incidence criterion above identifies which deck involutions realize it.

For the standard six labels E,E±x,E±z,E+x+z, the two actual triples exclude
two of the ten partitions. There are eight free deck involutions. The eight
nodes of X_6 form four orbits, so the quotient nodal model has four A1 nodes.
These are counts of embedded subgroups, not of nonisomorphic Enriques surfaces.

## 3. Seven roots: exactly the admissible Fano kernels

Since G_7 has dimension six, the nonzero pairwise-distinct image condition
requires dim(G_7/H)>=3. Hence **a free deck subgroup has order at most eight**.
At equality the seven image labels are all seven nonzero vectors of F2^3.
Their sum is zero, so every such bijection defines a quotient map from G_7.
Postcomposition by GL3(F2), of order 168, leaves the kernel unchanged.
The 7!/168=30 kernels are therefore exactly the 30 labelled Fano structures.

For a labelling alpha_i, the preimage of H in F2^7 is the kernel of the
3-by-7 matrix with columns alpha_i: the Hamming [7,4,3] code. It contains the
all-ones word. Its seven nonzero classes modulo that word are exactly the
seven Fano triples, each identified with its complementary quadruple.
The dual character inclusion is the simplex [7,3,4] code. These have different
typed roles, but there is an explicit identification: since seven is odd,
the canonical even-representative map [w] -> w+(sum w)1 sends G_7 to the
even-weight subspace of F2^7. It sends Hamming/<1> to the constant-weight-four
simplex subspace itself. Session60, section 2, proves this self-orthogonal
identification; it does not discard the distinction between a deck action
and its dual character inclusion.

Distinct nonzero labels automatically handle single and double crossings.
At a triple, the inertia maps injectively if and only if the three labels
are independent, equivalently their xor is nonzero. Therefore:

> All order-eight free deck groups are precisely the Fano kernels whose seven
> Fano lines avoid every actual concurrent triple of branch lines.

This is a complete classification of maximal-degree *free deck quotients*,
not a classification of every automorphism or every quotient of Y_7.

The incoming Fano-code packet, section 4, independently exhausts all 5,596
partial Steiner triple systems on seven labels and finds at least eight
Fano structures disjoint from each. Actual line-concurrence triples form
such a system. Together with the bound above, this proves that every
seven-root arrangement covered by T1 has maximal free deck degree eight.
The bound of eight *choices* is sharp in that enlarged combinatorial
universe; no characteristic-zero realization of its extreme is asserted.

For a free H, Y_7/H is smooth and minimal of general type. Etale division gives
K^2=16/8=2 and chi(O)=8/8=1. Since H^1(O) injects into that of Y_7, q=0 and
p_g=0. It is a numerical Campedelli surface, in the classical seven-line
(Z/2)^3 construction. Each actual triple gives one A1 point on X_7/H and
one (-2) curve on its resolution. For the hourglass, the 24 original nodes
form three H-orbits and yield three such curves.

The singularity boundary is classical and load-bearing: an independent
triple gives one A1, while a dependent triple gives **two 1/4(1,1)** points.
See [Alexeev–Pardini, sections 3.1 and 3.4, Lemma 3.4](https://arxiv.org/html/0901.4431v4).
Locally the offending all-three-sign flip at an A1 node lifts to scalar i
on its C2 double cover; it fixes the exceptional P1 after resolution.
Such a kernel is not free, even though it still has no invariant canonical
form. The inherited hourglass classification below uses additional central
component incidence and an explicit genus-one pencil; a singularity count
or canonical square alone does not give it.

## 4. The exact hourglass count is sixteen

Order the seven branch forms as

    0 E, 1 E+x, 2 E-x, 3 E+z, 4 E-z, 5 E+x+z, 6 E-x-z.

The actual triples are 012,034,056. The labels

    (1,2,4,3,6,5,7)

make all three independent and give a free kernel. The labels

    (1,2,3,4,5,6,7)

make all three dependent and give a hostile kernel.

There are 16 free kernels, 12 with exactly one bad triple, and two with
three bad triples; none has exactly two. One can obtain the count without
exhausting label permutations. Fix the Fano label of the central branch.
The three Fano lines through it partition the other six labels into pairs.
Compare their induced perfect matching with the actual pairs {12,34,56}.
Among the 15 perfect matchings on six objects, respectively 8,6,1 have
0,1,3 common pairs with a fixed matching. Each matching is realized by
7*(3!)*(2^3)=336 labelings. Dividing by 168 gives 16,12,2 kernels.
Two shared pairs force the third, explaining the missing case.

The checker independently enumerates all 1,395 rank-three subgroups of F2^6
and applies the inertia test directly; the same sixteen survive. This avoids
using the Fano description as both implementation and validation oracle.

The number sixteen also occurs as K_Y7^2. These are two different counts:
the latter is complete-intersection degree, while the former is an incidence
avoidance count. The next census gives a direct hostile to identifying them.

### Inherited minimal-model and isomorphism classification

[Session60, section 4](../session60_six_seven_geometry_20261025/seven-lines/RESULT.md)
proves the following stronger result for this specific three-spoke hourglass.
Write d for the number of bad actual triples.

| d | Embedded kernels | Singularities of the normal quotient | K^2 of its minimal resolution | Minimal-model type |
|---:|---:|---|---:|---|
| 0 | 16 | 3 A1 | 2 | Campedelli |
| 1 | 12 | 2 A1 and two 1/4(1,1) | 0 | Properly elliptic |
| 3 | 2 | Six 1/4(1,1) | -4 | Enriques after four contractions |

For d=1, the pencil through the dependent crossing pulls back to a connected
genus-one fibration, with O(2K)=h*O_P1(1) and P_(2m)=m+1. Its resolution
is already minimal. For d=3, four disjoint central (-1)-curves contract to
an Enriques surface with K^2=0 and nontrivial two-torsion canonical class.
Thus the resolution's K^2=-4 must not be assigned to that minimal model.
These are nonfree seven-root quotients, distinct from the free six-root
Enriques quotients in section 2. This classification is not automatically
a classification of the other seven-subsets below.

The sixteen free kernels form two D12 orbits, of sizes four and twelve.
The intrinsic bicanonical map makes that orbit test complete: there are
exactly two complex Campedelli isomorphism classes, with automorphism groups
C2×A4 and (C2)^3 respectively. Counts of embedded kernels remain different
from counts of surface classes. The analogous properly elliptic isomorphism
classification and a complete singular/multiple-fiber list remain open.

## 5. All six- and seven-subsets of the nine magic lines

Append labels 7 E+x-z and 8 E-x+z. Exact determinants recover precisely

    012,034,056,078,157,268,358,467

as the eight concurrent triples, with no four concurrent. Exhausting all
labelled coordinate subsets gives the following refined counts.

| Six-root triple count | Free involutions | Number of six-subsets |
|---:|---:|---:|
| 0 | 10 | 2 |
| 1 | 9 | 20 |
| 2 | 8 | 40 |
| 2 | 9 | 6 |
| 3 | 7 | 16 |

Two complementary concurrent triples forbid one partition rather than two;
this explains the second row with triple count two. Every one of the 84
six-subsets has a maximal free deck quotient of degree two.

| Seven-root triple count | Free order-eight subgroups | Number of seven-subsets |
|---:|---:|---:|
| 2 | 20 | 4 |
| 3 | 16 | 10 |
| 3 | 17 | 8 |
| 4 | 14 | 12 |
| 5 | 12 | 2 |

Every one of the 36 seven-subsets has a maximal free deck quotient of degree
eight. Subsets (0,1,2,3,4,5,6) and (0,1,2,3,5,6,7) both have three triples,
24 nodes, p_g=7 and K^2=16, but respectively 16 and 17 free Fano kernels.
The missing coordinate is the incidence pattern, not another scalar genus.

## 6. A genuine Paley connection, with orientation and normalization retained

Two Fano structures F and F' on these seven labels are disjoint precisely
when their associated kernels have trivial intersection: their nonzero
kernel elements are represented by their Fano triples. Both kernels have
dimension three, so disjointness is equivalent to G_7=H direct-sum H'.

There are 120 unordered disjoint Fano pairs. The exact finite enumeration
also proves the following labelled statement. The Paley tournament on Z/7,
oriented i→j when j-i belongs to {1,2,4}, has fourteen directed triangles,
partitioned into two disjoint Fano systems. Across its 240 labelled copies,
every disjoint Fano pair occurs exactly twice, from reverse orientations.
Thus an unordered pair forgets one orientation bit. The orientation
automorphism group has order 21; the incidence group of one Fano plane has
order 168. This does not identify those groups or supply a G2 action.

Relative to the actual hourglass triples, the 120 disjoint pairs have:

| Bad actual triples in the two Fano systems | Pairs |
|---|---:|
| (0,0) | 32 |
| (0,1) | 48 |
| (0,3) | 16 |
| (1,1) | 24 |

Consequently there are 32 pairs of complementary *free* Campedelli quotient
maps, or 64 labelled Paley orientations lying over them. No choice of such
an orientation is supplied canonically by the complex surface.

Here the combinatorics has an actual geometric consequence. Write
K=C(P2), L=C(X_7), and L_i=L^(H_i) for a complementary pair. Galois theory
gives L_1 intersect L_2=K and L_1 L_2=L. Since the extensions are Galois,
their tensor product over K is a field of degree 64, equal to L. For
Z_i=X_7/H_i, the natural map

    X_7 -> Z_1 fiber-product-over-P2 Z_2

therefore identifies X_7 with the **normalization of the component dominating
P2**. In this setting the whole fiber product has that component only:
the covers Z_i→P2 are finite flat (their invariant algebras are direct
summands of the finite locally free root-cover algebra), so their tensor
algebra is torsion-free and injects into its generic field. Hence it is
integral and X_7 is its normalization. The normalization statement follows
because X_7 is normal, finite over the fiber product, and has its function field.
Here X_7 and Z_i are the finite normal cover models; replacing them by their
smooth resolutions is not part of this finite fiber-product assertion.

Normalization cannot be dropped. Near a generic branch line both quotient
covers have ramification index two. After separating the unramified local
branches, the raw product contains the local equation u^2=v^2. Its two
branches must be separated by normalization. This is a geometric recovery
of a degree-64 root cover from two degree-eight quotients, not a rule for
lifting rational or integral points from either quotient.

## 7. Arbitrary binary relabelling and the geometric half-turn are different maps

The [support-reflection packet](../reflection_bridge_20261025/seven.md)
studies J(x)=x xor 7 for x!=7, with J(7)=7. It preserves three of the
seven standard Fano lines and destroys four. Applied to the good labels
(1,2,4,3,6,5,7), it gives

    (6,5,3,4,1,2,7).

Now 6 xor 5 xor 3=0: only the first actual triple becomes dependent. Under
this marking J induces the branch-index permutation (0 4)(1 5)(2 3), fixing
6. That permutation does not preserve the actual concurrent triples. This
arbitrary relabelling changes a free Campedelli quotient into the d=1
properly elliptic case. The test refutes preservation of freeness by an
arbitrarily marked J; it does not refute a geometric symmetry.

The actual half-turn A:(E,x,z)->(E,-x,-z) instead induces
(1 2)(3 4)(5 6), fixing the central index 0. It literally permutes the
branch forms and preserves all three actual concurrent triples. Under the
incidence-adapted marking

    (7,1,6,2,5,4,3),

J is exactly this A. As proved in
[Session60's reflection theorem](../session60_six_seven_geometry_20261025/reflection/BOOLEAN_REFLECTION.md)
and the [parallel Fano packet, section 6](../six_seven_geometry_20261025/fano_codes/RESULT.md),
A lifts to the root surface and conjugates every Fano kernel to a different
one in the same d-stratum. The thirty kernels form 8,6,1 pairs within
d=0,1,3 respectively. The adapted example exchanges the two d=3 kernels;
the same geometric A pairs the sixteen free kernels among themselves.
It induces an isomorphism Y/H -> Y/A(H), not a descent through a fixed
marked kernel. Its cycle type 1·2^3 cannot be a Fano collineation: a
nonidentity involution of GL3(F2) fixes three nonzero vectors, not one.

The new standalone integration checker verifies the two markings, the
literal linear action on branch coefficients, all 35 incidence decisions,
and the kernel conjugation and stratum preservation for all thirty kernels.
The needed sidecars are the marking and actual concurrent triples; the
zero/nonzero xor status of each such triple decides freeness.

## 8. Candidate atlas and stopping boundaries

| Proposed bridge | Source, map and preserved predicate | Lost datum / result |
|---|---|---|
| Six to seven | Root sign groups and point inertia; admissible partitions become admissible Fano kernels | PROVED maximal free deck degree 2→8, producing Enriques and Campedelli quotients; no automatic identification of their rational points |
| Simplex / cube | Simplex [7,3,4] character code embeds in the six-dimensional projective sign cube | PROVED dual quotient language; simplex here is a binary code, not an equivalence of Euclidean shapes |
| Seven canonical forms / seven Fano points | Canonical characters are the seven even complements of singleton labels | Canonical characters have binary rank six and form one seven-circuit, whereas Fano labels have rank three; cardinality alone gives no linear isomorphism |
| Canonical forms / Fano kernel | Restrict canonical characters to any of the thirty kernels | FINITE-EXACT regular representation minus the trivial character: no invariants and trace -1 on nonidentity elements; this does not detect freeness, since bad kernels have the same representation |
| Paley / paired Campedelli covers | Disjoint directed-triangle Fano systems give complementary deck kernels | FINITE-EXACT 32 both-free pairs; PROVED recovery after normalization; orientation is an extra bit |
| Hourglass mitre / Fano | Put a separate binary Fano labelling on branch lines | The old additive incidence is not preserved; the admissibility calculation succeeds precisely because zero-sum actual triples are avoided |
| LRC rank seven / seven-root cover | Both permit incidence and circuit encodings | OPEN connection: no map presently preserves LRC's phase, owner, or overlap invoice; matching the integer seven does not supply one |
| Six/seven / unrestricted tournaments | Seven admits Paley's two Fano systems, while the six-boundary port computation uses Hamiltonian reachability | The tested bridge preserves cyclic triangles, not port vertices, exact clone weights, or black-core reachability; no unrestricted A000568 recurrence or all-order connectivity consequence |
| Arbitrarily marked complement / quotient | Apply J to the chosen binary labels | REFUTED general preservation of freeness: the displayed marking moves d=0 to d=1; retain actual triple xor data |
| Geometric half-turn / quotient | Permute opposite branch pairs and conjugate the kernel | PROVED isomorphism between paired quotients within each d-stratum; no fixed Fano kernel and no arbitrary-marking identification |

The canonical-character calculation belongs to the companion arithmetic
lane; the checker here is an independent finite restriction audit. None of
these results proves Collatz convergence, an RH criterion, P versus NP, or
an implication among those questions. A future connection must specify a
map preserving their actual decision predicate, including arithmetic or
phase data erased by these finite geometric quotients.

## Reproduction and scope

From the repository root:

```sh
python3 -B 04-computation/six_seven_gale_cube_20261025/check_free_deck_codes.py --output /tmp/six_seven_free_deck.normal.json
python3 -B -O 04-computation/six_seven_gale_cube_20261025/check_free_deck_codes.py --output /tmp/six_seven_free_deck.optimized.json
cmp /tmp/six_seven_free_deck.normal.json /tmp/six_seven_free_deck.optimized.json
cmp /tmp/six_seven_free_deck.normal.json 05-knowledge/results/six_seven_gale_cube_20261025/free_deck_codes.json
python3 -B 04-computation/six_seven_gale_cube_20261025/check_reflection_markings.py --output /tmp/six_seven_reflection.normal.json
python3 -B -O 04-computation/six_seven_gale_cube_20261025/check_reflection_markings.py --output /tmp/six_seven_reflection.optimized.json
cmp /tmp/six_seven_reflection.normal.json /tmp/six_seven_reflection.optimized.json
cmp /tmp/six_seven_reflection.normal.json 05-knowledge/results/six_seven_gale_cube_20261025/reflection_markings.json
```

The checker uses only Python's standard library and explicit exception-based
checks. Its declared universes are 5,040 binary label permutations, all
1,395 rank-three subgroups of F2^6, all 84 six-subsets and 36 seven-subsets
of the specified nine forms, and all 240 labelled Paley tournaments. It also
checks all order-four subgroups in each six-subset, the positive good labels,
the three-dependent-triple hostile, and the modified-complement hostile.
Determinants recover actual incidences independently of the Fano enumeration.

The byte-identical normal and optimized certificate has SHA-256
`44506708332e41227250d9beab8a4afae9ab7d261874f08bdc756ebdc5f48f81`.
The standalone reflection integration audit passes **15,334 explicit gates**
in normal and optimized runs. Its byte-identical `reflection_markings.json`
has SHA-256
`c6b92ab849bec638c7a9d9bca164b09ca39bb9431328a3228c036634ac9d77b3`.
All subset counts are on labelled branch subsets and embedded subgroups.
The finite censuses themselves do not quotient by geometric isomorphism;
the two-class refinement above is separately inherited from Session60.
No exhaustive census of non-deck automorphisms of the root cover is claimed.
