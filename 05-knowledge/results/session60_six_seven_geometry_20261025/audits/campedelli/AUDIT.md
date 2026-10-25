# Independent audit: Fano kernels and the Campedelli/elliptic/Enriques trichotomy

**Verdict: PASS over C.** No mathematical counterexample or unresolved proof gate
remains for the statements below. The geometric proof here is independent of
quotient-divisor discrepancy calculations: it resolves the cover after blowing
up the base, and computes the canonical bundle on the resulting smooth cover.
The exact finite audit enumerates all 1,395 rank-three subgroups of C2^6, not
permutations of the seven labels. This is not a historical novelty claim.

The scope is precisely the seven distinct complex lines

    L0=E, L1=E+x, L2=E−x, L3=E+z, L4=E−z,
    L5=E+x+z, L6=E−x−z.

Do not transport the geometric trichotomy unchanged to characteristic 2 or 3.
Do not interpret the 30 kernels as 30 pairwise nonisomorphic abstract surfaces.
The subsequently audited d=0 isomorphism refinement is in Section 8.
They are 30 covering kernels over the fixed numbered plane arrangement.

## 1. Exact combinatorics and free-action gate

Write D=F2^7/<(1,...,1)> and let e_i be the seven inertia generators.
A rank-three quotient D/K -> F2^3 has distinct nonzero line labels precisely
when no e_i and no e_i+e_j (i != j) belongs to K. There are exactly 30 such
K. Equivalently, choose a bijection from the seven branch lines to F2^3\{0}
and quotient by GL(3,2), giving 7!/168=30.

For a label matrix A with columns all seven nonzero vectors, its binary kernel
is the [7,4,3] Hamming code. It contains the all-one vector. Each projective sign
class has a unique even-weight representative, and the even subcode is the
[7,3,4] simplex code, equal to the row space of A: each nonzero row evaluation
has weight four; the row space is self-orthogonal; both relevant spaces have
dimension three. The seven complements of its nonzero codeword supports are
the seven Fano lines, i.e. the dependent triples of labels.

The only geometric triple intersections are

    (0,1,2), (0,3,4), (0,5,6),

at [E:x:z]=[0:0:1], [0:1:0], [0:1:−1]. There are twelve ordinary double
intersections and no fourfold point. This follows either by direct substitution
or by the exact cross-product computation in the independent checker.

Fixing the Fano point 0, its three incident Fano lines define a perfect matching
of the other six labels. There are 15 matchings. Each matching has exactly two
Fano completions: after binary indexing of the three pairs, take the four
transversals of even or odd parity. The number of common edges with the actual
matching {(1,2),(3,4),(5,6)} is 0, 1, or 3, with counts 8, 6, 1. It cannot be 2,
because the last pair is then forced. Therefore the number d of geometrically
concurrent Fano lines has distribution

    d=0: 16 kernels; d=1: 12 kernels; d=3: 2 kernels.

Let X be the seven-root diagonal complete intersection in P6. A nonidentity
simplex sign has eigenspaces supported on a four-set and its complementary
three-set T. A projective fixed point must have all coordinates in one of those
eigenspaces. Support in the three-set would force four roots to vanish, which
is impossible. Support in the four-set means exactly the three branch lines
indexed by T vanish. It occurs if and only if T is one of the three actual
triple intersections. Thus the K action is free on the entire normal X exactly
when d=0, not merely free on the smooth open part. Consequently it is free on
the minimal resolution as well. For d>0 it fixes nodes; ignoring those nodes
would produce a false freeness claim.

Each universal triple fiber has eight A1 nodes, with local inertia rank three.
If its label triple is independent, K has trivial local stabilizer and the eight
nodes form one orbit. If dependent, the stabilizer has order two and the eight
nodes form two orbits. The independent checker verifies these coset counts for
all 30 kernels, and verifies the fixed-locus statement directly from coordinate
supports independently of the dependency test.

## 2. Smooth cover after one blowup at every triple point

Let b:Y -> P2 blow up the three triple centers. Write H=b*O(1) and E_i for the
exceptional curves, with H^2=1, E_i^2=−1, H.E_i=E_i.E_j=0 (i != j). Let G=F2^3
and a_j be the seven distinct nonzero labels. For triple i put s_i equal to the
sum of its three labels. Normalize the pulled-back G cover to obtain f:S -> Y.

The strict transforms of the original branch lines total

    D_strict ∼ 7H−3(E_1+E_2+E_3).

The exceptional E_i is branch precisely when s_i != 0. This follows from
pulling the defining root products across a local blowup: the exceptional
monodromy is the sum of the three old monodromies. Every branch component has
inertia order two. A new crossing has labels s_i and one of the original three
labels. If s_i != 0 those labels are distinct, since equality would force the
other two original labels to coincide. They are therefore independent.
Original ordinary double crossings also have independent labels. The branch
divisor is simple normal crossings, and its two local inertia labels are
independent at every crossing. The normalized cover S is therefore smooth:
locally the nontrivial connected part is (u,v)->(u^2,v^2), or just one square
root, with disjoint etale copies accounting for the unused group factors.
No further base blowup or hidden resolution is required.

The exact branch divisor is

    D ∼ 7H−2 sum_(independent i) E_i−3 sum_(dependent i) E_i.

Hurwitz and f*D=2 R give the exact LINEAR equivalence

    2K_S ∼ f*(2K_Y+D) = f*(H−sum_(dependent i) E_i).       (A)

This is not merely numerical equivalence. In particular

    K_S^2 = 2(1−d).

### Local exceptional fibers and singularity types

For an independent triple, E_i is branch with inertia <s_i>. Its reduced
preimage is a connected degree-four cover of E_i=P1 with residual group
G/<s_i>, branched at three points. It is rational by Riemann-Hurwitz, and
f*E_i=2A_i gives A_i^2=−2. Its contraction is the one A1 singularity.

For a dependent triple, E_i is unbranched. Its preimage is two disjoint
connected degree-four covers of P1, each branched at three points and hence
rational. Each has self-intersection −4. These are exactly the exceptional
curves in the minimal resolutions of the two 1/4(1,1) singularities.
To see the analytic type, the universal triple germ is the quadratic cone
w^2=u^2+v^2, equivalently C2/{±1}. The stabilizer acts by simultaneous negation
of the three quadratic coordinates, lifted by scalar multiplication by i on
C2. Its quotient is C2/μ4 with weights (1,1), not another A1.

Thus S is the minimal resolution of the plane quotient: its exceptional curves
are precisely (3−d) disjoint (−2)-curves and 2d disjoint (−4)-curves. These local
curves are not asserted disjoint from the central ramification curves.

## 3. The holomorphic invariants for every d

For each nonzero character chi of G exactly four of the seven labels satisfy
chi(a_j)=1. Let t_i be the number of these four lines through triple center i.
The character line bundle on Y is

    L_chi = 2H−sum_i floor(t_i/2) E_i,
    f_*O_S = O_Y direct_sum_(chi != 0) O_Y(−L_chi).        (B)

This follows directly by dividing the pulled-back double-cover equation by the
square of the exceptional factor to the power floor(t_i/2). Since 0<=t_i<=3,
−L_chi has the form M=−2H+sum_(i in I) E_i.
For every subset I, H.M=−2 and H.(K_Y−M)=−1, so h0(M)=h2(M)=0 because H is nef.
Moreover

    chi(O_Y(M)) = 1 + M.(M−K_Y)/2 = 0.

Hence h1(M)=0 as well. Since Y is rational, (B) proves

    p_g(S)=q(S)=0

for all 30 kernels, without assuming the singular quotient is smooth or
appealing to deformation invariance. The independent program checks all eight
possible subsets I and the integer Riemann-Roch intersections.

## 4. d=0: minimal numerical Campedelli surfaces

Equation (A) becomes 2K_S ∼ f*H. It is nef and big, so K_S is nef and big and S
is minimal of general type. Together with K_S^2=2 and p_g=q=0, this is precisely
a numerical Campedelli surface. The contracted quotient has exactly three A1
nodes. Its canonical bundle is ample as a rational class because its double
is the finite pullback of O_P2(1); it is the canonical model. In this case the
free degree-eight quotient of the minimal resolution of X is S.

This audit does not independently assert that X is simply connected or use
that assertion to identify pi_1(S). The numerical Campedelli conclusion and
free quotient statement need no such additional input.

## 5. d=1: minimal properly elliptic surfaces, proved by a pencil

Let E_b be the exceptional curve over the unique dependent triple. The divisor
F_Y=H−E_b is the basepoint-free pencil of lines through that center, even after
blowing up the other two centers. Its pullback gives a morphism

    h:S -> P1,        2K_S ∼ h*O_P1(1).

A generic pencil line meets precisely the other four branch lines. Their labels
are the complement of a Fano line and span all of G, so the induced degree-eight
cover of that P1 is connected. Each of the four branch points has inertia two.
Riemann-Hurwitz gives 2g−2=−16+4*4=0. Thus the generic fiber is a smooth connected
genus-one curve. Stein factorization gives h_*O_S=O_P1, and consequently

    h0(S,2mK_S)=m+1  for every m>=0.

In particular kappa(S)=1 (exactly, not merely a guess from K^2=0). Since K_S is
nef, S has no (−1)-curve and is already minimal. It is properly elliptic.
No abundance theorem or Enriques exclusion by numerics is needed.

For compatibility with the worker's central-curve proof: the central strict
line B has class H−E_1−E_2−E_3 and B^2=−2. The two independent exceptional
curves have the same nonzero label sum s. Over B, the residual monodromy has
rank one, so there are two rational degree-two components R_1,R_2, each with
square −2. Each independent exceptional (−2)-curve meets each R_j once.
Therefore the central pencil fiber is twice the reduced four-cycle I4:

    h*(central point)=2(R_1+R_2+A_1+A_2).

This explains why the worker obtains 2K_S ∼ 2F with F the reduced I4 fiber;
there is no factor-of-two contradiction with (A).

## 6. d=3: four blowups of an Enriques surface

Now (A) says 2K_S ∼ f*B, where B=H−E_1−E_2−E_3 is the central strict line.
Its inertia is <a_0>. All exceptional curves are unbranched, and B is disjoint
from the other strict branch lines. Thus its residual degree-four cover is
unramified. Since B=P1, it splits into four copies R_1,...,R_4, pairwise disjoint,
each mapping isomorphically to B. These multiplicities are essential:

    f*B=2(R_1+...+R_4).

Projection formula gives 2R_j^2=B^2=−2, hence each R_j is a smooth rational
(−1)-curve. Let eta:S -> T contract the four disjoint curves. The canonical
blowup formula and (A) yield

    K_S = eta*K_T + sum R_j,
    eta*(2K_T) ∼ 0.

Pullback on Picard groups is injective for these blowups, so 2K_T∼0 LINEARLY.
Blowing down preserves p_g=q=0. The line bundle K_T is not trivial, because a
trivial canonical bundle would have h0(K_T)=1. Consequently T is an Enriques
surface in the classical definition; K_T is nef and T is minimal. Exactly four
contractions have been performed, consistent with K_S^2=−4 and K_T^2=0.

The six exceptional (−4)-curves on S meet these four central curves in the six
edges of K4. Indeed central components are G/<a_0>; each dependent triple
<a_0,b> supplies the two cosets of that subgroup, each joining vertices [g]
and [g+b]. The three pairs exhaust the three nonzero directions of G/<a_0>.
After contracting the four central curves each old (−4)-curve becomes a
(−2)-curve, as an independent consistency check.

## 7. Sources, scope, and reproducibility

The proofs above are explicit calculations for this arrangement. The following
primary references verify the standard geometric input and terminology:

* V. Kharlamov and Vik. Kulikov, *Surfaces with DIF != DEF real structures*,
  Section 1.4, Lemmas 1.4--1.7 and Theorem 1.6, gives the blowup construction,
  exceptional monodromies, and the −2 versus two −4 exceptional curves.
  Section 1.5, Theorem 1.10 gives the Campedelli characterization of admissibly
  labeled seven-line covers. https://webdoc.sub.gwdg.de/ebook/serien/e/mpi_mathematik/2005/76.pdf
* R. Pardini, *Abelian covers of algebraic varieties*, J. reine angew. Math.
  417 (1991), 191--213/214, supplies the general abelian-cover framework.
  https://doi.org/10.1515/crll.1991.417.191
  The eigensheaf correction and local smoothness needed here are proved above,
  rather than relying on an inaccessible theorem statement.
* I. Dolgachev, *A brief introduction to Enriques surfaces*, Section 3, gives
  the classical definition q=p_g=0 and 2K=0. The implication used here is that
  definition, not a classification inferred from K^2 alone.
  https://sites.lsa.umich.edu/idolga/wp-content/uploads/sites/1334/2024/08/Kyoto13.pdf
  https://arxiv.org/abs/1412.7744

Reproduction:

    python3 independent_check.py
    python3 -O independent_check.py

Both completed successfully and generated byte-identical certificates. The
script uses explicit RuntimeError checks, so optimization cannot remove tests.
Its universe, exclusions, finite counts, and exact arithmetic are recorded in
independent_certificate.json. Its computations supplement, not replace, the
proofs of smoothness, canonical linear equivalence, and minimality above.

**First-checkpoint gate:** base arrangement, 30-kernel census, exact free locus,
and the 16 numerical Campedelli quotients pass independently.
**Extended-classification gate:** the 12 minimal properly elliptic resolutions
and two four-point blowups of Enriques surfaces also pass, via (A)--(B).
**Remaining blocker:** none for those scoped assertions. Global abstract
isomorphism-orbit counts, complete automorphism groups, fundamental groups,
and arithmetic descent are not established by Sections 1--7. The scoped d=0
isomorphism/automorphism upgrade in Section 8 is established separately.


## 8. Additional audited upgrade: exactly two Campedelli isomorphism classes

This section strengthens only the d=0 conclusion. The other d values still
have cover-orbit counts, not asserted abstract isomorphism classifications.

The arrangement automorphism group is exactly D12 of order twelve. Indeed its
central line is intrinsically the unique line containing three triple points,
so an automorphism is affine on E=1. It permutes the three parallel pairs and
therefore their three midpoint lines. The midpoint lines meet uniquely at the
origin, forcing the translation part to vanish. The resulting dual linear
map permutes {±x,±z,±(x+z)}. For the image of x there are six choices, and for
the image of z there are exactly two choices with their sum in that set;
these yield twelve matrices, all realized. This verifies the completeness,
not just the existence, of the symmetry worker's D12 subgroup.

For d=0 let q:S -> P2 be the resolved degree-eight plane cover. Since the nodes
are rational and the covering characters all select four lines,

    q_*O_S = O_P2 direct_sum O_P2(−2)^7,
    O_S(2K_S) = q*O_P2(1).

Projection formula therefore gives h0(S,2K_S)=3, and the COMPLETE bicanonical
map is exactly q. Every abstract isomorphism of these smooth minimal surfaces
must preserve the bicanonical linear system and descend to an element of
PGL3 taking their branch arrangements to each other. For the fixed arrangement
it descends to D12. Conversely, each arrangement symmetry taking one kernel
to another lifts by literally permuting the seven root coordinates and thus
induces a quotient isomorphism. Consequently two d=0 surfaces are isomorphic
if and only if their kernels are in the same D12 orbit.

The independent program has now been extended to apply the two realized
permutations to all 30 independently obtained kernels. It confirms orbit
sizes 4,12 for d=0, sizes 6,6 for d=1, and size 2 for d=3. Thus the sixteen
Campedelli kernels give EXACTLY TWO complex surface isomorphism classes.

For each d=0 kernel the descent sequence is

    1 -> F2^3 -> Aut(S) -> Stab_D12(K) -> 1.

The kernel consists exactly of deck automorphisms of the degree-eight function
field extension. Each stabilizer element lifts by a literal root-coordinate
permutation. These lifts form a complement, so the sequence splits. For the
orbit of size 12 the stabilizer is trivial, giving Aut(S)=C2^3 of order eight.
For the orbit of size four it is C3. Its nonidentity elements act nontrivially
on the quotient deck space, fixing exactly two of its eight elements (also
checked independently). The F2-representation splits as a one-dimensional
trivial module and the irreducible two-dimensional C3-module. Consequently

    Aut(S) = C2 × (C2^2 semidirect C3) = C2 × A4,
    |Aut(S)|=24

for the four-element orbit. No nonsplit extension ambiguity remains.

This bicanonical descent argument is self-contained from already audited
formulas. As a primary-source consistency check, Kharlamov--Kulikov Section
1.6, Lemma 1.11 and Theorem 1.13 assert the corresponding bicanonical and
descent statements for Campedelli arrangements, including nodal ones.

**Extended scope status:** the two d=0 isomorphism classes and their complete
automorphism groups pass. No abstract isomorphism classification for d=1 or
d=3, no fundamental-group calculation, and no arithmetic descent is claimed.
