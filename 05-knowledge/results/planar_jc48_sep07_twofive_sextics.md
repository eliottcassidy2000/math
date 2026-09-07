# One finite (2,5) cusp in the (4,6) polynomial-normalization class

Status: **PROVED + FINITE-EXACT; two independent analytic/source audits PASS.** JC(2) remains
**OPEN**. No assertion of external novelty is made.

## Precise claim

Let C be an irreducible complex affine plane curve with a birational polynomial
normalization t -> (U(t),V(t)), deg U=4 and deg V=6. Suppose its finite
singularities are exactly one unibranch (2,5) cusp and N>=2 ordinary nodes.
Then C cannot be the whole nonproperness support of a noninvertible complex
Keller map A2 -> A2.

The claim is about this declared curve class, not all rational sextics, all
nonproperness curves, or the degree of a Keller map. Dependencies are the
[odd-cusp actual passport](planar_jc48_sep06_odd_cusp.md),
[actual infinity-(2,9) plumbing exclusion](planar_jc48_sep06_boundary_plumbing.md),
[Nori budget](planar_jc48_sep06_resolution_budget.md), and the independently
audited [two-meridian infinity-(2,11) representative](planar_jc48_sep07_next_braid.md).
The normal forms, infinity-(2,13) proof, and deformation argument below are new
session deductions. The [exact source](../../04-computation/planar_jc48_sep07_twofive_sextics.py)
checks algebraic coefficients and the entire A6 boundary representation universe;
it does not replace the geometric deformation proof.

## Inheritance and decisive missing coordinate

The closest proved mechanisms are the (2,9) marked plumbing obstruction and the
two positive affine meridians for one (2,11) curve. A canonical hostile example
is the transitive A6 representation of that latter curve's infinity group:
infinity relations alone do not decide this case. The corrected near miss is
Nori equality at four nodes, which proves no exclusion. The least-used sidecar
is simultaneous embedded resolution with the actual infinity line marked. It
transports the successful affine-meridian certificate across the good family.

## Exhaustive normal form and four infinity types

Translate the finite cusp parameter and image to zero, make U monic, and
subtract the coefficient of t4 in V times U from V. Both linear coefficients
vanish, so

    U=t4+a t3+b t2,
    V=c t6+e t5+f t3+g t2,                 c != 0.

If b != 0, cancelling the t2 term and demanding order at least four gives
f=ga/b. Put d=g/b. The finite cusp is precisely (2,5) when

    e+2ad/b != 0.                                      (1)

Indeed V-dU has leading -d t4, and adding (d/b2)U2 cancels it; the fifth
coefficient is (1). If d=0 the same test applies directly. The use of U2 here
is a local singularity test, not an additional degree-preserving normal-form
operation.

If b=0 and a!=0, multiplicity two requires g!=0 and the orders three and two
give a (2,3) cusp, contradicting the hypothesis. If a=b=0, U=t4 and the
(2,5) condition is exactly gf!=0: V2-g2U has fifth coefficient 2gf.

At the unique point at infinity set z=1/t, X=U/V and Z=1/V. Their orders are
two and six; the line at infinity is Z=0 and has contact six. The seventh
coefficient of Z-c2 X3 is (2e-3ca)/c2. If this vanishes in the b!=0 case,
a cannot be zero: otherwise e=0 and both polynomials belong to C[t2],
contradicting birationality. Rescale t=a tau, U by a4 and V by ca6/2 to obtain

    U=t4+t3+beta t2,
    V=2t6+3t5+delta t3+delta beta t2,       beta != 0.    (2)

These affine coordinate rescalings preserve positive meridians and the whole
support question. Put

    rho=Z/X3-4,       k4=-6-24 beta,
    tau=rho/X-k4.

The coefficient of z in tau is 7+12 beta+8 delta. If it vanishes, then

    delta=-(7+12 beta)/8,                                (3)
    k5=2(3 beta+1)(20 beta+3),   eta=tau/X-k5,
    [z]eta=(3-20 beta)/4.

In this family (1) equals -7/(4 beta), so the finite cusp persists. At the
last exceptional value beta=3/20, k6=-7069/250 and

    [z](eta/X-k6)=3/100,
    [z13](Z-4X3-k4X4-k5X5-k6X6)=3/6400 != 0.

Thus the infinity branch has type (2,m), m in {7,9,11,13}. In the exceptional
U=t4 case, e!=0 gives m=7; e=0 gives ninth coefficient 2f/c2 !=0, hence m=9.
The arithmetic genus of an irreducible sextic is ten, so

| infinity m | finite nodes N=(17-m)/2 | Nori margin |
| --- | --- | --- |
| 7 | 5 | 2 |
| 9 | 4 | 0 |
| 11 | 3 | -2 |
| 13 | 2 | -4 |

This also bounds m<=13 independently under the stated N>=2 hypothesis.
The marked resolution has the original infinity line leaving the curve after
the third blowup: there rho=-c2 is nonzero. Subsequent even coefficients are
removed by local coordinate translations. A vanishing k4 or k5 does not make
a free center satellite: the other exceptional component is at infinity in
the transverse chart. There are then free double-point centers and the final
free/satellite pair. Thus the actual marked tree depends only on the listed m
within this normal form, not on nonvanishing auxiliary even coefficients.

For m=7 the already audited Nori inequality is strict. For m=9 the actual
marked tree and arrow are those of the existing plumbing theorem, whose
single-three-cycle A6 representations are all intransitive. This settles
those two cases.

## The infinity-(2,13) marked group

Starting with L of square one, blow up L, then L intersect E1, then L intersect
E2, then successively a free point on E3,E4,E5,E6, and finally E6 intersect E7.
The strict curve meets E8. The tree is

    L--E3--E2--E1,       E3--E4--E5--E6--E8--E7,

with weights in the order L,E1,...,E8

    (-2,-2,-2,-2,-2,-2,-3,-2,-1).

The negative intersection matrix has determinant -1. Its marked meridian
multiples in H1, for the positive curve meridian mu at E8, are

    (-6,-4,-8,-12,-10,-8,-6,-5,-10).

Use fibre generators l,a,b,c,d,g,h,e,f in that vertex order. The plumbing
relations, besides commutation across edges and [f,mu]=1, are

    l2=c, a2=b, b2=ac, c2=lbd, d2=cg,
    g2=dh, h3=gf, e2=f, f=he mu.                         (4)

The required odd-cusp passport has degree six, generic retained count three,
and monodromy A6 with mu a single three-cycle. We show every representation
of (4) with that meridian fixes a label. Relations eliminate to

    c=l2=a3, b=a2, d=(lb)^(-1)c2,
    g=c^(-1)d2, h=c^(-2)d3, f=c^(-5)d7=e2,
    mu=(he)^(-1)f.                                     (5)

The common square/cube c in A6 is of type identity, double transposition, or
five-cycle. If c is a double transposition, a=c, b=1, l has order four and
d=l^(-1); hence f=l3 has order four and cannot be a square in A6.

If c=1, then f=d7. The possible orders of d are 1,2,3,4,5. Orders 1,2,4,5
respectively make mu have order at most two, order four, f an impossible
square, or mu an element of order five. For order three, e=d2 and mu=d2;
therefore d is a single three-cycle. The subgroup generated by l and a is a
quotient of the spherical triangle group (2,3,3), hence of A4. Its permutation
image contains a single three-cycle, so it moves at most four of six labels:
a transitive degree-six A4 action is on cosets of an order-two subgroup and
has its three-cycles of type (3,3); a degree-four orbit together with a
nontrivial degree-two orbit is impossible since A4 has no quotient C2; and
two nontrivial degree-three orbits would again give type (3,3). All remaining
fibres and mu are in this same subgroup.

If c is a five-cycle, then l=c3, a=c2, b=c4, d=1, h=c3 and f=1. Now e is an
even involution and he mu=1. Let j be the label fixed by h. If mu moves j,
it contains j and two labels on the five-cycle. Cutting that five-cycle at
those two labels shows the product mu h has exactly two cycles on all six
labels. It therefore cannot equal the even involution e^(-1), which has four
or six cycles. Thus mu fixes j, and so do e and every fibre.

The exact independent control enumerates all360 elements of A6, all common
square/cube choices and every square root in (5), and then checks every
original relation (4). Exactly1120 assignments survive:40 fix three labels,
360 fix two, and720 fix one. This is a boundary representation census, not
a census of Keller maps. The analytic argument above explains the failure.
The actual infinity-complement epimorphism onto the affine complement then
contradicts the required A6 monodromy.

## Connected equisingular transfer of the infinity-(2,11) certificate

All remaining m=11 curves have the form (2),(3), with
beta in C minus {0,3/20}. Let B be the subset for which the normalization is
birational and the only finite singularities are the prescribed (2,5) cusp
and three ordinary nodes. It contains beta=1, the independently audited
representative U=t4+t3+t2, V=(16t6+24t5-19t3-19t2)/8.

B is Zariski open. Here are the needed details of this assertion. First remove
the parameters where a second normalization point maps to the finite cusp.
This is a closed parameter locus: near t=0, U=t^2(beta+t+t^2) with beta a unit,
so no additional zero can approach the fixed section t=0; t=infinity maps
to infinity, not to the finite cusp. The extra-preimage incidence is therefore
closed in the proper family P1 times the base, and its parameter image is
closed. This removal precedes lifting the normalization through blowups.

Now resolve the fixed finite cusp and infinity branch by their fixed sequences
of blowups, using the algebraic even-coefficient translations above. Their
centers are sections; excluding beta=0,3/20 keeps the final odd coefficients
nonzero. At each successive cusp center, the pullback ideal on the parameter
normalization is supported on t=0 with its fixed order and a unit, so it is
invertible. The same assertion holds at infinity with z=1/t. Thus the
normalization really lifts to the resolved surface family; a rational lift
over an unremoved extra-preimage locus is not being assumed. Remove next the
closed parameter locus where this proper normalization map fails to be
immersive. The resulting map is unramified and separated: it is the identity
on the base and has injective differential on every relative tangent line.
Its diagonal in the relative fibre product is therefore both open and closed.
The off-diagonal double-point incidence is closed in that fibre product and
proper over the base. In the threefold fibre product, the locus of pairwise
distinct points is likewise closed, by the three open-and-closed diagonal
conditions, and proper over the base. Nontransversality is a closed algebraic
condition on a double pair's two tangent lines. Hence its parameter image,
and the triple-incidence parameter image, are Zariski closed by properness.
Removing them leaves precisely ordinary double intersections away
from the resolved cusp and infinity. A nonbirational parametrization has a
positive-dimensional coincidence locus and cannot lie in this transverse
locus. Conversely every curve satisfying the declared good conditions lies
in this locus. The genus formula fixes the number of double points at three.
This algebraic incidence argument does not assume an auxiliary pair-elimination denominator
is a unit. In particular beta=1/4 must not be discarded just because one such
denominator vanishes. Thus B is nonempty and is the complex line with finitely
many points removed, and is path connected.

Over B, the three node loci form a finite etale multisection. Blow up this
whole multisection; no global numbering of nodes is required. Together with
the fixed cusp/infinity resolution, this gives a smooth proper family of
compact surfaces whose reduced total transform of C union L is a relative
simple normal-crossing divisor. Every stratum projects smoothly to B.
Such a family of pairs is locally differentiably trivial: in relative
normal-crossing charts lift a real vector field on B while fixing the divisor
coordinates, glue the lifts with a partition of unity, and integrate. The
lift still projects to the chosen base field and is tangent to every divisor
stratum. Properness gives the flow along compact subpaths. Its complex-normal
orientation preserves positive meridians up to conjugacy. All blowup centers
belong to the removed divisor, so the open fibre remains exactly A2 minus C.

Consequently every C in this good m=11 family has an affine-complement group
generated by two positive curve meridians, since the beta=1 representative
does. For a connected degree-d>1 unramified cover in which a generic curve
meridian fixes at least a retained sheets, put delta=d-a. A permutation moves
at most delta labels, so its disjoint nontrivial cycles can supply at most
delta-1 edges in a graph connecting its cycle orbits. Here a transitive
nontrivial action forces delta>=2; identity generators supply zero edges.
If r such meridians
generate a transitive action, the union graph is connected and hence

    r(delta-1) >= d-1.                                 (6)

For this whole-support odd-cusp passport d=6,a=3,delta=3. Two generators give
4>=5, a contradiction. The single infinity representation that motivated
this calculation remains valid; it simply does not extend through the
necessary affine braid relations.

## Scope, reproducibility, and next question

These independently audited four cases exhaust the stated (4,6), one
(2,5)-cusp, N>=2-node class. Equality in (6) can occur abstractly: two delta
cycles meeting at one label act transitively on 2delta-1 labels. Disjoint
three-cycles supply a hostile example to replacing transitivity by the
absence of a common fixed label. No argument here classifies multiple cusp
curves or higher-degree normalization pairs.

Run from the repository root:

    python3 04-computation/planar_jc48_sep07_twofive_sextics.py
    python3 -O 04-computation/planar_jc48_sep07_twofive_sextics.py

The always-active exact gates check coefficients, all four genus/Nori rows,
actual blowup weights and marked H1, and the full marked A6 census. Normal
and optimized outputs must agree with the frozen .out. The independently
certified two-loop paths remain in their separate immutable bundle.


## Independent audits and frozen pins

The [full analytic/source audit](planar_jc48_sep07_twofive_sextics_audit.md)
and the [separate marked-geometry audit](planar_jc48_sep07_twofive_sextics_geometry_audit.md)
both pass. The first repairs construction order by removing extra cusp
preimages before lifting the normalization, and proves algebraic openness
using the open-and-closed diagonal of the unramified separated map. The
second independently reconstructs the marked census and verifies that the
infinity13 case is nonvacuous. The beta=1/4 good-family control rejects an
auxiliary-denominator exclusion. These repairs were incorporated before
promotion; no previously proved theorem is retracted.

All3385 always-active gates pass normally and under optimization, including
both independent referee replays, with identical1890-byte output.

    source SHA256 02ef96f5a39fee7382fa81d0970fd06285c7bb3e180362bff234d4e2abad1e75
    output SHA256 2713e10ba1715d4ae3b9d5d0095be91fa901ea3f63a16c9d53778e15b072e472
    assignment SHA256 b35120fec36b3ded751aec092972d2e870466fb4d7defe0a193d41cc8412373f
