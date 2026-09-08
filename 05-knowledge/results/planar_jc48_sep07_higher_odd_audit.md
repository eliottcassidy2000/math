# Independent audit: exhaustive higher odd cusps and connected actual families

**Status: INDEPENDENT ANALYTIC + EXACT SOURCE AUDIT PASS for classification
and transport.** The representative braid certificates are a separate
required dependency for the whole-support exclusion in the
[primary note](planar_jc48_sep07_higher_odd.md). No numerical word is used
in this audit. JC(2) remains OPEN.

**Dependency closure, September 8:** the separate
[higher-braid audit](planar_jc48_sep07_higher_braid_audit.md) now accepts all
three actual representatives, including complete normal and optimized
replays. The previously conditional whole-support consumer is therefore
paid. No source or output bytes in this classification bundle changed.

I read the full proof and standalone source. The normal form preserves the
polynomial degree pair by parameter translation and affine target changes.
For b!=0, cancellation of the third and fifth odd terms forces f=ga/b and
e=-2ad/b. If a=0 both coordinates become even; birationality therefore
allows the normalized a=1,c=2 chart. I independently expanded the sixth
coefficient of V-dU+(d/b2)U2 as2+d(1+2b)/b2. Removing it by a multiple of
U3 gives seventh coefficient -(6b2+(4b+3)d)/b3. On its zero locus,
d=-6b2/(4b+3), the next even removal gives ninth numerator
4(1+6b)-24(b+2)=-44. At b=-3/4 the seventh is instead8, so dividing by
4b+3 loses no higher cusp. Thus the finite cusp is exactly seven or nine,
never eleven or higher.

The independent infinity formula is (2e-3ca)/c2. In the normalized chart
it becomes -(2d+3b)/(2b). Its zero gives d=-3b/2, where the finite seventh
is9/(2b2) and the infinity ninth is7/16. The finite-nine locus has nonzero
infinity seventh -9/[2(4b+3)]. Hence the only types are(7,7),(7,9),(9,7).
The literal source verifies that every lower term vanishes before the
claimed characteristic exponent; these are singularity tests, not
unjustified global polynomial normalizations. Rational sextic genus ten
then gives node counts4,3,3 within the declared inventory.

For b=0,a!=0, multiplicity two yields a(2,3) cusp, outside this domain.
For a=b=0, multiplicity two requires g!=0; absence of the fifth term
forces f=0, and birationality forces e!=0. The local difference U-V2/g2
has first nonzero term -(2e/g)t7, and the infinity seventh is e/2.
Thus the exceptional U=t4 chart belongs exactly to the(7,7) stratum.

I particularly checked the irreducible family used to retain that chart.
The defining polynomial2ag+eb2 is primitive and degree one in e over the
fraction field of C[a,b,g], so is irreducible by Gauss. Where b or g is
nonzero, one of its partial derivatives b2 or2g is nonzero. Its b-chart
and g-chart are smooth affine coordinate charts with one coordinate
inverted. The g-chart solves a=-eb2/(2g) and has actual polynomial V with
cubic coefficient -eb/2. The finite seventh from V as the order-two
coordinate agrees on overlap with -(b/g) times the b-chart coefficient.
The two nonvanishing conditions therefore glue. At b=0 it specializes
to the exceptional coefficients without a singular or missing chart.

The local cusp/infinity conditions define open subsets of this smooth
irreducible family. Intrinsic good-node conditions define a further
Zariski open subset. As in the previously audited twofive theorem, extra
cusp preimages are excluded first; either the quadratic U or V coefficient
is a unit near the fixed cusp section. The successive centre ideals then
pull back to parameter powers times units, allowing actual normalization
lifting. After removing nonimmersion, the unramified separated map has
open-and-closed diagonal. Nontransverse off-diagonal pairs and pairwise
distinct triples have proper algebraic incidence, so their parameter
images are closed. This retains exactly the good locus, including all
good exceptional curves and possible repeated projection values.

Every nonempty good chart is path connected: the complex affine line
between two good points is not contained in the bad algebraic set and
meets it in finitely many points, around which a path can pass. Both
charts, when nonempty, meet in the irreducible good locus. Their union is
therefore path connected. The other two strata are connected nonempty
opens of their displayed one-parameter lines. Nonemptiness is supplied
by the separately checked actual representatives, not by the local
coefficient identities alone.

Simultaneous blowups of the intrinsic node multisection and prescribed
cusp/infinity sections yield a smooth proper family of normal-crossing
pairs. The real horizontal field constructed in divisor charts, glued
tangent to each stratum, gives pair diffeomorphisms along compact base
paths. Complex normal orientation transports positive meridians. Every
blowup centre lies on the removed divisor, leaving the original affine
complement unchanged. Thus any independently certified two-positive-
meridian representative transfers through its entire good stratum.

The ensuing retained-sheet graph contradiction is correctly conditional
on that separate input. It needs d>1 and actual generic retained count
a with d<=2a. Then two generators contribute at most2(d-a-1)<=d-2 edges,
where transitivity needs at least d-1. It does not use an unproved
classification of degree-eight or degree-ten monodromy groups.

I replayed the standalone source normally and under optimization. Both65
always-active outputs agree byte for byte with the frozen output. It
imports no repository producer and includes both coordinate charts,
their overlap, the lost-denominator hostile, the even-map hostile and
all three genus rows. The finite source checks the local coefficients;
the analytic argument above pays the global connected-family transfer.

    python3 04-computation/planar_jc48_sep07_higher_odd.py
    python3 -O 04-computation/planar_jc48_sep07_higher_odd.py

Frozen SHA256:

    source edd6c8c77716bb695c4ed902586c19e68486119fddfbab1845c1bb8614515f73
    output c880a825cd246197e7d2a1d49f23d2cf0ad3b8065c0f9fa42c58f2a875f3d958

No mathematical correction was required. Final class status should be
promoted only after the separate higher_braid audit accepts all three
actual representatives and their rational path certificates.
