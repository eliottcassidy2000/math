# The nonordinary boundary of the three-critical (4,6) family

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The polynomial identities and named controls are exact. No actual
global braid computation, monodromy presentation, or Keller exclusion
is claimed for this boundary family. The braid-five counterexample in
§6 is an abstract finite permutation statement.

## 1. Complete boundary and inherited scope

Start with the complete degree-(4,6) three-critical normal form from the
[proved ordinary three-cusp family](planar_jc48_sep08_three_cusp_family.md):

    A(t)=(t²-1)(t-s),           s!=+-1,
    U=t4-(4s/3)t3-2t²+4st,
    V=t6+(6/5)(a-s)t5-(3/2)t4+2(-as²-a+s)t3+6as²t,
    U'=4A,       V'=6A(t²+at+as).                       (1)

Allow exactly one ordinary cusp jet to vanish. The three disjoint
boundary lines, over s!=+-1, are

    a=2,        a=-2,        a=-2s.                    (2)

Their nonordinary critical parameters are respectively e=-1,1,s.
The other two critical parameters still have ordinary cusp branches.

**Boundary classification.** On (2), the map is birational except at

    (s,a)=(-3,2), (3,-2), (0,0).                        (3)

Every birational boundary has local finite cusp branches of types
(2,5),(2,3),(2,3), and its unique projective infinity branch has type
(2,7). Thus a curve on this boundary whose cusp targets have no other
preimages and whose other singularities are ordinary nodes has exactly
three such nodes. The exceptions (3) are genuine degree-two covers;
they are not higher finite cusps or infinity-eleven sextics with
birational normalization.

The three lines describe one family up to invertible affine parameter
changes and invertible linear/affine target changes. One may use a=2,
s not in {1,-1,-3}. The intrinsic good locus is a nonempty Zariski open
subset of this irreducible parameter line. Two exact good controls are
s=2 and s=0, the latter retaining a harmless projection collision.

This is also exhaustive for birational polynomial normalizations of
degrees4,6 with three distinct cusp parameters, exactly two ordinary
cusps and one nonordinary cusp. Three distinct common critical
parameters exhaust the cubic U', so each is a simple root of U' and
U has local order two. The same derivative integration and genuine
linear target shear giving (1) then apply. It does not classify
repeated critical parameters or additional branches sharing a cusp
target; those lie outside the declared good locus.

The closest mechanisms are the exact ordinary jets and infinity
characteristic coefficients in the inherited family, and the local
even-term removal in the [higher odd-cusp classification](planar_jc48_sep07_higher_odd.md).
The latter's *one-cusp* theorem is not applied to a three-cusp curve.
The hostile is the even double-cover boundary, and the corrected near
miss is counting its parameter singularities as singularities of an
intrinsic sextic normalization. A second hostile is the length-five
braid pair in §6, which prevents transferring the new ordinary-cusp
closure proof to this boundary.

The five live concepts are derivative critical points, local odd
coefficients, birationality, infinity degree loss, and actual versus
abstract higher-cusp braid data. Root supplied the explicit braid-five
hostile; this lane derived the complete boundary geometry and the
two exact good literals. No external priority claim is made.

## 2. Exact local form: fifth cusp or quadratic cover

Write q(t)=t²+at+as. At the unique nonordinary critical e one has
a=-2e, hence

    q(e+h)=q(e)+h².

Set A1=A'(e), A2=3e-s. The distinct-critical hypothesis gives A1!=0.
The following are exact polynomials, not just asymptotic expansions:

    u=U(e+h)-U(e)=2A1 h²+(4/3)A2 h³+h4,
    w=V(e+h)-V(e)-(3/2)q(e)u
      =(3/2)A1 h4+(6/5)A2 h5+h6.                     (4)

Indeed A(e+h)=A1h+A2h²+h³, so w'=(3/2)h²u'. Removing the fourth
even term gives

    w-(3/(8A1))u²=-(4/5)A2 h5+O(h6).                 (5)

All lower coefficients vanish. If A2!=0, the branch has characteristic
orders two and five and is analytically a (2,5) cusp. The subtraction
in (5) is a local singularity coordinate test; it is not used to alter
the degree-(4,6) coefficient family.

If A2=0, both polynomials in (4) are even in h. Thus the entire map
factors through h², not just through a formal local power series.
The complete table is

| Boundary | Bad parameter e | Fifth coefficient | Vanishing value |
|---|---:|---|---:|
| a=2 | -1 | 4(s+3)/5 | s=-3 |
| a=-2 | 1 | 4(s-3)/5 | s=3 |
| a=-2s | s | -8s/5 | s=0 |

These are exactly (3). They have covering degree two: the field degree
[C(t):C(U,V)] divides both4 and6, and the even factorization makes it
at least two.

Conversely, if the map has degree two, its field involution is an
automorphism of P1 fixing the unique pole of U, so it is h -> -h after
an affine parameter translation. Both coordinate polynomials are even
about its fixed point. That point is a common critical parameter and
has zero ordinary jet. The other two jets on (2) are nonzero, so this
fixed point must be the declared bad e. The h³ coefficient in (4)
then forces A2=0. This proves that (3) is the entire nonbirational
boundary, without a generic birationality assumption.

The degree-two cover has a useful intrinsic interpretation. Put z=h².
Its image normalization is

    u=z²+2A1 z,       w=z³+(3/2)A1 z².

This is a cubic image, with an ordinary cusp at z=-A1 and a smooth
point at z=0. The two remaining critical h-parameters map to the
same intrinsic cusp. Thus retaining these exceptional parameters as
three cusps of a birational sextic would be a type error.

## 3. Infinity and a single irreducible parameter family

In the projective chart z=1/t, X=U/V, Z=1/V, one has orders two and
six. The seventh characteristic coefficient is the actual coefficient

    [z7](Z-X³)=[t11](V²-U³)=4(3a+2s)/5.              (6)

The t12 terms cancel and V³ has leading coefficient one; this explains
the coefficient equality without a formal reparameterization assumption.
On the three boundary lines, (6) is respectively

    8(s+3)/5,        8(s-3)/5,        -16s/5.

It vanishes exactly at the double covers (3). Every birational boundary
therefore has infinity type (2,7); there is no separate infinity-nine
or infinity-eleven birational stratum here.

The rational sextic arithmetic genus ten yields, on the declared
nodal good locus,

    10=2+1+1+3+N,       hence N=3.                     (7)

To see that the three starting lines represent one parameter family,
reflection t -> -t sends (s,a) to (-s,-a) and exchanges a=2 with a=-2.
For a=-2s, put

    t=alpha z+beta,
    alpha=(1-s)/2,       beta=(s+1)/2,
    s'=(s+3)/(s-1).                                    (8)

The old bad parameter s becomes z=-1, the old parameter1 becomes z=1,
and the old parameter-1 becomes z=s'. Rescale the target coordinates
by alpha^-4 and alpha^-6 and remove their constant terms. The quadratic
factor in V'/U' now has linear coefficient2 and constant

    K=(beta²-2s beta-2s²)/alpha².

Adding (3/2)(2s'-K) times the normalized U to the normalized V produces
exactly (1) with parameters (s',2). These are genuine invertible affine
parameter and linear/affine target operations. The exceptional values
s=0,1,-1 correspond to the cover or excluded repeated-critical
parameters; no additional chart is discarded. The source verifies
the two resulting full polynomial identities.

## 4. Intrinsic good locus and exact literals

Use a=2, s not in {1,-1,-3}. Exclude extra preimages of each prescribed
cusp target. The inherited original-coordinate divided resultants
remain valid on this boundary:

    Res_t((U(t)-U(e))/(t-e)², (V(t)-V(e))/(t-e)²).

The first divided polynomial is monic quadratic and its value at e
is U''(e)/2!=0. Thus the incidence is finite over the parameter line,
and it has no unremoved diagonal point. Up to the nonzero prefactors
already displayed in the ordinary family, the three resulting factors
at e=-1,1,s are

    Rminus=36(s+1)(8s²-27s-3),
    Rplus =4(s+3)²(8s-3),
    Rthird=(s+1)(s+3)²(11s-21).                       (9)

Remove their zeros. This excludes actual extra cusp-image branches or
coincident cusp targets, not merely coincident U-values.

Resolve the three finite cusp sections of fixed types5,3,3 and the
infinity branch of fixed type7. The units A1 and the fifth or third
coefficients pay fixed successive point-blowup chains, exactly as the
ordinary family's proper incidence argument, with the bad cusp's
fifth coefficient replacing its third. Infinity uses the nonzero
coefficient (6). No further even coefficient is divided out. The
normalization lifts after (9) has removed extra cusp-image preimages.
The lifted proper normalization is unramified; its diagonal is open
and closed in its self fibre product. Off-diagonal double and
pairwise-distinct triple incidences are consequently proper. Remove
the closed parameter images of tangencies and triple points. This
defines an intrinsic open good locus G, precisely where the prescribed
cusp targets have no extra branches and all remaining singularities
are ordinary nodes. Harmless projection collisions are retained.

The following exact literals show G is nonempty. Therefore G is
connected, being a nonempty Zariski open subset of the irreducible
parameter line. The same proper relative-SNC construction as in the
inherited family can transport marked complements *within G*, after
resolving the three nodes as a finite etale multisection. This does
not identify G with the ordinary-cusp locus or transport its D4
relations across the cusp degeneration.

For s=2,a=2 the literal curve is

    U=t4-(8/3)t3-2t²+8t,
    V=t6-(3/2)t4-16t3+48t.                             (10)

The critical parameters are exactly -1,1,2. Their target points are
(-19/3,-65/2), (13/3,63/2), (8/3,8); the target fibre gcds are
respectively (t+1)², (t-1)², (t-2)². Thus each prescribed cusp has
its sole normalization preimage. Put F(u,v)=Res_t(U-u,V-v), monic
of degree four in v. Exact discriminants are

    disc_t(U-u)=-(256/27)(3u-13)(3u-8)(3u+19),
    disc_v F=-(1048576/847288609443)
       *(u-3)²*(3u-13)³*(3u-8)³*(3u+19)^5
       *(2187u²-13962u+22139)².                        (11)

The residual cubic (u-3)(2187u²-13962u+22139) is squarefree and
avoids every critical U-value. At each of its three roots all four
parameter branches are analytic and U is a local coordinate. The
vertical discriminant has order two, so exactly one pair of analytic
V-values has a simple difference. This is exactly one transverse
ordinary node; a triple point or tangent pair would contribute a
higher order. At the critical projections the known cusp exponents
3,3,5 exhaust the discriminant. Thus there are no hidden additional
singularities, and (10) is an exact good representative.

The second literal, s=0,a=2, is

    U=t4-2t²,
    V=t6+(12/5)t5-(3/2)t4-4t3.

Its three distinct cusp targets are (-1,11/10), (-1,-21/10), (0,0),
again with the three squared linear fibre gcds. Here

    disc_t(U-u)=-256u(u+1)²,
    disc_v F=-(1048576/244140625)u³(u+1)^8
       *(9u+5)²*(625u²+434u+49)².                     (12)

The residual cubic is squarefree and avoids u=0,-1. At u=-1 the
distinct target cusps of types5 and3 account for exponent eight;
this is a harmless co-projection, not an extra branch of either cusp.
The same analytic discriminant argument gives exactly three nodes.

## 5. The next braid supplier remains unpaid

The literal (10) retains the U-polynomial of the certified ordinary
three-cusp representative and has six distinct critical projection
values: three cusps and three nodes. It has no extra smooth fold.
It is a cheap next input for the moving-centre rational Rouché engine,
with a new literal resultant and new paths. One cusp's local braid
has exponent five; the other two have exponent three and the nodes
have exponent two.

The source here proves geometry only. It does not supply a common
based access system, global braid words or actual retained subsets.
Numerical word scouts, if run next, must remain HEURISTIC until their
new rational path and local cluster obligations have been paid. In
particular no D4 relation is inherited across the boundary merely
because the polynomial coefficients specialize.

## 6. A minimal five-letter obstruction to the old local inequality

On five labels take

    sigma=(123),       tau=(345).

They satisfy the five-term alternating braid relation

    sigma tau sigma tau sigma = tau sigma tau sigma tau,

but not the ordinary three-term braid relation. Their moved supports
both have size three and intersect in exactly one label. Thus the
ordinary-cusp conclusion j>=ceil(t/2) fails: here j=1,t=3.

The correct odd-five conjugator is g=(sigma tau)². It satisfies
g sigma g^-1=tau, and maps A=Fix(sigma) onto B=Fix(tau). Nevertheless

    D=5,       k=|A|=|B|=2,       n=|A intersect B|=0.

This also refutes the suggested higher-cusp inequality 4n>=5k-D,
whose two sides are zero and five. It is a full-fixed-set example,
so missing actual fixed labels cannot repair that abstract implication.
It supplies no Keller map or actual global complement representation.

The ambient size five is minimal for failure of the ordinary
half-support bound under a length-five braid. On at most four labels,
moved count3 or4 already forces sufficient intersection by set sizes.
Moved count2 gives transpositions; two disjoint transpositions commute,
and the odd-five braid would then force them equal, a contradiction.
The source also exhausts all such small permutation pairs as a control.

The source-target connection is now precise: the local degeneration
replaces a three-term braid by a five-term braid. It preserves positive
meridian types and the common moved count, but destroys the ordinary
half-support inequality used in the preceding closure. The missing
sidecar is actual higher-cusp local/global monodromy in a common access
system. The explicit five-letter pair is the cheapest decisive test
against reusing the old inequality.

## 7. Reproduction and freeze boundary

Source: [planar_jc48_sep08_three_cusp_boundary.py](../../04-computation/planar_jc48_sep08_three_cusp_boundary.py).
Output: [planar_jc48_sep08_three_cusp_boundary.out](planar_jc48_sep08_three_cusp_boundary.out).
It uses exact SymPy polynomial arithmetic and integer permutations,
with all checks active under optimization.

    python3 04-computation/planar_jc48_sep08_three_cusp_boundary.py
    python3 -O 04-computation/planar_jc48_sep08_three_cusp_boundary.py

The source checks the complete local polynomials on all three boundary
lines, the exact infinity numerator coefficient, all three cusp-image
resultants, the genuine parameter/target relabeling identities, both
full literal resultants/discriminants and cusp fibres, and the explicit
braid-five hostile. Its only permutation bank is the complete D<=4
minimality control, containing 32 pairs satisfying the declared filters.
These controls do not claim a finite parameter census or an actual
global braid certificate.

The runs pass 122 always-active gates with semantic report SHA256
`5de3651c77f8f8be47c485a133457e80816441a96641cd0108dbc7c65049ece1`.
Normal and optimized runs are byte-identical. The source and output
are frozen with these pins:

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source `.py` | 8,722 | `a244aa123e124afd9ff33dc4b7f177cf629bc5d588fae635c9fa66b2ab31b3ee` |
| Output `.out` | 939 | `9f2cbbbbab0df8e4e152915cb3ba49d4a3def9811502f6f4b75d0d0edce1f24c` |

The [independent final proof/source audit](planar_jc48_sep08_three_cusp_boundary_audit.md)
accepts the complete geometry, all122 normal/optimized gates, both literal
controls and the minimal braid-five hostile. Actual braid data remain a
separate supplier; this geometry theorem does not claim a Keller exclusion.
