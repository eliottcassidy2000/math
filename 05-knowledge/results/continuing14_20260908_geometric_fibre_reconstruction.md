# Geometric generic-fibre reconstruction from degenerating quadruples

**Status: PROVED ANALYTICALLY + FINITE-EXACT; independently accepted
analytically by the root referee.** This strengthens the previously
audited C(T)-isomorphism obstruction. The earlier narrower statement was
correct; its frozen packet remains unchanged.

## 1. Strengthened theorem

Take any two members of the continuing13 mutation supplier, with
Q_i'=f_i^2, each f_i squarefree of degree q_i>=3, and the original
generic affine source curves

    U_i=Spec C(T)[z,1/f_i(z),1/(T-Q_i(z))].

Let L be an algebraic closure of C(T), with the same target parameter
T on both sides. Then

    U_1 x L is L-isomorphic to U_2 x L
    iff Q_2(az+b)=Q_1(z) for some a!=0,b in C.             (1)

Consequently the two degree-fifteen examples in
[the same-torsion pair](continuing14_20260908_same_torsion_different_fibres.md)
are nonisomorphic even as geometric generic affine source curves. Their
equal full pointed source torsion, component counts, pair degree thirty,
and ordinary puncture count twenty-three are unaffected by this stronger
separation. The argument does not classify other global first functions
or identify total fibrations from their generic fibres.

The C(T) proof distinguished residue degrees of closed punctures. That
information disappears after algebraic closure. The replacement is an
intrinsic cross-ratio count at the place T=infinity. No stable-model
existence or uniqueness theorem is needed.

## 2. Exact two-scale puncture configuration

Fix an extension to L of the valuation at T=infinity, trivial on C.
For a polynomial Q of degree N=2q+1, choose t in L with T=t^(-N),
and normalize nu(t)=1. Its punctures split into

    I = the q fixed finite roots of f;
    O = the N moving roots of Q(z)=T, together with infinity.

Every moving root has valuation -1. This follows by comparing the
unique dominant term of Q(z) when nu(z)<0; nonnegative valuation
cannot produce T of negative valuation. The scaled roots w=tz are
units satisfying

    t^N Q(w/t)-1=0.

The reduction of this polynomial is a_N w^N-1, with a_N!=0, and has
N distinct nonzero roots in C. To justify distinct reductions directly,
its discriminant is a unit. All its roots are units and all pairwise
differences have nonnegative valuation. The discriminant product then
forces every pairwise difference to have valuation zero. Thus all N
scaled moving roots have distinct nonzero reductions. The added point
infinity is a further distinct point of P1 in this scale.

In coordinate z the points of I have distinct finite constant reductions,
and all moving points and infinity reduce to infinity. In coordinate tz
the outer points have distinct reductions, none equal to zero, whereas
every inner point reduces to zero. This proves there are exactly the
two scales used below; repeated moving-root reductions or additional
collisions cannot have been suppressed.

## 3. A valuation-invariant count recovers the constant punctures

Call an unordered set of four distinct punctures *degenerating* if a
cross ratio lambda has reduction in {0,1,infinity}. Equivalently,

    nu(lambda)!=0 or nu(1-lambda)!=0.

The condition is independent of the ordering of the four points, since
the six cross ratios are related by lambda, 1-lambda and inversion.
It is preserved by every Mobius transformation over L, even one with
nonconstant coefficients of arbitrary algebraic degree over C(T).

Precisely the quadruples with two points in I and two in O degenerate.
For four inner points, or three inner and one outer, coordinate z has
four distinct reductions. For four outer points, or three outer and one
inner, coordinate tz has four distinct reductions. These quadruples
therefore have cross-ratio reduction outside {0,1,infinity}.

For a,b in I and c,d in O, use

    lambda=(a-b)(c-d)/((a-c)(b-d)).

If c,d are finite, the numerator has valuation -1 and the denominator
-2, so nu(lambda)=1. If one is infinity, the limiting expression gives
the same valuation. Thus all, and only, these two-plus-two quadruples
degenerate.

The number of degenerating quadruples containing a particular puncture
is consequently

    D_I=(q-1) binom(2q+2,2)       for a point of I,
    D_O=(2q+1) binom(q,2)         for a point of O.        (2)

For q>=2 these are distinct, with D_I/D_O=2(q+1)/q>1. Therefore I is
intrinsically recovered as the punctures with larger count. For the
degree-fifteen pair q=7, the counts are exactly 720 and 315.

This counting argument is the elementary substitute for the two-component
stable-limit picture with seven versus sixteen marks. It also explains
why the unequal cluster sizes matter. Equal cluster sizes can have a
nonconstant Mobius symmetry interchanging the clusters; an explicit
control is recorded in the engine.

## 4. Proof of the geometric reconstruction iff

An L-isomorphism of the affine curves extends to their unique smooth
projective completions P1_L and bijects their puncture sets. Their
geometric puncture counts are 3q_i+2, so q_1=q_2. Cross ratios are
preserved, as are the valuation and hence the counts (2). The extended
Mobius map therefore sends the q constant finite points of the first
configuration to the q constant finite points of the second.

Because q>=3, three distinct constant inputs have three distinct constant
outputs. Their unique Mobius map lies in PGL2(C), so the original map
has constant coefficients. Its image of infinity is constant and outer.
The sole constant outer puncture is infinity: no constant can solve
Q(z)=T. Therefore it fixes infinity and is phi(z)=az+b over C.

Choose any moving root z0 of Q_1(z)=T. Its image is a moving root for
Q_2, so Q_2(phi(z0))=Q_1(z0). But z0 is transcendental over C, since
T is transcendental. A polynomial over C vanishing at z0 must vanish
identically. This yields Q_2 o phi=Q_1. Conversely that identity maps
the derivative root sets and the moving root sets, and hence gives the
required geometric isomorphism. This proves (1).

The fixed-target requirement is retained; no field automorphism moving
T is allowed. The argument works inside a single chosen algebraic
closure and one valuation extension. It does not demand that an
isomorphism respect a model or a marking in advance: the finite counts
recover the needed part of the unlabelled puncture set intrinsically.

## 5. Consequence for the frozen explicit pair and the q=2 boundary

For the frozen pair Q_r=S((-1-r)z^3+r), where
S(w)=6w^5-15w^4+10w^3 and 6r^2-15r+10=0, the two polynomials have no
affine equivalence over C. The previously independently audited proof
uses centering, followed by the fact that an affine symmetry of S must
fix both critical points 0 and 1 because their values differ. The
independent exact invariant a_12^5/a_15^4 also has distinct values.
Since q=7, (1) applies and upgrades the nonisomorphism from C(T) to L.

For q=2, the counts still recover the two constant finite punctures,
but two prescribed constant images do not force a Mobius transformation
to be constant. For example z/(t z+1-t) fixes 0 and 1 while depending
on t. This is a precise stopping reason for this proof, not a
counterexample to geometric reconstruction in the q=2 mutation family.
The q=2 case remains unclaimed by this supplement. The earlier
reconstruction over C(T) continues to cover q>=2.

## 6. Finite controls and reproduction

The independent engine checks all four-element subsets of explicit
two-scale configurations for q=2,3,7, with three cross-ratio pairings
for each subset, using literal Laurent coefficients and determinant
valuations. It verifies the recovered counts point by point. These
are controls for the cluster lemma, not numerical substitutes for the
actual algebraic moving roots; the all-degree discriminant argument
above pays those roots. Exact discriminants for N=5,7,15 and the
equal-cluster and two-fixed-point hostiles are also checked.

There are 18,579 always-active exact gates. The complete finite universe
contains 9,255 quadruples and 27,765 cross-ratio pairings.

Run `python continuing14_20260908_geometric_fibre_reconstruction.py`
and the same command with `python -O`. The frozen output and certificate
are byte-identical with LF newlines. The engine writes its certificate
beside itself, or in 05-knowledge/results after repository relocation.
No previous frozen primary or audit packet is altered.
