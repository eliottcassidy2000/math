# Independent audit of the literal two-cusp sextic

**Status: independent analytic and exact-source audit PASS.**
This accepts the geometry, the bound of three positive meridian generators,
and the conditional actual Euler ledger in
[the primary note](planar_jc48_sep08_two_cusps.md). It does not identify the
full complement group or exclude a whole Keller nonproperness curve.
Those questions remain OPEN. No other member of the parameter family is
classified here.

## 1. Scope and source recovery

I read the complete primary proof and its standalone source, then ran the
source under normal Python and `-O`. Both outputs match the frozen output
byte for byte: 42 always-active gates, 505 output bytes. I also recovered
the exact node-page intersection and smooth-stratum constancy from
[the infinity supplier](planar_jc48_sep06_infinity.md), and the monic
quartic's outside-root section and fibre-generation argument from
[the audited next-braid theorem](planar_jc48_sep07_next_braid.md).
The old one-cusp budget is not an input to the new Euler identity.

The theorem is about the literal curve

    U=t^4-2t^2,
    V=t^6-3t^4/2+t^3/3-t.

The named THM-3844 comparison concerns a different quartic and transfers
no class-group or completion theorem. The candidate's three parameter
hostiles are controls, not a classification of their surrounding family.

## 2. Complete geometry, with independent elimination checks

The critical-parameter gcd is exactly t^2-1. Both critical branches have
U''=8 and nonzero ordinary-cusp determinant, and their image ideals have
only the stated double parameter root. Thus they are separate ordinary
cusps with no additional branch at either image, although both project to
u=-1.

The divided difference U(s)-U(t) has the two essential symmetric branches
p=0 and q=(p^2-2)/2, with p=s+t,q=st. The first gives the one unordered
pair with q=-3. The second gives the diagonal cusp pairs and three
non-diagonal pairs from 3p^3-2=0. Its discriminant 4-p^2 is nonzero at
all three roots. The intersection of the two first-difference branches
is not a second-difference zero. This exhausts the collision equations;
in particular no saturation discards the zero-sum node.

I independently replaced the tangency Gröbner calculation by the literal
factorization

    U'(s)V'(t)-V'(s)U'(t)
      =4(s-t)(s^2-1)(t^2-1)(1-6st(s+t)).               (A)

On the zero-sum pair, the last factor is 1 and
(s^2-1)(t^2-1)=4. On each of the other three pairs,

    (s^2-1)(t^2-1)=p^2(p^2-4)/4,
    1-6st(s+t)=6p-1.

Neither expression vanishes at a root of 3p^3-2. All four collisions are
therefore transverse. This independently confirms the stated Gröbner
support and excludes a tacnode or a critical colliding branch.

There is also a short exact alternative to the triple-image unit ideal.
Let r=3 rem_t(V-B,U-A), a monic cubic. Write

    rem_t(U-A,r)=c2 t^2+c1 t+c0,
    c2=9(A+1)^2+1,
    c1=3B-21A/2-9,
    c0=9A(A+1)/2-9B(A+1)-A.

Direct division gives these expressions, and

    c0+3(A+1)c1+3c2=3-A,
    c2+9(A+5)[c0+3(A+1)c1+3c2]=145.                 (B)

The coefficients cannot all vanish. Three distinct common roots of
U-A and V-B would force the monic cubic r to divide U-A, contradicting
(B). This proves that the four pairs have distinct node images. The
constant leading coefficient of r is essential: no parameter-dependent
division denominator has been discarded. I checked (A), the division,
and both identities in (B) separately from the producer's Gröbner call.

Finiteness follows from the monic equation for t over the image ring;
the finite collision list implies generic injectivity and birationality.
The degree-six projective parametrization has no basepoint: on the
finite chart its last coordinate is nonzero, while at its unique
infinity parameter its V-coordinate is nonzero. Therefore its image is
an irreducible sextic with normalization A1 on the affine chart.

For an independent infinity check set

    d(z)=1-3z^2/2+z^3/3-z^5,
    X=z^2(1-2z^2)/d(z),  Z=z^6/d(z).

Then

    Z-X^3=z^6[d(z)^2-(1-2z^2)^3]/d(z)^3,
    d(z)^2-(1-2z^2)^3=3z^2+2z^3/3+O(z^4).

It follows that Z-X^3-3X^4 has leading term 2z^9/3. Since X has order
two, the branch has its first characteristic odd order nine, hence
Puiseux pair (2,9). The line at infinity has contact six. The local
contributions 1+1+4+4=10 exhaust the arithmetic genus of a sextic.

I read and replayed the full actual resultant, its monicity and minimality,
its complete discriminant, and the three special-fibre factorizations.
The node projection polynomial is squarefree and avoids -1,0,3. At u=0
there is exactly one colliding pair, at t=0; U''(0)=-4 and V'(0)=-1 make
it a simple vertical fold at a smooth point. Its two other vertical
points have coordinates 2 plus or minus sqrt(2)/3, so they are distinct
and nonzero. At u=-1 the two cusp targets remain distinct despite their
shared projection. These checks pay the actual geometry needed below.

## 3. A single fold gives at most three positive meridians

For the monic quartic, choose a continuous real vertical section above a
Cauchy bound for all roots. It exists over the whole u-plane, including
critical values. Over the regular base, the complement is the four-point
fibre bundle. Its group is generated by fibre loops and section lifts of
base loops. Restoring the exceptional fibres kills those section lifts,
and arbitrary loops can be perturbed to miss the exceptional vertical
fibres. Hence the four positive fibre meridians surject onto the full
affine complement group.

An access path from the regular base to a small circle about u=0 transports
its simple smooth fold to a half twist on one embedded vanishing arc in
the four-punctured fibre. A braid change of positive geometric basis takes
that arc to a standard adjacent arc. Such a basis change preserves both
free generation and the property of being a positive curve meridian.
In this adapted basis the local relation is the fixedness of a single
half twist. For the convention

    H+(a,b)=(aba^(-1),a),

fixedness is precisely a=b; the inverse convention gives the same equality.
The relation is exact in the group because the outside-root section loop
is null after the critical fibre is restored; there is no residual inner
conjugation being mistaken for fixedness. Thus two of the four adapted
positive meridians are equal, and the other three generate the group.

This proves the positive-meridian generating bound without guessing an
access-path word. It does not determine the simultaneous conjugators of
the two cusp relations and four node relations. The two cusp cubes can
be put on disjoint local arcs, but this supplies no common global gauge
relative to the fold-adapted basis. There is no claim here of two
meridian generators or of a particular transitive permutation action.

## 4. The revised actual Euler ledger

This paragraph is conditional on the curve being the **whole irreducible
nonproperness curve of a nonautomorphic polynomial Keller map** of generic
degree d. Let a be its actual retained count on the whole smooth stratum,
whose constancy is supplied by the inherited local argument. Let n1,n2
be the actual cardinalities above the two cusps. For each node let omega
be the overlap of the two deleted subsets, each of size d-a. The exact
retained-sheet intersection lemma gives node fibre cardinality

    2a-d+omega.

Normalization is A1. Two unibranch cusps change neither chi_c(C)=1-N
nor chi_c(C2\C)=N, but removing them from the smooth stratum changes its
Euler characteristic to 1-2N-2=-1-2N. Constructible Euler integration
therefore gives

    1=dN+a(-1-2N)+sum_nodes(2a-d+omega)+n1+n2
     =-a+n1+n2+sum_nodes omega.                       (C)

This proves the displayed revised ledger with actual fibre counts;
no arbitrary fixed-point counts have been substituted. For this curve
N=4. In particular the one-cusp equality n+sum omega=1 is not inherited,
and its ensuing half-degree and small-cusp-passport exclusions cannot be
silently used here. The three-positive-meridian bound and (C) leave the
full group and whole-Keller exclusion OPEN.

## 5. Frozen replay and acceptance

Reproduction commands are exactly

    python3 -B 04-computation/planar_jc48_sep08_two_cusps.py
    python3 -B -O 04-computation/planar_jc48_sep08_two_cusps.py

All 42 checks are explicit and active under optimization. The fresh normal,
fresh optimized, and frozen output are byte-identical. Hashes independently
recomputed during this audit:

- Source, 4,663 bytes: `d3c1bffe31c3147557337066dcce1ee4dc97d19cd2f77832244fb1b0317f4af0`.
- Output, 505 bytes: `c1981a35fe9293387f60c356c3d4bad70ef58e25ef65ad1c35c7779c122831a5`.

No source/output changes were made. The sole requested prose clarification
is to put the whole-irreducible-Keller hypothesis and actual smooth-stratum
meaning of a explicitly before the primary's Euler ledger, as in Section 4
above. No mathematical formula requires correction. The literal geometry,
three-positive-meridian bound, and this conditional Euler identity pass.
