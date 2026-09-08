# Independent audit of the nonordinary three-critical boundary

**Status: INDEPENDENT ANALYTIC / SOURCE / FINITE-EXACT AUDIT PASS.**
This audits the candidate [three-cusp boundary proof](planar_jc48_sep08_three_cusp_boundary.md), its complete exact source, and both fresh replays. No correction to the mathematical statement or source is required. Promotion of the primary is owned by the coordinating agent.

The accepted result classifies a particular degree-(4,6) polynomial-normalization boundary and supplies two good literal curves. It does **not** compute their global braid groups, realize a Keller map, or exclude this boundary as a Keller nonproperness curve. In particular, no ordinary-cusp support estimate is inherited at the exponent-five cusp.

## 1. Inheritance and scope checked

I read the current proved [ordinary three-cusp family](planar_jc48_sep08_three_cusp_family.md), especially its complete normal form, original-coordinate cusp-image resultants, proper incidence construction, and marked transport in §§2, 4 and 6. Only those geometric mechanisms are inherited here. Its ordinary-cusp braid relations and its subsequent Keller consumers are not dependencies of the new boundary classification.

For three distinct critical normalization parameters, the cubic derivative of the quartic first coordinate has exactly those three roots, each simple. After an affine parameter change, target scaling and translation, they are −1, 1, s with s≠±1. Divisibility of the sextic derivative by this cubic leaves a quadratic quotient. The target shear changing its constant term to as is invertible and does not change either degree. Thus

    A=(t²−1)(t−s),  U'=4A,
    V'=6A(t²+at+as)

is complete for the declared class, without a bound on a possible second coordinate of a Keller map. Integrating with zero constants gives exactly the two primary polynomials.

At a root e of A, the ordinary jet is

    U''(e)V'''(e)−V''(e)U'''(e)
      =48 A'(e)²(2e+a).

Since A'(e)≠0, one bad jet occurs exactly on a=2, a=−2 or a=−2s. These lines are pairwise disjoint after s=±1 is removed. Consequently there is no omitted intersection giving two bad jets within the stated distinct-critical universe. Repeated critical parameters are a different problem and are not silently included.

The classification concerns local branches at the three critical parameters before the good-locus test. At a bad parameter outside the good locus another branch could still share a cusp target. The primary explicitly separates that incidence issue from the branch-type calculation.

## 2. Local fifth coefficient and all three covers

At the bad critical parameter e, a=−2e and q(e+h)=q(e)+h². With A1=A'(e) and A2=3e−s, direct integration gives the exact identities

    u=2A1 h²+(4/3)A2 h³+h⁴,
    w=(3/2)A1 h⁴+(6/5)A2 h⁵+h⁶,
    w=V(e+h)−V(e)−(3/2)q(e)u.

These were independently reconstructed from the derivatives, rather than accepted only as source literals. Subtracting 3u²/(8A1) removes all terms below order five and leaves coefficient −4A2/5. The local order-two coordinate can then be straightened by an invertible analytic parameter change; the nonzero odd order five gives characteristic pair (2,5). This local nonlinear target test is not a transformation of the global degree-(4,6) family.

The three fifth coefficients and their sole zeros are

| Boundary | Bad e | Fifth coefficient | Its zero |
|---|---:|---|---:|
| a=2 | −1 | 4(s+3)/5 | s=−3 |
| a=−2 | 1 | 4(s−3)/5 | s=3 |
| a=−2s | s | −8s/5 | s=0 |

At each zero A1 remains nonzero. Both entire translated coordinate polynomials become even in h, so the factorization through z=h² is global. The resulting image parametrization is

    u=z²+2A1 z,  w=z³+(3/2)A1 z².

It is birational because its coordinate degrees are two and three. Its common derivative gcd is z+A1; after translating z by −A1 and removing the tangent term this is an ordinary cubic cusp. The point z=0 is smooth because u'(0)=2A1≠0. The other two original critical h-values map to the same intrinsic cusp. These exceptional parameter maps are degree-two covers of cubic images, not sextic normalizations with extra high cusps.

The converse is also complete. The extension [C(t):C(U,V)] divides both four and six, so its degree is one or two. In degree two, the involution of C(t) fixes the unique pole of U. A nontrivial affine involution in characteristic zero is reflection about a finite point e, making U and V even there. Its fixed point is a common critical parameter with zero ordinary jet. Only the declared bad parameter has that property, and the h³ coefficient of U then forces A2=0. This proves that the three listed covers are the entire nonbirational boundary, rather than a list of known exceptions.

## 3. Infinity, relabeling and genus ledger

With z=1/t, X=U/V and Z=1/V have orders two and six. Exactly

    Z−X³=(V²−U³)/V³.

The degree-twelve numerator term cancels, so the order-seven coefficient is the degree-eleven numerator coefficient, namely 4(3a+2s)/5. No unproved replacement of z by an approximate parameter is involved. On the three boundary lines this gives respectively 8(s+3)/5, 8(s−3)/5 and −16s/5. Each vanishes only at the corresponding double cover. Every birational boundary therefore has precisely the infinity characteristic pair (2,7).

For a birational polynomial parametrization of degrees four and six, the projective curve has degree six and normalization P1. The pullback of a generic line has degree six, and there is no base point at infinity. Its arithmetic genus ten and the four unibranch delta invariants give, when all other singularities are nodes,

    10=2+1+1+3+N,

hence N=3. This ledger is not applied to the degree-two exceptional parametrizations.

Reflection sends (s,a) to (−s,−a), so it identifies the first two boundary lines. On a=−2s the affine parameter change

    t=αz+β,  α=(1−s)/2,  β=(s+1)/2,
    s'=(s+3)/(s−1)

sends the bad critical parameter to −1 and the other two to 1 and s'. Its denominator is a unit on the declared base. The inverse formula for s is the same fractional transformation; s=0 corresponds to the excluded cover s'=−3, and no additional allowed parameter chart is discarded. After scaling U,V by α⁻⁴,α⁻⁶ and removing constants, the quotient in V'/U' has linear coefficient two. Adding (3/2)(2s'−K)U, with the primary's explicit K, corrects its constant term. The source verifies the resulting full U and V identities, not merely their jets. All these are genuine affine parameter and invertible linear/affine target operations.

## 4. The intrinsic good locus is actually open

On a=2, s∉{1,−1,−3}, the divided first-coordinate polynomial

    (U(t)−U(e))/(t−e)²

is monic quadratic, and its value at t=e is U''(e)/2=2A'(e)≠0. Thus simultaneous vanishing of the divided U and V polynomials is a finite incidence over the parameter line with no surviving diagonal point. The three exact resultant factors are

    36(s+1)(8s²−27s−3),
    4(s+3)²(8s−3),
    (s+1)(s+3)²(11s−21),

up to the nonzero prefactors explicitly retained in the source. Removing their zeros excludes exactly extra preimages of the prescribed cusp targets, including coincident cusp targets. It does not remove mere equality of first-coordinate values at distinct target points.

The fixed (2,5),(2,3),(2,3) cusp types and (2,7) infinity type admit the same successive relative point-blowup construction used in the inherited proved family. Here the order-two coefficients and the fifth, third and seventh odd coefficients are units. The change from order three to order five at the bad cusp changes its fixed resolution chain; it does not allow transport across the degeneration from an ordinary cusp. No even coefficient is assumed nonzero or divided out.

After the cusp-image exclusions, the normalization lifts through those centres. Its only critical normalization parameters were the three resolved finite cusps and infinity, so the lifted proper normalization is unramified. Its diagonal in the relative self fibre product is consequently open and closed. The off-diagonal double incidence, and the pairwise-distinct triple incidence obtained by removing the three diagonals, are proper over the base. Nontransverse tangent pairs are closed on the double incidence. Their images and the triple images are therefore closed. Removing them gives exactly the claimed intrinsic good locus: every other branch is immersed, and the only remaining multibranch singularities are transverse doubles.

Both literals below belong to this locus. It is thus a nonempty Zariski open subset of an irreducible parameter line, hence path connected over C. Resolving its nodes as a finite étale multisection gives the inherited proper relative-SNC complement transport within this locus. Positive meridians and their neighborhoods can be transported along a chosen path. This is a mechanism for transporting a future actual supplier; it is not already a braid computation. In particular it supplies no identification with the ordinary three-cusp good locus or its D4 presentation.

## 5. Both exact literal controls and the discriminant argument

The source independently forms each resultant from the literal U and V, checks its substitution identity, and proves that it is monic of degree four in the target second coordinate. It then computes both complete discriminants, the residual squarefreeness and coprimality gates, and every actual target fibre gcd. I read those computations and replayed them in both Python modes.

For (s,a)=(2,2), the cusp target points are

    (−19/3,−65/2), (13/3,63/2), (8/3,8),

with sole preimages −1,1,2 and fibre gcds their squared linear factors. The residual node polynomial is

    (u−3)(2187u²−13962u+22139).

It is squarefree and avoids every critical U-value. At a residual root, all four parameter branches are analytic in U because the parameter discriminant is nonzero. The vertical discriminant has order exactly two. Thus precisely one difference of analytic V-branches has order one: this is one transverse ordinary node. Two pairs, a triple collision, or a tangent pair would have larger discriminant order. At the three critical projections the cusp contributions of orders five, three and three exhaust the complete discriminant; no additional singularity is hidden at those fibres.

For (s,a)=(0,2), the cusp targets are

    (−1,11/10), (−1,−21/10), (0,0),

again with the three squared linear fibre gcds. The residual polynomial is

    (9u+5)(625u²+434u+49).

It is squarefree and avoids u=0,−1. At u=−1, the two distinct target cusps contribute five plus three to the discriminant exponent eight. Their equal U-values are harmless because their V-values differ and neither target has another normalization preimage. The same analytic-branch argument gives exactly three further nodes. Thus this literal genuinely validates the projection-collision clause rather than assuming a generic projection throughout.

These discriminants exhaust all possible finite singularities: away from parameter critical values and vertical collisions the curve is locally a separated analytic graph. The unique infinity branch was independently classified above.

## 6. The exponent-five hostile and its minimality

Use right-to-left permutation composition and

    σ=(123),  τ=(345).

Their product is (12345), so g=(στ)² conjugates σ to τ and sends Fix(σ)={4,5} to Fix(τ)={1,2}. The five-term braid holds, while the three-term braid fails. Both supports have size three, their intersection has size one, and the full retained sets have k=2 and n=0. Hence both the ordinary half-support inequality and the proposed inequality 4n≥5k−D fail. The latter compares zero to five. This is already a full-fixed example, so retaining more fixed labels does not repair it.

The analytic ambient-minimality argument is complete. In at most four labels, two supports of common size three or four intersect in at least half their size by cardinality. A moved support of size two is a transposition; two disjoint such transpositions commute, and an odd five-term braid between commuting permutations forces them equal. The identity has no bad overlap. These exhaust the possible moved counts.

As an independent control I enumerated **all** ordered permutation pairs in S2, S3 and S4 without first imposing the producer's equal-moved-count filter, applying each five-map word successively to each label. There are respectively 2, 6 and 24 braid-five solutions, all diagonal. This recovers the producer's total 32 and verifies that the filter omitted no small braid-five witness. The proof of minimality is not inferred from that finite count.

No global representation or actual sheet retention is claimed for the five-letter pair. It establishes precisely the loss of the local ordinary-cusp inequality under the changed braid exponent. Reusing the ordinary three-cusp closure on this boundary would therefore be invalid.

## 7. Reproduction, exact pins and acceptance

Commands run independently in `/tmp/math-wt-planar-jacobian-sep06`:

    python3 04-computation/planar_jc48_sep08_three_cusp_boundary.py
    python3 -O 04-computation/planar_jc48_sep08_three_cusp_boundary.py

Both fresh runs exited successfully and their bytes equal the frozen output: **122 always-active gates, 939 bytes**. The semantic report digest is

    5de3651c77f8f8be47c485a133457e80816441a96641cd0108dbc7c65049ece1

| Audited artifact | Bytes | SHA256 |
|---|---:|---|
| Primary `.md`, pre-promotion | 15,711 | `b656218a1efc96851cfe83b0be5ed0b2df9901239feb0e36ddda83694e1ca041` |
| Source `.py` | 8,722 | `a244aa123e124afd9ff33dc4b7f177cf629bc5d588fae635c9fa66b2ab31b3ee` |
| Frozen `.out` | 939 | `9f2cbbbbab0df8e4e152915cb3ba49d4a3def9811502f6f4b75d0d0edce1f24c` |

Additional independent calculations reconstructed the boundary polynomials by integrating the derivative normal form, recovered all fifth and seventh coefficients, verified the exceptional cubic images' sole intrinsic critical point, and checked the unfiltered small permutation universe above. They supplement the full literal source replay; they are not a parameter census or an actual braid certificate.

**Final acceptance: PASS.** The three-line exhaustion, exactly three degree-two covers, birational cusp and infinity types, affine relabeling, intrinsic nonempty good locus, both exact literals, and minimal exponent-five hostile are correctly scoped. Global braid data and any Keller consequence for this nonordinary boundary remain OPEN.
