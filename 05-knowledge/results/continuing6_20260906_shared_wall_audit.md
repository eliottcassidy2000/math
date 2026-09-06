# Independent referee: an admissible original-phase denominator wall

**Status: PASS.** I accept the exact C/D shared-root charts, the eight positive singular pairs, and the selected admissible C-wall witness in `continuing6_20260906_shared_wall.md`. No mathematical repair is requested. This is an anchored coefficient-model result, not an integer-support Laurent realization. The witness has a negative original response; it obstructs discarding a singular chart fibre, not the desired sign theorem.

The independent source imports no producer implementation. It reconstructs the shared-root eliminations and performs the witness arithmetic in the exact quotient ring QQ[r]/P_C. Different rational root brackets verify the full root geometry. Normal and optimized runs pass **136 always-active exact gates**, with identical raw LF output.

## 1. Both charts, including their excluded denominators

I independently solve B(r)=C(r)=0 and B(r)=D(r)=0 for y,z. The solutions agree with the producer in both directions; the equations are linear with nonzero coefficients for the retained r>0 domain. Substituting s=u/r into the unchanged original phase gives exactly

    g_C = x u^2(5u^2-24u+9)/(9r^2)-H_C/(33r),
    g_D = x u^2(23u^2-60u+12)/(12r^2)-H_D/(11r).

The verifier compares the full H_C,H_D polynomials, not merely the x coefficients. Away from the respective quadratic wall, the displayed rational solution for x is legitimate. On the wall, a nonzero constant excludes a phase; a zero constant leaves x free in the algebraic chart. Dividing there would delete that fibre.

Reducing H modulo its coefficient quadratic gives a linear equation in u. Its coefficient is coprime to the corresponding primitive quartic P_C or P_D in r. The rational u reconstruction is therefore unique at every quartic root. I independently reconstruct both quartic resultants, prove squarefreeness, count four positive roots for each, and check every published isolating interval. Substitution of the reconstructed u into its quadratic vanishes modulo the quartic; the linear remainder vanishes identically. Hence the original H also vanishes: there are no extraneous resultant roots or omitted reconstruction-denominator cases.

Both coefficient quadratics have two positive roots, so all eight recovered pairs have r,u,s>0. Only one of the eight pairs is claimed admissible. The largest r in each chart exceeds71/10 and is excluded even by the anchors: Cauchy's inequality for the other four beta roots gives 5r^2-26r-67<=0, whereas this polynomial is already positive and increasing for r>=71/10. No conclusion about admissibility of the other five unchecked pairs is inferred.

## 2. The exact selected point

The selected r is the unique P_C root in

    (3528453826758,3528453826759)/10^12.

The source verifies the exact opposite endpoint signs and Sturm count, then refines this same interval by **40 rational sign bisections**. It never replaces r by a decimal approximation. Let u be the stated rational function of r, s=u/r, and

    x=(3/5)r(4r^2-45r+150)+1/10,
    y=(14/9)xr-(7/3)r^2(r^2-12r+45),
    z=r^2/18.

All these quantities are evaluated in the quotient ring. The shared equations and original phase vanish exactly. The u interval lies in(0,1), selecting (12-3sqrt(11))/5 from the two quadratic roots. The exact interval tests show x,y,z,s>0. Thus the zero-beta boundary is not being reused accidentally.

The defining polynomial quotient need not be treated as a numerical approximation or a universal ordering of all conjugates. Every inverse used by the arithmetic exists modulo P_C; the one selected real embedding is fixed by its isolating interval.

## 3. Independent root and interlacer proof

Synthetic division by the known v-r factor is performed in the exact quotient ring for B and C, with zero remainders checked. I evaluate their quotient polynomials, and D itself, at these independent rational brackets (all endpoints have denominator1000):

| Polynomial | Brackets | Endpoint signs |
|---|---|---|
| B/(v-r) | (18,20), (669,671), (2450,2453), (6330,6334) | +/-, -/+, +/-, -/+ |
| C/(v-r) | (387,390), (1950,1954), (6129,6133) | -/+, +/-, -/+ |
| D | (181,184), (1492,1495), (3382,3386), (5938,5942) | +/-, -/+, +/-, -/+ |

Every sign is certified by a rational interval over the entire refined r interval. Each set has as many disjoint sign changes as the degree, proving that all roots are real, simple and contained in those brackets. Inserting the exact r root proves:

- B has five distinct positive roots;
- C has exactly one shared B root, its third root equal to the fourth beta root r, and its other roots strictly interlace;
- D strictly interlaces B and shares no beta root.

The coefficient anchors give beta sum13 and square sum59. These proofs retain every root, including the known common factor. There is no cancellation that would silently weaken the original geometry.

## 4. Original phase and full response

The verifier independently rebuilds the length14 odd/even binomial coefficients and both Laurent products directly in QQ[r]/P_C. Their coefficientwise product retains all response exponents, including q_(-1)=28 and the mixed top terms. Evaluating the first Hadamard polynomial at -s gives zero exactly, independently of the previously substituted chart expression.

Evaluating the complete Q(-s) gives a reduced polynomial of degree at most three in r. Rational Horner enclosure on the 40-times refined r bracket proves

    -5634725 < Q(-s) < -5634723.

This is a separate numerical-proof path from the producer's direct multivariate native interval evaluation. Reduction can enlarge coefficients and worsen interval conditioning, which is why exact refinement is done before applying the final field enclosure. No floating-point evaluation is used in the sign certificate.

The response is negative at the selected original phase. The statement does not require, and the audit does not infer, a sign at every phase of this witness or every point of the singular fibre.

## 5. Consequence and precise remaining question

The point is in the full weak C/D model with simple positive beta roots and z>0, yet the original phase's x coefficient and constant both vanish. Moreover D is strict, so this point is not recoverable by silently moving to a D-shared chart. Generic division on both shared charts would lose this valid C-boundary point.

There is a small qualitative consequence of the certified strict brackets. Keep r,u,s fixed and vary x, with y,z given by the C chart. The phase remains identically zero along this affine fibre. The strict endpoint signs for B/(v-r), C/(v-r) and D persist for x in some open interval around the certified value, as do z>0 and the negative response. Thus the admissible wall is locally a genuine one-parameter family, not an isolated point caused by numerical coincidence. No explicit interval endpoints or global sign on that fibre are claimed.

The remaining object is the x-dependent original response on admissible portions of these finite singular(r,u) fibres, together with the generic shared-root charts. The present result validates retaining the singular fibres. General shared-positive-root sign remains open, and no integer-support realization or LRC consequence is attached to this coefficient-model witness.

## 6. Reproduction and frozen pins

Run `continuing6_20260906_shared_wall_audit.py` ordinarily and with `python -O`. The `.out` and `_optimized.out` files are byte-identical raw LF captures. Adjacent dependencies and the standard sibling results layout are supported.

- Audit source SHA256: `1c19836e3381a709604b2d0d27b273b5de1efbc1c1cc1901d7119e2fa6cecba4`.
- Audit output SHA256: `a9dcb1e6357ca7a3d2072f1d24edc7a44f4c788d360a8d48b8c17b0c383de2cd`.
- Producer source: `ef604a3c9276171b844f6e260bc788a17160ad6e85a70f3fbd131ee4f3c362c4`.
- Producer output: `4803b7577514138cbe0fca4ef339b0e2cb690ba421e2276eb823846ef0d33572`.
- Producer certificate: `6adab7b43f30296cdac428596ce7c046d163751b6d48fea4db02fe17a100dce4`.
- Audited prepromotion report: `7050e407746e3cb5b23e4bdaa72f90cc0af487070e255c8f49b27a83c091820b`.

All referee files are outside the repository. No producer or maintained file, namespace, or Git state was changed. Parent owns promotion and integration.
