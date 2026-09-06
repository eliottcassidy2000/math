# Independent audit of the explicit two-anchor Laurent tail

**Status: PASS: independent analytic audit and exact standard-library controls.** The producer may be promoted to PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED within its stated model. No mathematical correction is required. This proves e4>1/100 and an explicit negative-response tail at the same original first-row zeros. It does not prove negativity over the remaining finite phase interval or actual Laurent noncancellation in general.

Audited producer: `C:/w/continuing2_20260906_anchor/continuing2_20260906_effective_anchor.md/.py/.out`. Read against the incoming proved qualitative result `05-knowledge/results/open_frontier_sep06_quadratic_anchor.md`, including both exact anchors, both contiguous interlacers, the Laurent shifts, and the lower carry. Producer pins at review:

* Source: `3d67e31dbdcd6226da0a086cbbfb4a8d142db24984fe5e58508a370e19bc6927`.
* Output: `a6c7117385f0ae039f18e5a08215bf1974db212b6f7ca9948a1376bbc5c969a8`.
* Report before audit promotion: `14eeea20c250bd53337ca1dd1049f6b14a5eae6098e64adc132f24c8b3bca023`.

## 1. The model and its endpoint signs are correctly typed

The five ordered roots are nonnegative, satisfy sum13 and square sum59, and give the monic polynomial B. Both displayed monic quartics C,D must weakly interlace B. If the ordered B roots are a<=b<=c<=d<=e, weak interlacing implies C(d)<=0 and D(e)>=0. This remains true at repeated or zero roots: the quartic root in [d,e] gives the former sign, while every quartic root is at most e for the latter. No unjustified continuity passage or strict-root assumption enters.

The proof does not use arbitrary independently selected real-rooted polynomials in place of the prescribed contiguous C,D. Dropping D preserves the exact boundary (0,0,3,5,5), where e4=0. The second interlacer is therefore a necessary retained predicate for this argument. No actual-source coefficient realization of every model shape is assumed.

## 2. Independent proof review of the fourth-coefficient floor

Assume epsilon=e4<=1/100. Cauchy gives d+e<=sqrt118<11, so a+b+c>2 and c>2/3. If d<=1, the remaining largest root is at least9, contradicting square sum59. Thus d>1 and e>=d>1. The nonnegative term bcde is at most epsilon, hence b<=3epsilon/2; therefore t=a+b<=3epsilon and f=e5<=3epsilon^2/2.

The exact formal elimination of e3 at d yields

    C(d)=d^2(d-5)^2/3 -5epsilon/21 +2f/(3d).

The nonnegative f term and C(d)<=0 imply |d-5|<1/10. All divisions are valid because d>1. In particular c<16/5, so 5-c>9/5. Eliminating the two small roots from the remaining sum and pair sum gives

    (5-c)(5-e)=-3t+2delta+t^2+delta t+delta^2-ab,
    delta=d-5.

The bound ab<=t^2/4 and the preceding bounds yield absolute value at most2433/8000<9/25. Dividing by the strictly positive5-c proves |e-5|<1/5.

At e, the other exact elimination gives

    D(e)=e^2(7e^2-67e+157)/12 -23epsilon/84 +5f/(12e).

Writing e=5+eta shows the quadratic factor is -3+3eta+7eta^2<=-53/25 for |eta|<1/5. Since e>24/5, the negative leading term is at most-2544/625; the positive f term is at most1/76800. Dropping -23epsilon/84 is an upper bound. Thus D(e)<=-7815143/1920000<0, contradicting weak interlacing at the largest root. Every sign, inequality direction, and endpoint is valid, proving the strict floor throughout the closed class.

## 3. Independent carry compiler and tail proof

The referee constructs O,E directly from their14th-binomial coefficients and multiplies them as ordinary polynomials. This independently recovers

    O^2+z^-1 E^2=sum_j binom(28,2j+2) z^j,

at every exponent, including the Laurent endpoint. It then constructs beta,C_raw,D_raw with their actual z^-1 shifts, forms beta^2+2z C_raw D_raw, and applies the Hadamard product in the z exponent only. Coefficient variables e3,e4,e5 remain formal throughout. This differs from assuming the28th-binomial multiplier or importing the producer's SymPy expressions.

The resulting Q has exact support -1 through8, with q[-1]=28. The beta-square coefficient at exponent-2 is1 but is killed by the multiplier; its inclusion in the later coefficient-mass upper bound is harmless overcounting. The top coefficients are exactly the producer's q8 and q7, including the crossing term(6/49)e4^2. All coefficient polynomials are nonnegative.

The directly reconstructed P gives the stated original-root relation

    e5 s=(12/7)e4-e3/s+10/s^2-1/(11s^3).

It applies to every positive original zero, with no assumption of simplicity or that the zero is the smallest positive one. The term involving e3 and the final negative term can be dropped in an upper bound because all elementary coefficients and s are nonnegative/positive as required. The leading cancellation coefficient is exactly-38346750. Hence the nonpositive e4 e5 remainder can also be dropped, preserving the negative e4^2 term.

Coefficientwise GC,GD<=GB, AM-GM on product(1+a_i), and the binomial maximum give the stated total positive coefficient bound. For s>=1, every term of exponent j<=6, including j=-1, contributes at most q_j/s in absolute value after division by s^7. The analogous assertion for s<1 is false, and the referee includes the literal s=1/2 hostile; the producer properly restricts its use.

Combining the strict floor with these bounds gives

    Q(-s)/s^7 < -2607579/7000
                +(17194291313831660508/390625)/s.

The exact ratio of the latter constant to the former is734462481750246368/6215625, strictly less than118163898523. Thus the stated inclusive tail endpoint is sound. The integer is not asserted minimal or optimal. Strictness also survives if the coarse final ratio were attained, because the e4 floor itself is strict.

## 4. Independent controls, scope, and freeze

The standalone referee uses only the standard library and its own sparse rational Laurent-polynomial arithmetic. It imports no producer or repository mathematics. It verifies formal identities, every multiplier coefficient, the complete response support and coefficient signs, the exact rational inequalities and constants, a literal product for the strict root tuple(1,3,9,22,30)/5, and ordered interlacing signs. The repeated boundary tuple is checked directly and D(5)=-25/4 confirms the necessary second-row distinction.

Normal and optimized runs are byte-identical LF and pass124 always-active gates. The universal statement is supplied by the proof above; these exact finite controls validate its algebra and hostile boundaries, not a sampled substitute for its quantifiers.

    python -B 04-computation/continuing2_20260906_effective_anchor_audit.py
    python -B -O 04-computation/continuing2_20260906_effective_anchor_audit.py

* Audit source SHA256: `f8646767fc4e36756abce07c8f12d5d82faab514873a921e74ec4157b4b59a42`.
* Audit output SHA256: `bd0a9c51514ac7aeb18f821ae8dbcda7c5c4540280fa29a10052f976d2a8efa6`.

Promotion is justified for this explicit quantitative refinement of the incoming model theorem. The interval0<s<118163898523 remains an unclosed continuous sign problem, not a finite list of speed rows. No additional root-family, all-h transport, or general noncancellation conclusion is inferred.
