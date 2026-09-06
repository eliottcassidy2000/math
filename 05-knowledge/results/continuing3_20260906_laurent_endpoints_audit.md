# Independent proof and coefficient audit: endpoints39,45,51,57,63

**Status: PASS.** The five named actual Laurent families may be promoted to PROVED ANALYTICALLY + COMPLETE EXACT CERTIFICATE + INDEPENDENTLY AUDITED, relative to the real-root supplier in THM-4436 and its explicitly cited finite-preserver input. No mathematical repair was requested. The statement is exactly h=6,...,10, not an all-h result or an arbitrary endpoint-at-most63 classification.

Audited producer stem: `continuing3_20260906_laurent_endpoints`. The full proof, source, output and complete3170-coefficient JSON were read. Review pins before status promotion:

* Source: `3b5e9127cd02881afbd2a3cfe8041c990cd2a8688a4eeac86a4cb389c110eae8`.
* Report: `0901f1e9d1d83a2ec08f2e1a8f031082ddf346f33715f2aeaafb7a258ef50d54`.
* Output: `f5ea177b5b3d3d49bd29dd677a35b8244a17ca5c38eb90e44867f5684b56c6e4`.
* Complete JSON: `b55bec0ba7f1063396af1a8db9be725f279fd01a403adb94b9b026b9490d8e62`.

## 1. Source coefficients and first-return quantifiers

For N=6h+3, g>=3h+2, x=g-3h-1>=1, the support(-N,2g-N,3g-N) has exactly one negative and two distinct positive exponents. At mass m, the charge equation is

    -Nm+2g n_beta+3g n_gamma=0.

If gcd(g,N)=1, it implies g divides m. At mass g the nonnegative solutions are exactly(x+j,3h-3j,1+2j),0<=j<=h. At mass2g they are exactly(2x+e,6h-3e,2+2e),-1<=e<=2h. Solving the charge equation for n_beta and enforcing nonnegativity gives these complete sets; there is no omitted congruence or positivity filter. In particular e=-1 is a real channel, since2x-1>=1.

For tau=alpha gamma^2/beta^3 and X=alpha^x beta^(3h) gamma, the first monomial at j=0 is X and its jth variation is tau^j. The analogous second variation is X^2 tau^e, including the negative exponent. Direct multinomial division gives exactly the producer's p_x,q_x and its nonzero scales binom(g,2h+1) and(2g)_(4h+2). The latter has positive length4h+2 and positive lowest factor2x+1; its cancellation with the literal coefficients is legitimate for every x>=1.

The theorem's first-detection conclusion follows only after combining this arithmetic with noncancellation: no mass below g or strictly between g and2g can return at all. This distinction is correctly preserved. Removing gcd(g,N)=1 invalidates the first-return conclusion. At g=3h+3 the literal triple(1,h-1,1) returns at g/3=h+1. The same divisibility equation proves no earlier positive mass. The independent checker verifies the earlier empty fibres and this attained hostile in all five dimensions.

## 2. The real-root supplier is used within its exact scope

Read **THM-4436**, `01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md`, including its declared proof relative to a cited finite real-root preserver. Its complete factorial-row theorem applies with A=2,B=3,r=0,z=1 and the present integer x>=1. Multiplying its polynomial by the positive normalization does not change roots. It supplies exactly h distinct strictly negative roots of p_x. It does not itself prove independence of the two moments; the quotient response below supplies that missing step.

The certificate can be interpreted for real x>=1, but this audit does not extend the real-root supplier to noninteger x. The producer correctly uses its real-parameter positivity only to assert negativity at any real root that exists; the all-roots negative/simple supplier is invoked at integer x. A root rho of p_x is never zero because p_x(0)>0 for x>=1.

The h cancellation phases are h distinct values of the quotient coordinate tau, not h individual coefficient triples. Every phase is attained by beta=gamma=1,alpha=rho. All coefficients remain nonzero and X is nonzero. Choosing all coefficients1 gives the noncanceling mass-g alternative. Thus both first-detection values are attained for every allowed h,g.

## 3. Polynomial response, inverse carry, and degree bound

The source begins in Q(x)[t]/(p_x), where monicity is clear and t is invertible because p_0 is nonzero. Its transition to a polynomial coefficient representative needs the exact cancellation q_(-1)/p_0. Separating the even factors in the2h+1-term falling factorial cancels x+1,...,x+h and leaves

    [2^(h+1)(3h)!/((6h+3)!(2h+1)!)]
       x product_(j=0)^(h-1)(2x+2j+1).

This is in Q[x] with degreeh+1. It is not a claim that the inverse of t lies over Q[x] before multiplication by the carried coefficient. The sign in t^-1=-sum_(j=1)^h p_j t^(j-1)/p_0 is correct.

Give x and t weight1. The first coefficient p_j has degreeh-j. Every nonnegative q_e has degree2h-e, and monic division by p preserves the total weight at most2h. The inverse contribution at t^j has degree(h+1)+(h-j-1)=2h-j. Therefore the exact remainder R_j has degree at most2h-j. Multiplication by t^j gives matrix entry(i,j) degree at most2h+j-i. In a k-by-k principal determinant, the row and column index sums cancel, bounding the characteristic coefficient c_k by2hk. This argument is uniform in h and independent of a finite parameter bank.

The producer's polynomial Faddeev-LeVerrier construction and terminal Cayley-Hamilton identity are consistent with the intended characteristic polynomial. The independent path below does not reuse that algorithm, its monic long-division routine, or its symbolic factorial rows.

## 4. Independent complete identity verification

The referee reads the complete JSON as data and checks the positive rational content, strictly positive integer coefficients, exact degree/list length, primitive integer content, and all five family labels. It then reconstructs actual multinomial fibres by a different algorithm: enumerate n_gamma, solve the original charge equation for n_beta, and enforce all three nonnegativity constraints. The rows are normalized from these literal coefficients, rather than constructed from the producer's falling-factorial expressions.

At each integer parameter, the quotient response is built by Horner evaluation in the companion operator. The t inverse is checked literally by multiplication. All matrix denominators are cleared to produce an integer matrix. Its characteristic polynomial is computed with **Berkowitz**, distinct from the producer's symbolic polynomial Faddeev-LeVerrier method. The integer coefficients are then divided by the corresponding powers of the clearing denominator and compared exactly with Horner evaluation of the saved shifted-coefficient polynomials.

The complete identity bank is:

| h | x values | Highest coefficient degree | Characteristic comparisons |
|---|---|---|---|
|6|1,...,73|72|438|
|7|1,...,99|98|693|
|8|1,...,129|128|1032|
|9|1,...,163|162|1467|
|10|1,...,201|200|2010|

In total this gives665 exact quotient matrices,5640 characteristic coefficient comparisons and18915 complete literal source fibres. Of the665 parameter values,397 meet the final first-return gcd condition; the other values remain legitimate algebraic identity controls. Every equality agrees with the complete JSON.

For each c_k, the number of distinct x values exceeds its proved degree bound2hk. The saved polynomial and independently reconstructed characteristic coefficient must therefore be identical as polynomials. This is a complete bounded-degree identity certificate, not a negative sign scan or a conjecture inferred from finite roots. No numerical roots or floating-point evaluations enter.

## 5. Positivity gives the same-root noncancellation

Every c_k has positive rational content and strictly positive coefficients in x-1. It follows that c_k(x)>0 for every real x>=1. The monic characteristic polynomial is then strictly positive for every real z>=0, including z=0.

At a real p_x root rho, evaluation gives an algebra homomorphism from the quotient algebra to R. Cayley-Hamilton applied to multiplication by R_x yields C_x(q_x(rho))=0. The real number q_x(rho) cannot be nonnegative and is therefore strictly negative. This remains valid without simplicity, although the supplier later gives simplicity for the actual integer parameters. The complete lower carry and the same original first root are retained throughout.

For the actual complex coefficients, the raw doubled moment is a nonzero complex scalar X^2(2g)_(4h+2) times this strictly negative real normalized value. Its raw sign is not asserted, but its nonvanishing follows. Combined with Section1, this proves first detection exactly g or2g and excludes a common first/doubled zero in each of the five named families.

## 6. Frozen referee and remaining boundary

The independent source imports only the standard library and SymPy's integer Berkowitz implementation, with no producer or repository mathematical imports. Normal and optimized outputs are byte-identical actual LF and pass14789 always-active exact gates. The default run is the full identity audit; the optional --quick flag is explicitly diagnostic and is not the frozen output.

    python -B 04-computation/continuing3_20260906_laurent_endpoints_audit.py
    python -B -O 04-computation/continuing3_20260906_laurent_endpoints_audit.py

* Audit source SHA256: `94a1ab60946e16920f0c53104f104087cb376bb3679890ced9e7cd82c3871658`.
* Audit output SHA256: `5f138b8db1ecf383c1984ff1f7e12f549e087aa1d5a79288af8cd9d4896b3c46`.

Promotion is justified for the exact statement h=6,7,8,9,10 and every allowed integer g. The all-h characteristic-positivity conjecture, arbitrary support shapes with endpoint at most63, and the finite-phase sign problem for the separate anchored semialgebraic model remain OPEN. Nothing in this audit supplies a positivity-preserving recurrence between h values.
