# Independent hostile audit of MRB contours, deformation, and inverse clocks

Verdict: PASS. No correction to the submitted mathematical claims was found. The qualifications concerning the absence of a proved modular bridge, the restricted perturbation class, signed reconstruction, and numerical versus certified decisions are necessary and correctly stated.

This audit did not import or replay any submitted script. It independently derived the estimates, evaluated a different Euler enclosure with direct integer-binomial weights, and tested the actual nearest-value midpoint inequalities. Only files in `session31-analytic-audit` were created or changed; no git operation was performed.

## 1. Files reviewed

- `session31-mrb/MRB_CONTOUR_AUDIT.md`
- `session31-mrb/certify_mrb.py` and its recorded `.out`
- `session31-root/ENTIRE_MRB_DEFORMATION.md`
- `session31-root/check_entire_deformation.py` and `.json`
- `session31-root/GREEDY_INVERSE_CLOCKS.md`
- `session31-root/check_greedy_inverse.py` and `.json`

The mathematical and implementation claims were audited. Assertions about literature-wide priority or repository-search coverage were not independently re-audited.

## 2. Euler identity and original remainder are correct

Write T(g)=sum_{j>=0}(-1)^j g(M+j), with Delta g(x)=g(x+1)-g(x). Directly shifting partial sums, taking the limit, and using g(x)->0 gives

T(g) = g(M)/2 - T(Delta g)/2.

Iteration gives the stated finite Euler expansion and remainder (-1)^N 2^(-N) T(Delta^N f). There is no missing factor of two and no assumption of complete monotonicity. For N>=1 the differentiated tail is absolutely summable under the bound below, so this operation is legitimate.

For x>=M>=4 on the circle |z-x|=x/2, its argument is at most pi/6, its modulus is between x/2 and 3x/2, and log|z| is positive. Thus, with L(x)=log(3x/2)+pi/6,

|Log z/z| <= 2 L(x)/x.

The function L(x)/x is decreasing on this range, so |exp(Log z/z)-1| <= 2 exp(b) L(x)/x with b=2L(M)/M. Cauchy's derivative estimate gives

|f^(N)(x)| <= 2^(N+1) N! exp(b) L(x)/x^(N+1).

The N-fold unit-cube representation for Delta^N, followed by the decreasing majorant and the elementary bound sum_{j>=0} u(M+j) <= u(M)+integral_M^infinity u(x) dx, gives exactly

2 N! exp(b) { L(M)/M^(N+1) + [L(M)/N+1/N^2]/M^N }.

The integral term 1/N^2 is correct. The submitted code evaluates every nonrational quantity in the interval context and adds the symmetric analytic error interval. The integer sign, prefix endpoint M-1, N input values, and denominator powers are correct. The printed bound 2.653901274338...e-84 is consistent.

## 3. Independent, stronger enclosure

The independent script uses M=384, N=128 and Cauchy circles |z-x|=2x/3. Since asin(2/3)<pi/4, set

L(x)=log(5x/3)+pi/4, and b=3L(M)/M.

The same proof, with the different radius, gives

B = 3 N! (3/4)^N exp(b) { L(M)/M^(N+1) + [L(M)/N+1/N^2]/M^N }.

The script computes B < 1.511975e-132, using 180-decimal interval arithmetic. It does not compute a difference triangle. Instead it expands the finite Euler transform algebraically as

sum_{k=0}^{N-1} (-1)^k f(M+k) W_k / 2^N,

where W_k=sum_{j=k}^{N-1} binom(j,k) 2^(N-j-1) is computed as an exact Python integer. This independently checks the submitted recursive implementation as well as its parameters.

The new C enclosure lies strictly inside the submitted printed C enclosure. It certifies D>0 and q^4-D>0. The bounds are rigorous analytic-plus-interval certificates, subject to the usual correct-execution assumptions for mpmath's directed interval primitives; they are not claims of independently machine-checked formal proof.

## 4. Entire perturbation and theta obstruction

For h_K(z)=exp(-z) product_{j=1}^K(z-j), all finite-prefix and branch-jump assertions hold exactly. It is an entire polynomial times an exponential. It has polynomial growth on vertical lines, exponential decay on the positive real axis, and hence is admissible for the shifted alternating Abel-Plana formula. It does not obey the stronger uniform half-plane O(log|z|/|z|) condition, which the source correctly disclaims.

Using x=-exp(-1), the alternating perturbation sum is

K! sum_{n>=K+1} binom(n-1,K) x^n = K! [x/(1-x)]^(K+1)
= (-1)^(K+1) K!/(e+1)^(K+1).

The shift is nonzero. Consequently arbitrary real shifts preserve every stated datum, but the theorem cannot exclude a special identity for the unperturbed full analytic profile. Direct numerical finite sums for K=1,...,8 independently agree with the exact expression to better than 1e-160.

With conventional theta definitions,

theta_3(0,q)+theta_4(0,q)-2 = 4 sum_{r>=1} q^(4r^2).

Therefore the unweighted positive even-theta tail is strictly greater than q^4, while the certified defect is strictly less than q^4. The sign obstruction is correct, and it does not apply to arbitrary signed or normalized combinations.

## 5. Contour identities and algebraic tail

The shifted alternating Abel-Plana formula yields C=integral_0^infinity Im f(1+it)/sinh(pi t) dt, since f(1)=0. The imaginary-part convention has the correct sign. A direct formal expansion gives

f(1+u)=u-u^2+u^3/2+u^4/6-3u^5/4+O(u^6),

so g(t)=t-t^3/2+O(t^5). Expanding the positive 1/sinh kernel and using absolute integrability gives the odd-Laplace decomposition.

For odd ell, rotating the integral of f(z) exp(i pi ell z) along the rectangle based at z=1 yields

integral_1^infinity f(x) exp(i pi ell x) dx
= -i integral_0^infinity f(1+it) exp(-pi ell t) dt.

Taking real parts gives the submitted F_ell formula and its sign. Integration by parts has no endpoint term because sin(pi ell)=0 and f(x)->0. Laplace scaling gives 2/(pi^2 ell^2)-6/(pi^4 ell^4)+O(ell^-6), and summing the odd reciprocal-square tail gives 1/(2pi^2 J)+O(J^-2). Thus the natural odd-mode tail is genuinely algebraic.

For the vertical-tail bound, writing r=|1+it|>=1 gives |Log(1+it)/(1+it)| <= [log r+pi/2]/r <= pi/2, which justifies the source's global |f| <= exp(pi/2)-1 estimate. The small-circle warning is valid because the negative-axis boundary values themselves grow as x^(-1/x) near zero.

The AGM identities, modulus-versus-parameter convention, A(8) normalization, and honeycomb algebraic identity were checked directly and are consistent.

## 6. Entire deformation: proof correct, numerical root now certified

For h_k=log(k)/k, H=log(3)/3 is the integer maximum. The linear term is conditionally convergent, but the terms of degree at least two converge absolutely and locally uniformly since sum h_k^2 is finite. This yields the entire function and the coefficient sign

sum_k (-1)^k h_k^j = (-1)^(j+1) eta^(j)(j).

The remainder estimate follows from the exponential Taylor remainder and h_k^(J+1)<=H^(J-1)h_k^2. It is valid also for J=1 and complex t. The source does not incorrectly apply the imaginary-part integral to complex t: its integral discussion explicitly restricts to real t. The derivative endpoint values are correct.

The independent interval method extends to f_t(z)=exp(t Log z/z)-1 with analytic tail B_t obtained by replacing exp(b) by |t| exp(|t|b). For the t derivative h(z)exp(t h(z)), the corresponding factor is exp(|t|b), without |t|. Evaluating interval t in [0.999999,1.000001] therefore gives a uniform certificate for D_8'(t):

-0.089603809 < D_8'(t) < -0.089598986.

The endpoint defect signs are strictly positive and strictly negative, respectively. Continuity and strict monotonicity prove there is exactly one real zero in this interval. Two tighter exact-decimal endpoint evaluations certify that zero lies strictly between

1.000000211135851595066050802731745766168538790567401655840
and
1.000000211135851595066050802731745766168538790567401655842.

This 2e-57-wide bracket strengthens the submitted numerical-only qualification; that qualification was conservative, not erroneous. Nothing here proves global uniqueness of a zero or a modular identity.

## 7. Greedy conventions, inverse asymptotics, and rounding bounds

The logarithmic derivative is correct. For x>=3, the elementary sine/cosine bounds imply pi cot(pi/x)>=x-pi^2/(2x)>log x-1; the last difference is positive at 3 and increasing. Conversely pi cot(pi/x)<x and log x>1. Combining this with A(x)<=2/x proves 0<-A'(x)<=2/x^2.

Expanding exp(-log(x)/x) and sin(pi/x) independently gives

A(x)=2/x-2log(x)/x^2+[(log x)^2-pi^2/3]/x^3+O((log x)^3/x^4).

Reciprocation gives the stated X(x) expansion. Substituting x=X-L+c/X, log x=L-L/X+O(L^2/X^2), and 1/x=1/X+O(L/X^2) cancels the 1/X term exactly when c=L-L^2/2-pi^2/6. Since the derivative of X(x) tends to 1, the remaining O(L^3/X^2) residual gives an inverse error of that same order. This is an asymptotic statement, not an integer enclosure.

For x in [m,m+1], each nearest-amplitude error is bounded by half the adjacent amplitude gap, and integrating the derivative gives the bound 1/m^2. The uniform fractional-part estimate follows from -A'(u)=2/x^2+O(log x/x^3) for |u-x|<=1. Taking the minimum is Lipschitz, so the error remains additive uniformly, including arbitrarily near integers. The warning against dividing by min(theta,1-theta) without a lower bound is essential and correct.

The lower bound A(x)>=c/x uses exp(-log(x)/x)>=exp(-1/e) and the decreasing sinc function on [0,pi/3]. It yields x(r')>=c/r' and hence the all-step next-index inequality exactly as written. The ±1 losses correctly distinguish floor x from the selected neighboring integer.

A concrete adversarial example is x=100.499: rounding x gives 100, but the nearest amplitude is A(101). Thus the source's nonlinear-rounding warning is substantive.

## 8. Four nearest indices independently certified

For each selected k, the new interval calculation proves

r - [A(k)+A(k+1)]/2 > 0,
and
[A(k-1)+A(k)]/2 - r > 0.

Strict monotonicity then gives a globally unique nearest integer. It also separately brackets each continuous inverse, so the over/under labels are not inferred merely from rounded roots.

The certified index sequence is

8,
105719004,
77401040908104410,
24845630402060641702669818431121505.

The first continuous inverse lies in (8,9), and each later inverse lies between its displayed index minus one and that index. The signed residual is negative at the first step and positive at the next three steps. Thus the signed reconstruction through these terms is

C = A(8)-A(105719004)-A(77401040908104410)
    -A(24845630402060641702669818431121505)-r_4,

with r_4>0. The smallest of the eight certified midpoint margins is greater than 2.278349784e-70, well above the independent C uncertainty. The fourth remaining magnitude is approximately 1.3921088628035451e-69. This certificate upgrades the later two numerical decisions; it does not certify a fifth integer.

## Reproduction

Run from the workspace root:

    python session31-analytic-audit/independent_checks.py

Files `independent_checks.json` and `.out` contain the detailed interval evidence and summary. The script has no network, source-file mutation, or git operation.
