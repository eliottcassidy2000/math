# Decimal prime densities: a negative square, not a positive spectral measure

Status: PROVED from exact prime counts; FINITE-EXACT arithmetic independently replayed. No claim of priority over the moment-problem literature. No implication that RH is false or that all spectral approaches fail.

## 1. The sequence and its native coordinate

Let delta_k be the proportion of primes among positive integers with exactly k decimal digits:

    delta_k = [pi(10^k)-pi(10^(k-1))]/[9*10^(k-1)], k>=1,

where pi(1)=0. Powers of ten are composite, so open versus closed upper endpoints do not change these counts. The first terms are 4/9, 21/90, 143/900. They count primes, not repunits, and the data alone did not specify the sequence until this interpretation was identified.

Define a linear functional L on real polynomials by L(x^n)=delta_(n+1), n>=0. The following exact identity is a decisive hostile test:

    q(x)=x(x-1)^2(x+1)(2x-1),
    L(q(x)^2) = -2447633/22500000000 < 0.

In expanded form q has coefficients (0,-1,3,-1,-3,2), and q^2 has coefficients

    (0,0,1,-6,11,0,-21,18,5,-12,4).

Thus this certificate uses only delta_3, delta_4, delta_5, delta_7, delta_8, delta_9, delta_10, delta_11. The script `check_prime_density_square.py` evaluates it using exact fractions. It also checks the positive control integral_0^1 q(x)^2 dx=1/396.

## 2. Exact spectral obstruction

There is no positive Borel measure mu on the real line with finite moments of all orders such that

    delta_(n+1)=integral t^n dmu(t), n>=0.

Otherwise L(q^2)=integral q(t)^2 dmu(t)>=0, contradiction. Equivalently there is no self-adjoint operator A and vector v in the required power domains with

    delta_(n+1)=<v,A^n v>, n>=0.

The spectral theorem would make L(q^2)=||q(A)v||^2 nonnegative. This excludes an arbitrary real spectrum, not merely a positive operator. It does not exclude signed or off-diagonal matrix elements, a different observable, a regularization, nonlinear transforms, or models which retain additional data. It also does not assert the same obstruction for every sufficiently late shifted tail.

The natural shifted generating series M(z)=sum_(n>=0)delta_(n+1)z^n therefore cannot equal <v,(I-zA)^(-1)v> as a genuine moment-generating germ with these moments. If a proposed resolvent has a different spectral coordinate or only a formal asymptotic expansion, those hypotheses must be checked separately.

## 3. Why its asymptotic singularity alone would mislead

For any fixed integer base b>=2, the prime number theorem gives

    delta_(b,k) ~ 1/(k log b).

Subtract the two prime-counting asymptotics at b^k and b^(k-1); their ratio stays fixed, so their leading difference is nonzero. For the ordinary generating series D_b(r)=sum_(k>=1)delta_(b,k)r^k, elementary Abelian comparison with sum r^k/k gives

    D_b(r) ~ log(1/(1-r))/log b, r -> 1 from below.

This is nonrational: a rational function near r=1 is bounded or has an integer-order pole, not logarithmic divergence. Nor is D_b(q), with q=exp(2*pi*i*tau), a nonzero holomorphic level-one modular form of positive integral weight: its zero constant term would make it a cusp form, and the S-transformation would force exponential decay times a power along tau=iy ->0, contradicting the divergence.

But a logarithmic singularity is entirely compatible with continuous self-adjoint spectrum. Multiplication by t on L^2([0,1],dt), observed through the constant vector, gives

    <1,(I-zA)^(-1)1> = integral_0^1 dt/(1-zt) = -log(1-z)/z.

The obstruction to the actual prime densities is the explicit negative square, not the logarithm. A branch singularity also differs from an isolated eigenvalue pole.

## 4. A stronger positive comparator: the logarithmic-integral block model

For k>=2 define

    delta^li_(b,k) = (1/(b-1))*integral_1^b du/[(k-1)log b+log u].

This is exactly [li(b^k)-li(b^(k-1))]/[(b-1)b^(k-1)], with the integral taken away from the singularity at1. Its shifted tail is a positive Hausdorff moment sequence:

    delta^li_(b,n+2)=integral_0^1 x^n w_b(x) dx,
    w_b(x)=(b*x-1)/[(b-1)log(b*x)] >0.

The value at x=1/b is defined by continuity as1/(b-1). To prove the identity, insert

    1/[(n+1)log b+log u]
       = (1/log b)*integral_0^1 x^(n+log(u)/log(b)) dx

and integrate u. Thus the smooth prime-density model really has positive continuous-spectrum realization; it is not just visually similar to one. The actual arithmetic densities need their own positivity test, and fail it. The parallel recovery report supplies an additional negative-square certificate for the actual shift delta_(n+2), matching the smooth model's index.

## 5. Sources, provenance and limits

Counts through10^7 were independently sieved in the recovery lane. Counts10^8 through10^11 are recorded in Andrew Granville, *A good new millenium for the primes*, Table1, printedp2: https://dms.umontreal.ca/~andrew/PDF/PrimesICM1.pdf . The proof above treats those exact counts as cited arithmetic inputs, not a claimed local exhaustive sieve through10^11. The separate source-check report records the table inspection and reproduction details.

Prior repository mechanisms: `01-canon/theorems/THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go.md` uses Hankel positivity/rank to test faithful carriers; `01-canon/theorems/THM-3056-product-gamma-reciprocal-hypergeometric-and-hankel-reversal.md` already uses a negative Hankel determinant to distinguish positive coefficients from positive moment realizations. This session applies that established obstruction to the user's decimal prime-density sequence.

The finite negative-square identity is unconditional. The prime-number-theorem asymptotic is a separate classical input. No zero-location conjecture, finite Weil test, or unproved coercivity estimate is assumed.
