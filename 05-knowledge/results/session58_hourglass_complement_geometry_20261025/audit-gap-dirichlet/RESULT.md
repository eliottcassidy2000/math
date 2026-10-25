# Independent audit: candidate Dirichlet series

Verdict: PASS for the producer's communicated factorization and asymptotic statements, conditional only on the already-audited Gaussian indexing and balanced-box counting theorem. No RH, generalized RH, zero-free half-plane assumption, or claim about existence of actual additive gap relations is used.

## 1. Exact coefficients and ancestry

For a center root e, total N(e) is the balanced-box count from its split-prime exponents. Define P(e)=0 if e=1 or if e has a prime factor 2 or congruent to 3 modulo 4; otherwise P(e) is the primitive full-maximum box count. A candidate's unique maxima vector gives its unique primitive center d dividing e. Missing coordinates are identically zero; insertion restores the original candidate. Thus N(e)=sum_{d|e}P(d), as an actual disjoint ancestry partition, not merely a formal identity. Center roots, not squared centers, are the summation variable.

Writing X=(L,V,U,W), use weights (1,-6,-3,14), polynomials

- Q_L(t)=1+18t+5t², r_L=3, c_L=10;
- Q_V(t)=1+4t−t², r_V=3, c_V=3;
- Q_U(t)=1+3t, r_U=2, c_U=2;
- Q_W(t)=1+t, r_W=2, c_W=1.

The total component A_X has split-prime factor Q_X(t)/(1−t)^r_X, and factor (1−t)^−1 at every inert prime and at 2. The primitive component B_X has split-prime factor Q_X(t)/(1−t)^(r_X−1), and factor 1 at every inert prime and at 2. Consequently A_X=ζ B_X on Re(s)>1, and

F(s)=[A_L−6A_V−3A_U+14A_W−6ζ]/48,
G(s)=[B_L−6B_V−3B_U+14B_W−6]/48.

The constant −6 in G is essential: the four component coefficients at e=1 sum to 6, whereas P(1)=0. The functions N and P themselves are not being asserted multiplicative. These are finite linear combinations of Euler products, not single Euler products for N and P.

## 2. Convergence, continuation, and poles

For Re(s)>1 all component Euler products and Dirichlet series converge absolutely: polynomial growth in each exponent gives local factors 1+O(p^−σ), uniformly away from σ=1, and the prime sum converges. The signed combinations converge absolutely as well.

For χ=χ_4,

A_X=ζ^(c_X+1)L(s,χ)^c_X H_X,
B_X=ζ^c_X L(s,χ)^c_X H_X,

where

H_X=(1−2^−s)^c_X × product over p≡1 mod 4 of Q_X(p^−s)(1−p^−s)^(2c_X+1−r_X) × product over p≡3 mod 4 of (1−p^−2s)^c_X.

This checks 2 separately: the leading H_2 is (1−2^−s)^10, not exponent 11. At inert primes both total and primitive factorizations give precisely the stated factors. Every split local factor of H has constant term 1 and zero linear term. Thus H converges locally uniformly on Re(s)>1/2, giving a holomorphic function there. Finite local factors may have complex zeros; no global nonvanishing claim is justified or needed. All factors at s=1 are strictly positive, and their deviations from 1 are summable, so H_X(1)>0.

Positive integer powers of ζ and L require neither logarithmic branches nor division by a possibly vanishing L-value. F and G therefore continue meromorphically throughout Re(s)>1/2 with no possible pole except 1. Their leading components have exact pole orders 11 and 10; every other component has smaller order. The common leading Laurent coefficient is

C=(π/4)^10 H_L(1)/48 >0,

H_L(1)=2^−10 product over p≡1 mod 4 of (1+18/p+5/p²)(1−1/p)^18 × product over p≡3 mod 4 of (1−1/p²)^10.

## 3. Tauberian normalization, checked against the original source

Source: Hubert Delange, Généralisation du théorème de Ikehara, Annales scientifiques de l'École Normale Supérieure, series 3, 71 (1954), 213–242; Theorem III, §5.2.1, printed pp. 235–236. https://www.numdam.org/item/ASENS_1954_3_71_3_213_0.pdf

I independently inspected both the extracted text and rendered pages 235–236. The excluded exponents in Theorem III are the nonpositive integers, not the positive integers. Its standing transform is the ordinary Laplace transform of a nonnegative nondecreasing function. Set α(t)=sum_{n≤exp(t)}N(n). Then its transform is F(s)/s for Re(s)>1. Apply Theorem III with a=1 and ω=11; use G(s)/s and ω=10 for primitive counts. Both transforms are holomorphic at all other points of Re(s)=1 and have the required local pole representation. Division by s leaves the leading coefficient C unchanged. Therefore

sum_{e≤x}N(e) ~ C x(log x)^10/10!,
sum_{e≤x}P(e) ~ C x(log x)^9/9!.

This uses nonnegativity of the actual candidate counts, not of their signed component combination. Lower Laurent terms can be absorbed into the holomorphic numerator multiplying the leading pole. The citation is a theorem application, not a claim of a power-saving error term.

## 4. Optional all-triple comparison

Let M(e)=product over split primes p^a||e of (2a+1). All distinct unordered gap triples number T(e)=(M³−9M²+23M−15)/48. Apply the same ancestry convention to define primitive T_prim. The leading cubic local series is

sum_{a≥0}(2a+1)^3 t^a=(1+23t+23t²+t³)/(1−t)^4.

Its c is 13. Put

K(s)=(1−2^−s)^13 product over split p of (1+23t+23t²+t³)(1−t)^23 product over inert p of (1−t²)^13,
C_all=(π/4)^13 K(1)/48>0.

The total and primitive all-triple poles are 14 and 13. Their summatory leading terms are C_all x(log x)^13/13! and C_all x(log x)^12/12!. Hence the aggregate candidate/all-triple ratios are asymptotic to

(C/C_all)(13!/10!)(log x)^−3,
(C/C_all)(12!/9!)(log x)^−3,

respectively. These are ratios of aggregate counts, not averages of per-center fractions. No actual Schur-count asymptotic follows.

## 5. Independent finite verification

check.py does not import or execute the producer's checker or use the balanced closed formula. It enumerates every signed-orbit box triple needed for all e≤5000, checking repeated absolute maxima and primitive maxima directly. It independently expands rational local Euler series by repeated geometric convolution and compares every coefficient to direct enumeration. It verifies balanced and all-triple divisor ancestry, both ramified/inert stripping controls, the exceptional unit coefficient, and the cubic all-triple series through exponent 100.

Normal and optimized Python runs pass and are byte-identical. No assertions are used. There are 11 distinct exponent boxes in this bound. Totals: sum N=6025, sum P=2592, sum T=10498, sum T_prim=5010. Selected controls: N(65)=P(65)=4; N(130)=N(195)=4 while P(130)=P(195)=0; N(1105)=150 and P(1105)=138.

Publishable audit artifacts are RESULT.md, check.py, certificate.json, check.out, and check-optimized.out. The downloaded complete paper and rendered source pages are private inspection artifacts only and must not be copied wholesale into publication.

## 6. Final producer manuscript review

After the independent derivation and checks, I read the completed producer RESULT.md in full. Its additional check H_W(s)=1/[(1+2^−s)ζ(2s)] is correct on Re(s)>1/2. The exact abscissa 1 follows from the positive leading singularity and absolute convergence to its right. The consequence sum P/sum N ~10/log x is correct. Its final section for actual signed gap relations is carefully limited to ancestry, domination, and absolute convergence: multiplication of gaps by a common positive square preserves and reflects the signed equation, so that structural convolution is valid under the explicitly fixed unordered-triple convention. It does not transfer the candidate Euler products or asymptotics to actual relations. No corrections are required.
