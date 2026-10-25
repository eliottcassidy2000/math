# Exact Dirichlet series for valuation-filtered Gaussian gap candidates

Status: PROVED from the previously proved gap-triple enumeration and disjoint primitive ancestry. The asymptotic passage uses Delange's original Tauberian theorem with its hypotheses verified below. These are counts of necessary valuation candidates, not counts of actual signed gap relations, magic hourglasses, or magic squares.

## 1. Objects and normalization

For an integer center root e≥1, write a_p=v_p(e) for its split primes p≡1 (mod 4). Define N(e) to be the number of unordered triples of distinct positive Gaussian gaps satisfying the repeated-absolute-maximum valuation filter. Let P(e) count primitive such triples. The Gaussian-box enumeration and unique maxima-vector ancestry give

    N(e) = [∏L(a_p) − 6∏V(a_p) − 3∏U(a_p) + 14∏W(a_p) − 6]/48,
    N(e) = Σ_{d|e} P(d),

where L(a)=12a²+8a+1, V(a)=2a²+4a+1, U(a)=4a+1, W(a)=2a+1. Products here and below indexed by split primes use just the primes dividing e unless otherwise stated. Empty products equal one, so N(1)=P(1)=0. P(e)=0 if 2 or any prime ≡3 (mod 4) divides e. For nontrivial split-only e,

    P(e) = [∏(24a_p−4) − 6∏(4a_p+2) − 3·4^k + 14·2^k]/48,

where k is the number of distinct prime factors. The e=1 exception cannot be omitted: this last raw expression gives 1/8 there.

Source for these exact combinatorial facts: predecessor `session58-tropical-gap-triples/RESULT.md`, §§2 and 8. They supply nonnegativity without relying on cancellation in the displayed formulas. The variable bounded in all summatory statements below is the center root e, not the central square e².

## 2. Euler products and exact heredity

For X∈{L,V,U,W}, let x_X(e)=∏X(a_p). Set

    A_X(s)=Σ_{e≥1} x_X(e)e^(−s).

Put t=p^(−s). Its split-prime Euler factors S_X(t) are

    S_L(t)=(1+18t+5t²)/(1−t)^3,
    S_V(t)=(1+4t−t²)/(1−t)^3,
    S_U(t)=(1+3t)/(1−t)^2,
    S_W(t)=(1+t)/(1−t)^2.

These follow by summing the quadratic or linear polynomial X(a) over a≥0. At every inert prime and at 2, all four Euler factors equal (1−t)^(−1), since those prime exponents do not change x_X(e).

Each local factor is 1+O(p^(−σ)) for σ>1, uniformly outside finitely many primes, and each local power series converges absolutely there. Therefore A_X has an absolutely convergent Dirichlet series and Euler product for Re(s)>1. In that half-plane,

    F(s):=Σ N(e)e^(−s)
        =[A_L(s)−6A_V(s)−3A_U(s)+14A_W(s)−6ζ(s)]/48.

Absolute convergence follows, for example, by the triangle inequality against the five convergent series above. F itself is NOT an Euler product: N is not multiplicative and N(1)=0.

Define B_X(s)=A_X(s)/ζ(s). Its Euler product has factor 1 at 2 and inert primes, and split factors

    B_L,p=(1+18t+5t²)/(1−t)^2,
    B_V,p=(1+4t−t²)/(1−t)^2,
    B_U,p=(1+3t)/(1−t),
    B_W,p=(1+t)/(1−t).

Their positive-exponent coefficients are respectively 24a−4, 4a+2, 4, and 2, and the constant coefficient is one. Consequently these also converge absolutely for Re(s)>1. Exact divisor ancestry, or the coefficient formulas including e=1, gives

    G(s):=Σ P(e)e^(−s)
        =[B_L(s)−6B_V(s)−3B_U(s)+14B_W(s)−6]/48,
    F(s)=ζ(s)G(s).

The last subtraction in G is the constant Dirichlet series −6, supported at e=1, NOT −6ζ(s). This corrects the primitive empty-support anomaly. The relation F=ζG is an exact arithmetic convolution identity, not an analogy with another generating function or a claim that G itself factors multiplicatively.

## 3. Zeta and character factorization

Let χ=χ_4 be the real Dirichlet character modulo 4: χ(2)=0, χ(p)=1 for split odd primes, and χ(p)=−1 for inert primes. To avoid confusion, write D(s)=L(s,χ_4) for the Dirichlet L-function, keeping L(a) for the quadratic row-count polynomial.

For X=L,V,U,W respectively set

    c_X = 10,3,2,1;
    r_X = 3,3,2,2;
    Q_X(t) = 1+18t+5t², 1+4t−t², 1+3t, 1+t.

Then

    A_X(s)=ζ(s)^(c_X+1) D(s)^c_X H_X(s),
    B_X(s)=ζ(s)^c_X D(s)^c_X H_X(s),

where

    H_X(s)=(1−2^(−s))^c_X
       × ∏_{p≡1 mod4} Q_X(p^(−s))(1−p^(−s))^(2c_X+1−r_X)
       × ∏_{p≡3 mod4} (1−p^(−2s))^c_X.

Verification is prime-by-prime. At a split prime ζ^(c+1)D^c contributes (1−t)^(−(2c+1)), leaving the required Q/(1−t)^r. At an inert prime it contributes (1−t)^(−(c+1))(1+t)^(−c), and multiplication by (1−t²)^c leaves (1−t)^(−1). At 2 only ζ contributes; multiplying its factor (1−t)^(−(c+1)) by (1−t)^c leaves exactly (1−t)^(−1). Thus the special prime-2 normalization is essential.

For split primes the exponent 2c+1−r is respectively 18,4,3,1, exactly the linear coefficient of Q. Hence each split H factor is 1+O(|t|²). The inert H factors have the same bound. The products therefore converge absolutely and locally uniformly in Re(s)>1/2, proving that H_X is holomorphic there. This statement allows zeros of finitely many local factors; no globally nonvanishing logarithm is asserted. At real s=1 all factors are positive, and absolute convergence of their deviations from one proves H_X(1)>0.

In particular the leading factor is

    H(s):=H_L(s)
      =(1−2^(−s))^10
       ∏_{p≡1 mod4}(1+18p^(−s)+5p^(−2s))(1−p^(−s))^18
       ∏_{p≡3 mod4}(1−p^(−2s))^10.

As a check, H_W(s)=1/[(1+2^(−s))ζ(2s)] in this half-plane. D is entire, ζ has only its simple pole at 1, and all exponents used are nonnegative integers. Thus the displayed products extend A_X and B_X meromorphically throughout Re(s)>1/2 with no poles except possibly at 1. Character zeros create no poles and require no zero-free-region argument.

At s=1 the A_X pole orders are exactly 11,4,3,2; the remaining ζ term has order 1. The B_X pole orders are exactly 10,3,2,1, and the final constant has no pole. Since D(1)=1−1/3+1/5−⋯=π/4 and H(1)>0, no lower-order terms can cancel the leading pole. Both F and G have the same leading Laurent coefficient

    C = H(1)(π/4)^10/48 > 0,

at orders 11 and 10 respectively. In particular F and G have abscissa of absolute convergence exactly 1, not merely at most 1.

## 4. Verified Tauberian passage

Primary source: Hubert Delange, “Généralisation du théorème de Ikehara”, Annales scientifiques de l'École Normale Supérieure, 3e série, 71 (1954), 213–242; Theorem III, §5.2.1, printed pp.235–236 (PDF pages 23–24, zero-based).

- Bibliographic record and DOI: https://archive.numdam.org/articles/10.24033/asens.1023/
- Original paper: https://www.numdam.org/item/ASENS_1954_3_71_3_213_0.pdf
- Exact source pin: Theorem III, local form (11), conclusion immediately following (12); the Laplace-transform and monotonicity convention is in §1, printed p.213.

The integer-pole specialization needed here says that if a nonnegative nondecreasing α(t) has Laplace transform f(s), holomorphic across Re(s)=1 except at 1, and f(s)=(s−1)^(−m)g(s) locally with g holomorphic and g(1)>0, then α(t)~g(1)e^t t^(m−1)/Γ(m). Positive integers m are allowed; the excluded integers in the original theorem are nonpositive.

To apply it, take α_N(t)=Σ_{e≤exp(t)}N(e). This is nonnegative and nondecreasing, and for Re(s)>1 Tonelli followed by absolute convergence gives

    ∫_0^∞ exp(−st) α_N(t) dt = F(s)/s.

For complex s the identity follows by absolute convergence, and for positive real s by direct integration of each indicator t≥log(e). The same identity holds with P,G. The factor 1/s is holomorphic and equals one at s=1. The established continuation to Re(s)>1/2 verifies holomorphy at every other point of the boundary line. Locally take g(s)=(s−1)^m F(s)/s or (s−1)^m G(s)/s; these are holomorphic, with value C>0. Thus every hypothesis is checked, and Delange gives unconditionally

    Σ_{e≤x} N(e) ~ C x(log x)^10/10!,
    Σ_{e≤x} P(e) ~ C x(log x)^9/9!.

Only these leading equivalents are claimed. No remainder estimate or full lower-order polynomial is inferred from the pole order alone. The factor C is the leading Laurent coefficient, not the ordinary residue of the higher-order pole. A further immediate consequence is ΣP/ΣN ~ 10/log x; this is an aggregate candidate-weighted primitiveness ratio.

## 5. Comparison against all unordered gap triples

This is an additional counting corollary, not an existence theorem. Put w(e)=∏(2a_p+1). There are (w(e)−1)/2 distinct positive gaps, so the count of all unordered triples of distinct gaps is

    T(e)=binom((w(e)−1)/2,3)
        =[w(e)^3−9w(e)^2+23w(e)−15]/48.

Let R(e) count primitive triples without the valuation filter. The same unique maximum-vector partition gives T(e)=Σ_{d|e}R(d); R vanishes at centers divisible by 2 or an inert prime. This uses only the Gaussian box and is valid without any balance predicate.

The leading multiplicative function w(e)^3 has split local series

    Σ_{a≥0}(2a+1)^3 t^a = (1+23t+23t²+t³)/(1−t)^4.

Its root Dirichlet series factors as ζ(s)^14 D(s)^13 J(s), with

    J(s)=(1−2^(−s))^13
       ∏_{p≡1 mod4}(1+23p^(−s)+23p^(−2s)+p^(−3s))(1−p^(−s))^23
       ∏_{p≡3 mod4}(1−p^(−2s))^13.

Again J is holomorphic with an absolutely convergent product for Re(s)>1/2, and J(1)>0. For w² the corresponding split series is (1+6t+t²)/(1−t)^3, giving root pole order 5, and for w the order is 2. The constant sequence contributes ζ. Dividing by ζ for primitive ancestry reduces all positive pole orders by one and replaces the constant-sequence Dirichlet series by 1. Nonnegativity of T and R and the same continuation verify the identical Delange argument. With C_all=J(1)(π/4)^13/48,

    Σ_{e≤x}T(e) ~ C_all x(log x)^13/13!,
    Σ_{e≤x}R(e) ~ C_all x(log x)^12/12!.

Consequently the two surviving fractions satisfy

    (ΣN)/(ΣT) ~ [13!/10!] [H(1)/J(1)](π/4)^(−3)(log x)^(−3),
    (ΣP)/(ΣR) ~ [12!/9!]  [H(1)/J(1)](π/4)^(−3)(log x)^(−3).

The factorial ratios are 1716 and 1320. These are ratios of aggregate counts, not averages of the fraction at each center. The filter preserves actual signed gap relations but discards their numerical additive test; no density or existence result for actual hourglasses follows except the corresponding necessary upper bound.

## 6. Exact finite validation and hostile controls

`check.py` uses only the Python standard library and no assertions. Both ordinary and `python -O` runs produce identical output. Through every e≤100000 it checks nonnegative integral N,P, the full divisor ancestry, and independent Euler-factor coefficient reconstruction. Rational local series are checked through exponent 100, and cancellation of every linear H term is checked exactly.

Selected coefficients (e; N(e); P(e)):

    1;0;0       5;0;0       25;0;0
    65;4;4      130;4;0     195;4;0
    325;15;11   845;15;11   1105;150;138
    4225;54;28  71825;1592;730

The first nonzero coefficient is e=65. Through 100000 the exact sums are ΣN=515553 and ΣP=178798. The 130 and 195 cases explicitly test removal of prime 2 and inert scaling in the primitive series. The constant subtraction is checked at e=1. These finite identities validate the independently proved algebra; they do not prove an asymptotic numerically or assert that the asymptotic is already accurate at this cutoff.

Artifacts: `check.py`, `certificate.json`, `check.out`, `check-optimized.out`. No complete third-party article has been copied into this deliverable.

## 7. What remains true for actual additive relations

Fix the same counting convention: each unordered triple of distinct positive Gaussian gaps is counted once if it admits a signed relation ±d_1±d_2±d_3=0. Denote the actual-relation count by A(e), and its primitive subset by A_prim(e). This does not count oriented hourglass arrangements with extra multiplicities.

The unique common-root ancestry used above also restricts to this subset. Scaling all associated roots by g multiplies each gap by g², so it preserves and reflects the signed equation. Dividing by the unique common root divisor preserves the unordered, distinct, positive triple convention. Therefore

    A(e)=Σ_{d|e}A_prim(d),
    0≤A(e)≤N(e),   0≤A_prim(e)≤P(e).

Their Dirichlet series converge absolutely for Re(s)>1 by domination, and satisfy the exact relation F_A(s)=ζ(s)G_A(s) there. No Euler product, continuation across Re(s)=1, leading pole, or leading asymptotic for F_A or G_A is established. This separates a structural ancestry identity valid for actual relations from the explicit coefficient formulas proved only for the necessary valuation filter.
