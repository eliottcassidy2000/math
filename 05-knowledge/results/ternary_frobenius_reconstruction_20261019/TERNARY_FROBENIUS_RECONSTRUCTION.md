# Frobenius defects in ternary period reconstruction

Status: original proved algebraic results and separately labelled finite-exact computations. Independently reviewed, 2026-10-19; see `review/PERIOD_PROOF_REVIEW.md`.

## 1. Source, scope, and inheritance

Vasily Ilin and Ken Ono, *Ternary polynomial relations and periods of cusp forms*, arXiv:2610.19309v1 (13 October 2026), https://arxiv.org/pdf/2610.19309.

We use only the following definitions from the paper. On polynomials over a field k, put

- C f(x,y,z) = f(y,z,x);
- A f(x,y,z) = f(-y,x-y,z), so A^3 = I;
- N = I + A + A^2;
- W_n(k) = {homogeneous f of degree n: Cf=f, Nf=0};
- U_n(k) = k[x,y]_n + k[y,z]_n + k[z,x]_n;
- rho(f) = f(x,y,0), and Phi(q)=q(x,y)+q(y,z)+q(z,x).

CITED: Theorem 1.1 and Corollary 1.3 of Ilin–Ono say, in characteristic zero, W_n is contained in U_n and restriction identifies W_n with normalized binary period polynomials, with inverse Phi. Odd degrees vanish. Their Section 7.4 explicitly shows that positive characteristic is different: a GL_3(F_p)-invariant product times a degree-ten solution gives an all-interior solution of degree p^3+9, and an explicit degree-eleven example works at p=2. Their Lemma 5.3 gives a complete characteristic-two symmetric kernel as a rank-eight free module over the three Dickson invariants.

The following results are additional elementary derivations, not claims made by that paper. They distinguish two different failures: odd supported solutions already occur in degree p, whereas a nonzero kernel of rho means genuine loss of reconstruction.

Repository inheritance: `01-canon/theorems/THM-4067-seminormal-period-kernel-and-figure-eight-completeness-obstruction.md` separates period cancellation from common-target reconstruction. The map here is literal polynomial restriction; its missing sidecar is ker rho. This is a precise parallel in the logic of reconstruction, not a map to Keller pairs and not progress on the Jacobian conjecture. The hostile example is an all-interior solution invisible on every coordinate plane. The underused operation is reduction modulo p followed by an explicit additive-polynomial lift.

Live concept board: norm kernels; boundary restriction; additive/Frobenius polynomials; rational-to-polynomial descent; integral lifting. No tournament is imposed on the polynomial coefficients.

## 2. Complete odd supported classification

PROVED. Let p>=5 be prime, k a field of characteristic p, and n>=1 odd. Then

    dim_k(W_n(k) intersect U_n(k)) = floor(log_p n).

Here floor(log_p 1)=0. More explicitly, for a power a=p^r<n put

    Q_a(X,Y) = X^(n-a) Y^a - X^a Y^(n-a).

Every supported solution is alternating and uniquely has the form

    f = Phi(sum_{a=p^r<n} c_a Q_a).

If n is a power of p, the coefficients c_a are arbitrary. If n is not a power of p, their sum must be zero. Consequently:

- if n=p^h, a basis is Phi(Q_1), Phi(Q_p), ..., Phi(Q_(p^(h-1)));
- otherwise, writing h=floor(log_p n), a basis is Phi(Q_(p^r)-Q_1), r=1,...,h.

The non-power case with h=0 has an empty basis.

### Proof

Since 3 is invertible, Nf(0,0,z)=3f(0,0,z) shows that all pure powers vanish. Every remaining monomial of f in U_n uses exactly two variables. Cyclicity therefore gives f=Phi(q) uniquely, with q=f(X,Y,0) zero on both axes. Restriction of Nf=0 implies N_2 q=0.

Because 2 is invertible and exchanging x,y preserves both defining equations, split f into symmetric and alternating parts, with corresponding binary symmetry q(v,u)=epsilon q(u,v), epsilon=+1 or -1. After the three terms of N_2q(x,y) are removed, direct expansion gives

    N Phi(q) = r(y,z) + epsilon r(x,z) + r(x-y,z),
    r(u,z) = q(u,z) + epsilon q(-u,z).

For epsilon=+1, putting x=y gives 2r(x,z)=0, so q is odd in each variable. A symmetric q of odd total degree with that property is zero.

For epsilon=-1, the residual equation says that r(u,z) is additive in u. A polynomial additive in one variable in characteristic p is a sum of powers u^(p^r), with coefficients in k[z]. This follows directly from the binomial coefficients, or Lucas's digit criterion. Homogeneity and odd n make each complementary exponent n-p^r even. Alternation then forces exactly the displayed linear combination of Q_a. The terms are independent: the odd exponent a distinguishes them, while n-a is even.

For a=p^r and n-a even, Frobenius gives (X-Y)^a=X^a-Y^a. A direct expansion now yields the same norm for every Q_a:

    N_2 Q_a = X^n - Y^n - (X-Y)^n =: F_n(X,Y).

Furthermore, the residual ternary norm is zero because its odd part is additive. Thus N Phi(sum c_a Q_a)=(sum c_a)F_n. Finally, F_n=0 exactly when n is a power of p, again by the binomial/Lucas criterion. Counting available powers a<n gives floor(log_p n) in both cases. This proves necessity, sufficiency, uniqueness, and the dimension formula.

### First example and precise failure

For every p>=5,

    f_p = Phi(X^(p-1)Y - XY^(p-1))

is a nonzero alternating degree-p solution. Its integral norm is

    X^p - Y^p - (X-Y)^p,

which vanishes modulo p. This solution still has two-variable support and is recovered by Phi from its restriction. It defeats odd-degree vanishing and the normalized two-term binary-period condition, but it does NOT itself defeat injectivity of rho. The distinction matters.

## 3. A six-term rational construction from two additive polynomials

PROVED. Let k have odd characteristic, and let L,M in k[T] be additive polynomials with zero constant term. Define

    H_(L,M)(x,y,z) = det [[x^2,L(x),M(x)],
                          [y^2,L(y),M(y)],
                          [z^2,L(z),M(z)]] / (xyz).

This is a polynomial: each entry in every row is divisible by that row's variable. It is alternating, hence cyclic, and N H_(L,M)=0.

For the norm identity, write D(u,v)=L(u)M(v)-L(v)M(u). This is alternating and additive in each variable. Over the rational function field,

    H = x D(y,z)/(yz) + y D(z,x)/(zx) + z D(x,y)/(xy).

Put P=D(x,y), Q=D(y,z), R=D(z,x). Additivity gives

    A H = y(R+Q)/((x-y)z) - (x-y)Q/(yz) - zP/(y(x-y)),
    A^2 H = (x-y)R/(xz) - x(Q+R)/((x-y)z) + zP/(x(x-y)).

The coefficients of P,Q,R each sum to zero in H+AH+A^2H. The rational identity is therefore a polynomial identity after the established divisibility.

Taking L(T)=T^(p^a), M(T)=T^(p^b), a,b>=1 distinct, gives the six-term alternant

    Alt(1,p^a-1,p^b-1),

of degree p^a+p^b-1, where Alt(i,j,k)=det[t_r^(i),t_r^(j),t_r^(k)] for (t_1,t_2,t_3)=(x,y,z). All exponents are positive, so rho(H)=0. It is nonzero because the three exponents are distinct. In particular, every odd prime has a nonzero all-interior solution of degree p^2+p-1.

The same rational calculation also works in characteristic two, where signs coincide. Taking L=T^4, M=T^8 yields Alt(1,3,7), the characteristic-two degree-eleven example (equal to xyz D^2 from Section 7.4 of Ilin–Ono). At p=3, Alt(1,2,8) is a degree-eleven example.

## 4. Linear-in-p zero-restriction construction

PROVED. For every prime p>=5, define, over k of characteristic p,

    2 G_p = Alt(1,3,3p-1)
            - Alt(1,p+1,2p+1)
            - 2 Alt(1,p+2,2p)
            - Alt(2,p,2p+1)
            - Alt(2,p+1,2p)
            + Alt(3,p,2p)
            - 2 Alt(p,p+1,p+2).

Then G_p is nonzero, homogeneous of degree 3p+3, alternating, and satisfies N G_p=0 and rho(G_p)=0. Thus the restriction map fails to be injective by degree 3p+3 for every prime p>=5. This is an upper bound, not a universal minimum.

### A universal rational certificate

Introduce independent variables (x,y,z,X,Y,Z). For a pair (t,T), write m_(a,b)(t,T)=t^a T^b, allowing a=-1. Define an alternant of three such functions by evaluating them at (x,X),(y,Y),(z,Z). Let F be half of

    Alt(m10,m30,m(-1,3))
    - Alt(m10,m11,m12)
    - 2 Alt(m10,m21,m02)
    - Alt(m20,m01,m12)
    - Alt(m20,m11,m02)
    + Alt(m30,m01,m02)
    - 2 Alt(m01,m11,m21).

Let A act simultaneously on both triples:

    (x,y,z,X,Y,Z) -> (-y,x-y,z,-Y,X-Y,Z).

Then the rational identity (I+A+A^2)F=0 holds over Q. An entirely integral certificate is obtained by putting P=2xyz F:

    (x-y)P - x(A P) + y(A^2 P) = 0.                 (* )

This is a finite six-variable polynomial identity with integer coefficients. The primitive certificate stored in `rational_kernel_d3_k3_eps-1.json` is the negative of this displayed P; the overall sign has no effect on the norm identity, and the verifier checks the sign explicitly. The accompanying standalone verifier constructs P from the seven displayed alternants, expands the three linear substitutions, and checks that every coefficient of (* ) vanishes. Thus the identity does not depend on a finite range of primes or degrees.

Substitute X=x^p,Y=y^p,Z=z^p. Frobenius makes the six-variable action compatible with the original three-variable action. The only denominator in F is a coordinate variable in the term T^3/t, which becomes t^(3p-1); hence the specialization is polynomial and is exactly G_p. Every exponent in every alternant is positive, so all coordinate-plane restrictions vanish.

Nonzero: for p>=5 the seven unordered exponent triples in the formula are mutually distinct, and each contains three distinct exponents. The coefficient of x y^3 z^(3p-1) in 2G_p is 1. No cancellation can remove it. Homogeneity and alternation are immediate.

The same construction can be substituted with any odd-characteristic additive polynomial L: set (X,Y,Z)=(L(x),L(y),L(z)). It again gives a polynomial solution divisible by xyz. It may be zero for special L; for L(T)=T^p, p>=5, the preceding coefficient proves nonzero. Adding a linear term cT to L does not change the construction, as a further universal identity verified by the certificate.

## 5. A symmetric degree-(2p+12) family

PROVED by an exact universal polynomial certificate. Define the weighted homogeneous degree-seven polynomial H(x;s,t), with weights deg x=deg s=1 and deg t=2, by

    H = 175 s x^6
      + (-70 s^2 + 315 t) x^5
      + (-175 s^3 + 525 s t) x^4
      + (136 s^4 - 754 s^2 t + 639 t^2) x^3
      + (21 s^5 - 14 s^3 t - 231 s t^2) x^2
      + (-34 s^6 + 263 s^4 t - 519 s^2 t^2 + 108 t^3) x
      + 11 s^7 - 111 s^5 t + 287 s^3 t^2 - 68 s t^3.

For a prime p put

    S_p = xyz sum_cyclic [
      (y^2-z^2)^2 H(x;y+z,yz)
      (y^(p-1)-x^(p-1))(z^(p-1)-x^(p-1)) ].

For every prime p>=11, S_p is nonzero, homogeneous of degree 2p+12, symmetric, and satisfies N S_p=0. Every monomial uses all three variables, so rho(S_p)=0.

The displayed polynomial is identically zero for p=3,5,7. This is a degeneration of this construction, not absence of other solutions at those primes.

### Universal proof and transparent representation

Set B_xy=xY-yX, B_xz=xZ-zX and let

    F = sum_cyclic [(y^2-z^2)^2 H(x;y+z,yz) B_xy B_xz / x].

This rational function is homogeneous of bidegree (12,2) in the independent triples (x,y,z) and (X,Y,Z). It is symmetric under simultaneous coordinate permutation and is invariant under the gauge change

    (X,Y,Z) -> (X,Y,Z) + c(x,y,z).

Put P=xyz F. Then P is the 363-term integral polynomial saved in `rational_kernel_d12_k2_eps1.json`. The independently implemented verifier establishes both its equality with the compact formula above and the exact identity

    (x-y)P - x(A P) + y(A^2 P) = 0,

where A acts simultaneously on both triples. Division in the rational function field proves N F=0. Frobenius specialization X=x^p,Y=y^p,Z=z^p gives exactly the displayed S_p, which is manifestly a polynomial divisible by xyz.

There is also a useful matrix description. Write

    P = R_x X^2 + R_y Y^2 + R_z Z^2
        + B_xy^coefficient XY + B_xz^coefficient XZ
        + B_yz^coefficient YZ,

where the superscript distinguishes coefficient names from wedge expressions. Gauge invariance forces

    B_xy^coefficient = (z^2 R_z - x^2 R_x - y^2 R_y)/(xy),

and cyclic analogues. Thus one diagonal coefficient determines the whole symmetric quadratic form. Here

    R_x = y^2 z^2 (y^2-z^2)^2 H(x;y+z,yz).

This is the exact representation-theoretic content being used: a quadratic form annihilating the original-vector line and therefore descending to the auxiliary-vector quotient by that line. A literal identification with a cusp-form representation has NOT been established.

### All-prime nonvanishing, with explicit coordinate bounds

Every monomial of F has exponent a_i in [-1,12] on each original variable, and exponent b_i in {0,1,2} on each auxiliary variable, with sum b_i=2. For p>=17, the encoding

    (a_i,b_i) -> a_i+p b_i

is injective: if a_i-a'_i is a multiple of p, its absolute value is at most 13<p, so it is zero. The coefficient of

    x^1 y^12 z^(-1) Z^2

in F is exactly 11. Hence its specialized monomial x y^12 z^(2p-1) has nonzero coefficient 11 modulo every prime p>=17. The two boundary cases are checked explicitly by the same exact coefficient dictionary:

- p=13: coefficient of x y^12 z^25 is 11 modulo 13;
- p=11: coefficient of x y^6 z^27 is 10 modulo 11.

This proves nonzero for every p>=11; it does not extrapolate from a finite prime list.

### Why this family was found, and what is minimal

The p=11,13,17,19 first-kernel witnesses compress into this common rational bidegree-(12,2) representation. A finite exact rational ansatz searched all P with denominator xyz, the stated bidegrees and symmetry, and the condition that a variable absent from the auxiliary monomial has original numerator exponent at least two. This is precisely the condition ensuring every Frobenius specialization has all-interior support.

For auxiliary degree two and original rational degree d=0,...,11, the ansatz has zero kernel over Q in both symmetry types. At d=12 its symmetric kernel has dimension one and its alternating kernel remains zero. At auxiliary degree three, the first solution through d=3 is the alternating bidegree-(3,3) family in Section 4, again with dimension one. The statements about these bounded ansatz spaces are FINITE-EXACT: full column ranks modulo 100003 certify absence over Q, and the integral kernel vectors certify the one-dimensional cases. They are not lower bounds on every modular polynomial.

A separate linear-auxiliary search finds no solutions for either symmetry and d=0,...,30. This is a bounded negative control, not an all-degree vanishing theorem.

## 6. Finite-exact threshold census

The universe is all homogeneous cyclic polynomials of the stated degree over F_p, not only symmetric or alternating candidates. A cyclic-monomial orbit basis contains each distinct monomial once, including size-one orbits. The norm matrix is expanded exactly with integral binomial coefficients, then reduced modulo p. Restriction-kernel matrices retain exactly the orbits with all three exponents positive. Their nullity is dim ker(rho|W_n).

Initial exact minima:

| p | first degree with nonzero ker rho | nullity there | certified lower range |
|---|---:|---:|---|
| 2 | 11 | 1 | every n=0,...,10 |
| 3 | 11 | 1 | every n=0,...,10 |
| 5 | 18 | 1 | every n=0,...,17 |
| 7 | 24 | 1 | every n=0,...,23 |
| 11 | 34 | 1 | every n=0,...,33 |
| 13 | 38 | 1 | every n=0,...,37 |
| 17 | 46 | 1 | every n=0,...,45 |
| 19 | 50 | 1 | every n=0,...,49 |

For p=5,7 the first kernel is spanned by the degree-3p+3 construction. For p=11, degree 34 is strictly earlier than its degree-36 instance. The first kernels at p=11,13,17,19 are exactly the spans of S_p. No formula for the general minimum is proved. The p=11 first witness is symmetric and much denser, so the alternating sparse construction is not a candidate for the whole minimum.

The initial production census also covers all n=0,...,36 for p=2,3,5,7,11,13. A fully independent standard-library implementation re-expands the substitutions literally and uses sparse column elimination; it verifies every degree up to each displayed minimum, and both universal certificates. Its checks run identically under Python -O. A separate supported-space test matches Section 2 in every odd degree through 37 for p=5,7,11,13. These tests validate implementations; Section 2 is proved in all degrees.

### A decisive exceptional prime: 29

FINITE-EXACT COUNTEREXAMPLE. The proposed general formula

    delta_p = min(3p+3, 2p+12)

is false. At p=29 there is an alternating all-interior solution of degree 42, well below 70. The explicit certificate and an independent literal-substitution verifier are in `EXCEPTIONAL_PRIME_29.md` and `verify_exception_29.py`. Its 138 terms group into 23 alternants, every exponent is even, and its norm is identically zero modulo 29. This is a polynomial identity, not a test on finite-field points.

Consequently the universal constructions prove upper bounds, not a universal equality. The p=29 phenomenon requires an additional arithmetic explanation. The open problem is now to classify the exceptional primes and degrees below those bounds, not to prove the disproved formula. A finite extension of the production census finds first degrees 58,42,74 for p=23,29,31 respectively; these new minimality claims require the separate independent review recorded with the package.

## 7. Connection contract and limits

Source: integer norm-relation matrices and characteristic-zero restriction theorem.
Target: norm kernels over F_p and explicit additive-polynomial constructions.
Map: coefficient reduction; additionally, universal rational identities specialize through Frobenius.
Preserved: literal cyclicity, norm vanishing, homogeneity in the monomial cases, and coordinate-plane invisibility where stated.
Destroyed: odd-degree vanishing, two-term period normalization, and injective reconstruction. Polynomial function evaluation on F_p^3 also destroys nonzero polynomial information; these identities are in polynomial rings, not functions on finite point sets. The characteristic-zero detector partial_x partial_y partial_z must also not be imported: in characteristic p it misses interior monomials whenever an exponent is divisible by p. The computations use coefficient support, not that derivative test. Hasse coefficient operators retain the missing p-adic exponent digits.
Sidecar: retain ker rho and the additive-polynomial exponents. For lifting back to characteristic zero, retain integral divisibility/lift obstructions; mod-p existence does not lift automatically.
Decisive hostiles: the degree-34 symmetric kernel at p=11 defeats a universal-minimum claim for the degree-3p+3 family; the degree-42 alternating kernel at p=29 defeats equality with the improved two-family upper bound. The paper's all-degree mod-four obstruction for odd alternating integral lifts, Theorem 5.2, is another explicit warning against lifting modular solutions naively.

No result here produces a characteristic-zero counterexample to Ilin–Ono, a modular cusp form from a modular polynomial alone, a Keller map, a conclusion about JC(2), or a proof of a cellular-automaton conjecture.


## 8. Reproduction and artifacts

From this package directory, run:

    python3 -B probe_norm.py
    python3 -B explore_supported.py
    python3 -B extend_thresholds.py
    python3 -B rational_kernel.py
    python3 -B verify_periods_independent.py
    python3 -B -O verify_periods_independent.py

The production search needs NumPy and SymPy. The independent verifier uses only the Python standard library and imports no production algebra or matrix code. Output files live beside the scripts. `modular_census.json`, `extended_thresholds.json`, and `independent_verification.json` include exact universes, dimensions, and witnesses. The two rational-kernel JSON files are complete integer coefficient certificates, not sampled evaluations.

The standalone theorem proofs are Section 2 and Section 3. Sections 4 and 5 use explicitly specified finite polynomial identities verified by coefficient expansion, which proves their all-prime specializations. The threshold table and bounded rational-ansatz uniqueness are finite-exact results.
