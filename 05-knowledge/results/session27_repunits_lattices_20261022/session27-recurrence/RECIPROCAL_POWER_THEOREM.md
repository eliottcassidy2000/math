# Reciprocal-cubic power compression: exact universal theorem and recovered seven-check proof

Status: proved algebraic synthesis, independently verified. No novelty claim. In particular, the seven-check Ramanujan application is explicitly published in Hirschhorn (1996), not a new improvement discovered in this session. No modular-form assertion is made.

## 1. Statement over arbitrary commutative rings

Let R be any commutative ring with identity, t an element of R, and d a nonnegative integer. Write E for the forward shift. Suppose finitely many R-valued sequences u^(j) satisfy

C_t(E)u^(j)=0, where C_t(x)=(x+1)(x²−tx+1)=x³+(1−t)x²+(1−t)x+1.

Thus u_(n+3)=(t−1)u_(n+2)+(t−1)u_(n+1)−u_n. For any homogeneous polynomial P of degree d with coefficients in R, the sequence v_n=P(u_n^(1),…,u_n^(m)) is annihilated by

F_(d,t)(x)=(x−(−1)^d) ∏_(r=1)^d [x²−(−1)^(d−r)L_r(t)x+1],

where L_0(t)=2, L_1(t)=t, and L_r(t)=tL_(r−1)(t)−L_(r−2)(t).

The degree is 2d+1; the leading coefficient is 1 and constant coefficient is (−1)^(d+1), a unit in every such R. The theorem includes integer sequences, every modulus (prime or composite), characteristic 2, rings with zero divisors or nilpotents, and t=±2. No denominators are inverted in the conclusion.

For sequences initially given on n≥0, the recurrence has a unique extension to all integers because its constant coefficient is 1. The resulting v has the same property. Any 2d+1 consecutive zero values force v to vanish at every integer index. Similarly, two sequences known to satisfy F agree everywhere if they agree on 2d+1 consecutive indices. This is a sufficient universal test, not a claim that every observable has that minimal order.

The coefficients of P must be independent of n, and P must be homogeneous of the stipulated degree. Fixed shifts of the input sequences are allowed, since they still satisfy C_t. Different t parameters or arbitrary nonhomogeneous expressions are not covered by this single F.

## 2. Proof, including repeated-root specializations

Work first over A=Z[T], with T an indeterminate, and set

C(T) = [[0,1,0],[0,0,1],[−1,T−1,T−1]].

For d≥1 let M_d(T)=C(T) tensor … tensor C(T), with d factors; for d=0 use the 1×1 matrix [1]. Its entries lie in A. We prove the polynomial-matrix identity

F_(d,T)(M_d(T))=0 in Mat_(3^d)(A).

Embed A into Q(T), and then into a splitting field of x²−Tx+1. Its roots are λ and λ^(−1). Together with −1 these are distinct over this function field, so C(T), and hence M_d(T), are diagonalizable. Each tensor eigenvalue has the form

(−1)^a λ^(b−c), with a+b+c=d and a,b,c≥0.

Put r=b−c. Then a≡d−r mod 2, so this eigenvalue is μ_r=(−1)^(d−r)λ^r, where −d≤r≤d. Conversely every r in this interval occurs: choose b=max(r,0), c=max(−r,0), a=d−|r|. The equality of the signs follows because r and |r| have the same parity.

The central eigenvalue is μ_0=(−1)^d. The pair μ_r, μ_(−r) has product 1 and sum (−1)^(d−r)(λ^r+λ^(−r)). The last parenthesis is exactly L_r(T), by its initial values and recurrence. Thus F is the product of one linear and d quadratic factors vanishing on every eigenvalue of M_d. Diagonalizability proves F(M_d)=0 over the splitting field.

Crucially, each entry of F(M_d) already belongs to Z[T]. Injectivity of Z[T] into that field means every such entry is the zero polynomial. Apply the unique ring map Z[T]→R sending T to t. The same matrix identity therefore holds over R. This proves the repeated-root and positive-characteristic cases without incorrectly diagonalizing a specialized matrix. It also proves validity when t itself is nilpotent or a zero divisor.

For a degree-d monomial in the given sequences, take the tensor product of their d state vectors (allowing repeated factors). Their simultaneous evolution is M_d(t), and their pointwise product is obtained by applying the tensor product of the first-coordinate functionals. The matrix identity gives F(E)=0 on this product. Linearity gives the result for P. Constants give the d=0 case.

The companion matrix has determinant −1, so its inverse exists integrally and after every specialization. The same holds for its tensor powers. Both the original and resulting recurrences therefore work backwards. Monicity propagates a zero block forward; the unit constant term propagates it backwards. This completes the proof.

### Degeneracy checks and why multiplicities cannot be discarded

At t=−2, L_r(−2)=2(−1)^r, hence F_(d,−2)(x)=(x−(−1)^d)^(2d+1).

At t=2, L_r(2)=2, hence F_(d,2)(x)=(x−1)^(d+1)(x+1)^d.

These multiplicities are part of the specialization. Replacing F by the product of its distinct specialized factors would generally be invalid because the companion matrix is not diagonalizable there. No minimal-order claim at special t or in positive characteristic is needed here.

## 3. Sharpness of the universal degree in characteristic zero

Take real t with |t|>2 and choose the real quadratic root λ with |λ|>1. The 2d+1 numbers μ_r above are distinct: their absolute values |λ|^r are distinct as r varies. The sequences e_r(n)=μ_r^n are linearly independent by the Vandermonde determinant on n=0,…,2d.

Every e_r is a homogeneous degree-d monomial in the three solutions (−1)^n, λ^n, λ^(−n), by the explicit choice of a,b,c above. Every degree-d observable is a linear combination of these e_r. Consequently the vector space of all degree-d observables has dimension exactly 2d+1, with shift eigenvalues μ_r. A common annihilator must vanish at every μ_r and therefore has degree at least 2d+1. F attains that bound.

For a rational t with |t|>2 the same minimal-common-annihilator statement holds over Q by extension of scalars: a smaller rational annihilator would remain a smaller annihilator in the splitting field, which is impossible. Rational initial-state basis solutions span the solution space, so the argument concerns observables of actual rational solutions as well, not only a formal eigenbasis.

This establishes sharpness for the entire observable space. It does not say every specific P, initial triple, or residual needs 2d+1 terms. Indeed the Ramanujan residual below is identically zero, and the cubic difference is just (−1)^n, of minimal order 1 in characteristic zero.

One useful structural interpretation over the splitting field: the three eigenmode sequences X=λ^n, Y=λ^(−n), Z=(−1)^n obey XY=Z². Thus homogeneous degree-d monomials reduce modulo this quadratic relation to a 2d+1-dimensional space. The independent μ_r establish that there are no further homogeneous relations in degree d for |t|>2. This explains why merely using the symmetric-power dimension binomial(d+2,2) still overcounts.

## 4. Independent Ramanujan application at t=83, d=3

The inherited generating functions have denominator 1−82z−82z²+z³. Their recurrence is C_83. Directly computed Lucas traces are L_1=83, L_2=6887, L_3=571538. Therefore

F_3(x)=(x+1)(x²−83x+1)(x²+6887x+1)(x²−571538x+1)

= x⁷−564733x⁶−3889880904x⁵+322812677460x⁴+322812677460x³−3889880904x²−564733x+1.

Each of a³,b³,c³ is annihilated by F_3. Also (−1)^n is annihilated, because x+1 divides F_3. Thus the residual a_n³+b_n³−c_n³−(−1)^n is annihilated. This does not require treating the subtracted exponential as an unjustified degree-zero term.

Independent exact numerator/denominator convolution gives:

| n | a_n | b_n | c_n | a_n³+b_n³−c_n³ |
|---|---:|---:|---:|---:|
|0|1|2|2|1|
|1|135|138|172|−1|
|2|11161|11468|14258|1|
|3|926271|951690|1183258|−1|
|4|76869289|78978818|98196140|1|
|5|6379224759|6554290188|8149096378|−1|
|6|529398785665|543927106802|676276803218|1|

The seven residuals are exactly zero. Since F_3 is monic with constant term 1, they prove the identity for all n∈Z. The backward initial triple is (−9,12,10), giving (−9)³+12³−10³=−1, hence 9³+10³=12³+1=1729.

## 5. Reproducible exact controls

Run `python verify_reciprocal_powers.py` in this directory; output is preserved in `verification.txt`. All 229 assertions passed. They comprise:

- Polynomial identities on symmetric-coordinate matrices over Z[T], d=0,…,3.
- Integer symmetric-coordinate matrices at t∈{−3,−2,−1,0,1,2,3,83}, d=0,…,6.
- Modular tests at t∈{−2,0,2,3,83}, d=1,…,4, modulo 2,4,6,8,9,25. Composite rings are intentional controls, not just field tests.
- Full tensor matrices at t∈{−2,2,3,83}, d=1,…,3.
- Independent seven initial residuals, coefficient factorization, backward extension and 1729 check.
- Nonzero cube sequences individually tested against the seventh-order recurrence, avoiding a vacuous test solely on an already zero residual.
- Hostile controls: perturbing a numerator breaks the identity; deleting x+1 from the proposed universal annihilator fails on the generic t=83 symmetric-coordinate matrix.

The proof establishes all parameters and rings. The computations audit implementations and special cases; finite tests alone are not the universal proof.

## 6. Prior-art and scope ledger

1. **Hirschhorn 1996, directly inspected primary author PDF:** M. D. Hirschhorn, *A proof in the spirit of Zeilberger of an amazing identity of Ramanujan*, Mathematics Magazine 69 (1996), 267–269. https://web.maths.unsw.edu.au/~mikeh/webpapers/paper40.pdf . This paper explicitly gives the seven spectral modes and seven-case sufficiency proof for the very same generating functions, then checks n=0,…,6. Thus the seven-check application is recovered published mathematics, not new session research.

2. **Hirschhorn 1995, directly inspected primary author PDF:** *An amazing identity of Ramanujan*, Mathematics Magazine 68 (1995), 199–201. https://web.maths.unsw.edu.au/~mikeh/webpapers/paper39.pdf . Gives explicit eigenmode expressions and a quadratic-form construction from the second-order sequence h_(n+2)=9h_(n+1)+h_n.

3. **Han–Hirschhorn 2006, directly inspected primary author PDF:** *Another look at an amazing identity of Ramanujan*, Mathematics Magazine 79 (2006), 302–304. https://web.maths.unsw.edu.au/~mikeh/webpapers/paper128.pdf . Gives the quadratic-form explanation and a 3×3 matrix formula, and explicitly includes negative indices. The square-eigenvalue relation explains 83=9²+2. These are established antecedents, not new constructions here.

4. **General recurrence-product literature:** Alecci–Barbero–Murru, *Some notes on the algebraic structure of linear recurrent sequences*, Ricerche di Matematica 74 (2025), 431–447, online 2023. https://link.springer.com/article/10.1007/s11587-023-00826-5 . The inspected introduction records the classical ring-level Hadamard-product closure and points to Cerruti–Vaccarino, *R-Algebras of linear recurrent sequences*, J. Algebra 175 (1995), 332–338. The older paper was not directly inspected here. This background rules out any claim that product closure itself is novel.

What this deliverable supplies: a self-contained uniform reciprocal-cubic formula at arbitrary degree; a polynomial-matrix proof valid under all commutative-ring specializations; an exact generic sharpness statement; reproducible controls; and corrected attribution of the concrete seven-check case. The bounded literature inspection does not establish whether the identical all-degree formula has appeared before, so no priority claim is made.

The traces L_r can be called reciprocal Lucas traces (or characterized simply by the displayed second-order recursion). This is an elementary recurrence/eigenvalue statement; neither the label Lucas nor proximity to Ramanujan supplies a modular form, an automorphic transformation law, a tournament invariant map, or a K3-period identification.
