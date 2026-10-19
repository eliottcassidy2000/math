# Elliptic interpolation of the Pythagorean AP seed, with an infinite positive Berggren family

**2026-10-19. PROVED structured construction.** The main content is an exact elliptic recurrence through the specified AP triangle, a non-torsion theorem for primitive Pythagorean seeds, and a positive family tied to an outer Euclidean/Berggren transformation. It is not a new label-size existence theorem for unrestricted sparse square-sum paths: elementary greedy extension already gives stronger numerical bounds there. Consecutive-interval coverage and preservation of all common-center AP petals remain separate constraints.

**Literature status.** Ordinary consecutive-interval square-sum existence is already established: Robert Gerbicz records paths for every n≥25 and cycles for every n≥32 in his January21,2018 [OEIS A090461 contribution](https://oeis.org/A090461), linking his deterministic construction. Any interval-coverage obstacle mentioned here concerns extending this particular structured construction or preserving additional prescribed data, not an open claim of ordinary existence.

## 1. The requested seed has a positive non-torsion elliptic continuation

On the integral curve

    E: y²=x³-482x²+6195836x-299007352,
    P=(482,51840),                                      (1)

put a_n=x(nP), n≥1. The sequence begins

    482, 3362, 6242, 158, 57442, 8198402/7921, ... .       (2)

The successive square-root witnesses begin

    62, -98, 80, 240, -21522/89, ... .                   (3)

Thus the first five terms already give an unscaled integer square-sum path. The extra original edge remains: 482+6242=82².

**Theorem 1.** All a_n are defined, positive, and pairwise distinct rational numbers. Every a_n+a_(n+1) is a rational square. Every finite prefix can be multiplied by one integer square to give a positive integer square-sum path containing a square-scaled copy of the original triangle. The required label bit length is O(k³ log(k+1)) for a k-term prefix; standard fixed-curve canonical-height theory sharpens this to O(k³), but the conservative bound already follows from the repository's elementary division-polynomial estimate.

Scaling is essential to this statement. At k=6 the least square-denominator scale is H=89, so the integer prefix is H² times (2); it no longer starts with the exact unscaled numbers. Elementary greedy methods can instead preserve those numbers exactly, but do not preserve this elliptic recurrence.

### Why consecutive sums are squares

For a generalized short equation y²=x³+a_2x²+a_4x+a_6, the chord/tangent slope s between Q and P satisfies

    x(Q)+x(P)+x(Q+P)=s²-a_2.

Here a_2=-482 and x(P)=482. They cancel, giving

    x(nP)+x((n+1)P)=s_n².                               (4)

Use the tangent at n=1. Later denominators never vanish: x(nP)=x(P) would imply nP=±P, impossible for n≥2 once P has infinite order.

### Positivity and distinctness

The derivative of the cubic in (1) is

    3x²-964x+6195836,

whose discriminant is -73420736<0. The cubic is strictly increasing, and its value at zero is -299007352<0. Its unique real root is therefore positive. Every finite real point on E has x at least that root, so every finite multiple has positive x.

The curve is nonsingular: a strictly increasing cubic has no repeated root here, or use the explicit discriminant -14863179074355855360000. Good reduction at7 and11 gives respectively10 and12 points, so rational torsion has order at most2. The point P has nonzero y and is not two-torsion, hence is non-torsion. The standard prime-to-residue-characteristic torsion injection suffices: the7-primary part injects at11 and vanishes, the11-primary part injects at7 and vanishes, and the remaining order divides gcd(10,12)=2.

Finally, x(nP)=x(mP) would imply nP=±mP. For positive distinct n,m this contradicts non-torsion. This proves positivity and distinctness independently of finite orbit checks.

## 2. A general three-term interpolation lemma

Let A,B,C be pairwise distinct rational numbers with

    A+B=s_1², B+C=s_2²,

where the signs of the rational square roots may be chosen. Put

    Y=-(s_1+s_2)(B-A)/2,
    U=2Ys_1-A², V=Y²-AU,
    E: y²=x³-Ax²+Ux+V, P=(A,Y).                         (5)

Assume the cubic is nonsingular and Y≠0. Then

    x(P)=A, x(2P)=B, x(3P)=C,
    x(nP)+x((n+1)P) is a rational square                 (6)

whenever the points are defined.

**Proof.** Equation (5) puts P on E. Its tangent slope is (A²+U)/(2Y)=s_1, giving x(2P)=s_1²-A=B. The chord slope from P to2P is -s_1-2Y/(B-A)=s_2, so x(3P)=s_2²-B=C. The general cancellation (4) proves the last assertion. ∎

Neither nonsingularity nor global positivity follows for an arbitrary input triple. For example A=14,B=11,C=110 with s_1=5,s_2=11 gives (x-2)²(x-10), a singular cubic. The Pythagorean specialization below has a separate nonsingularity proof.

For the requested triple, choose s_1=62,s_2=-98. Then Y=51840 and (5) gives exactly (1). The opposite relative sign of the two square roots is meaningful: the other choice gives a different interpolating curve and does not supply the same positivity argument.

## 3. Every primitive Pythagorean AP triangle supplies a non-torsion section

Let (a,b,c) be a primitive Pythagorean triple with positive legs, arranged with a<b, and suppose c²>4ab so the AP triangle is positive. Define

    A=2(c²-4ab), B=2c², C=2(c²+4ab),
    s_1=2(b-a), s_2=-2(a+b),
    Y=16a²b, U=4Y(b-a)-A², V=Y²-AU.                     (7)

The three pair sums are squares: A+B=s_1², B+C=s_2², A+C=(2c)². Formula (7) is exactly the interpolation (5).

**Theorem 2.** The curve E_(a,b):y²=x³-Ax²+Ux+V is nonsingular and P=(A,Y) has infinite order for every primitive Pythagorean seed. Global positivity of all later x-coordinates is an additional condition, proved for the family in Section4.

### Nonsingularity

The discriminant of the monic cubic factors exactly as

    -65536 a⁵b⁴(a³+6a²b+14ab²-2b³).                     (8)

Its last factor could vanish only if a/b were a rational root of t³+6t²+14t-2. The rational-root candidates ±1,±2 all fail. Thus (8) is nonzero. The elliptic discriminant is16 times this cubic discriminant.

### A sixth-multiple denominator proves non-torsion

Exact chord-and-tangent formulas give

    x(6P)=2(a⁴+4a³b+5a²b²-12ab³+4b⁴)/(a+2b)²,
    y(6P)=16a²b²(a²+6ab+16b²)/(a+2b)³.                  (9)

Put d=a+2b. The numerator of x(6P) is congruent to64b⁴ modulo d. Since gcd(a,b)=1, every odd prime p dividing d therefore gives

    v_p(x(6P))=-2v_p(d).

For a primitive Pythagorean triple, the even leg is divisible by4, so v_2(d)≤1. Exactly one leg is divisible by3, so 3 does not divide d. If a is even, the numerator in (9) has 2-adic valuation3, clearing the denominator's entire 2-part. Thus the reduced denominator is exactly

    denominator(x(6P)) = (d/gcd(d,2))² > 1.

Multiplying by9 and subtracting the integer3A does not change this denominator. Consequently X=9x(6P)-3A is nonintegral.
The rational change

    X=9x-3A, Z=27y

puts E into the integral short model

    Z²=X³+(81U-27A²)X+(729V+243AU-54A³).

Thus6P has a nonintegral x-coordinate on an integral short model. By the classical Nagell–Lutz integrality theorem it cannot be torsion, and neither can P. This uses an established theorem, not a torsion census. A primary teaching reference is Andrew Sutherland's [MIT18.782 Lecture24, Section24.3](https://ocw.mit.edu/courses/18-782-introduction-to-arithmetic-geometry-fall-2013/f37bbbd7fb8c7ee4e077cc5aa1d2a220_MIT18_782F13_lec24.pdf).

### The first two extra integer labels and a golden positivity threshold

The first five x-coordinates simplify to

    A, B, C,
    D=2(b²-a²-4ab),
    F=2(a²+17b²+4ab).                                  (10)

The third and fourth slopes are2b and6b. Therefore C+D=(2b)² and D+F=(6b)². They are distinct: D<A<B<C<F. The new label D is positive exactly when

    a/b < sqrt(5)-2 = phi^(-3).

This is an exact quadratic positivity condition, not an identification of the Gold and square-sum problems. The original AP condition is weaker, a/b<2-sqrt(3). Thus a positive AP triangle need not admit this particular two-integer extension with positive D, even though unconstrained greedy extension is always possible.

## 4. An explicit positive family on an outer Pythagorean parameter ray

For an even integer T≥0, set

    m=6T+5, n=5T+4,
    a=m²-n²=11T²+20T+9,
    b=2mn=60T²+98T+40,
    c=m²+n²=61T²+100T+41.                              (11)

The identity6n-5m=-1 proves gcd(m,n)=1. The parameters have opposite parity, so (a,b,c) is primitive. At T=0 it is (9,40,41), giving exactly482,3362,6242.

This is explicitly tied to the outer Euclidean/Berggren map

    L=[[2,-1],[1,0]], L⁴=[[5,-4],[4,-3]],
    (m,n)^T=L⁴(2T+1,T)^T.                              (12)

The T=0 source vector(1,0) is a boundary pair; its image(5,4) is the ordinary primitive seed. No claim that the boundary pair itself is a nondegenerate primitive triple is needed.

For the coefficients A,U,V in (7), all coefficients of the following polynomials in T are strictly positive. Ascending coefficient lists are

    A:
    482,2944,6484,6176,2162;

    3U-A²:
    18355184,170125952,686630080,1575122048,2244442016,
    2032158592,1140210880,361836928,49605104;

    -V:
    299007352,8974023936,88616121872,465807661952,
    1540102485064,3461402023936,5494932857312,
    6256099417344,5094172876424,2902177296128,
    1100255644432,249542740864,25624171192.

These exact identities prove A>0, U>A²/3 and V<0 for every T≥0. The derivative3x²-2Ax+U is everywhere positive. The cubic has one positive real root, and every finite real point has positive x. Theorem2 gives non-torsion. Hence every even T supplies an infinite distinct positive rational square-sum path through its AP triangle.

For the subfamily T=154k, there is also an independent non-torsion proof using only the fixed reductions7 and11: all curve coefficients and marked-point coordinates reduce to the T=0 data, with point counts10 and12. This confirms the original congruence-family construction; Theorem2 removes the need for that congruence restriction.

## 5. Uniform exact scaling and bit complexity

For an integral monic equation y²=x³+a_2x²+a_4x+a_6, a rational point's reduced x-denominator is a square. Indeed, at a denominator prime with v_p(x)=-e<0, the leading term x³ has strictly lower valuation than the other cubic terms. Thus2v_p(y)=-3e, forcing e even.

Write x(nP)=u_n/q_n² in lowest terms and put

    H_k=lcm(q_1,...,q_k),
    L_n=H_k² x(nP), 1≤n≤k.                             (13)

The L_n are distinct positive integers. Their adjacent sums are H_k²s_n². A rational number whose square is an integer is itself an integer, so these sums are actual integer squares. The first three labels are H_k² times the prescribed AP triangle.

For family parameter T, let L=1+log(T+2). The integral short model has coefficient and point weighted heights O(L). The universal division-polynomial coefficient bound proved in [the output-height packet, Section4](square_sum_contact_output_height_20261019.md) therefore gives

    h(x(nP))+h(y(nP))=O(n²[L+log(n+1)]).

This is a uniform argument from fixed-degree coefficient formulas and division-polynomial norms, not an extrapolation from finite tests. Summing denominator heights yields

    log H_k=O(k³[L+log(k+1)]).

Consequently every label in (13) has O(k³[L+log(k+1)]) bits. Sequential exact additions, gcd/lcm arithmetic, and certificate extraction take polynomial bit time in k and L.

For the T=0 seed, q_1=...=q_5=1 and q_6=89. Thus the first five labels are unscaled integers, while the six-term scaled prefix uses H_6=89. No assertion of infinitely many integral multiples on the fixed curve is made.

## 6. What this structure adds, and what remains

The elementary sparse-path construction already appends integer labels with O(k²) numerical size and can preserve the exact seed. A closure construction likewise gives polynomial-size cycles. Those stronger unconstrained existence bounds are recorded separately; they are the correct baseline.

The elliptic construction adds different information:

- the exact requested AP triangle occurs as P,2P,3P on one curve;
- every primitive Pythagorean seed produces a non-torsion section;
- explicit multiples give a golden-unit positivity threshold for a two-label continuation;
- a concrete outer-parameter family keeps the entire real orbit positive;
- all subsequent contacts are controlled by one group law, with exact rational height bounds.

These are structural links among Pythagorean triples, square-sum graphs, and elliptic dynamics. They do not supply a Hamiltonian path on a consecutive interval, nor preserve the common-center AP pairing of every later label. The interval-packing problem remains the meaningful square-sum obstacle after the unconstrained baseline is recognized.

## 7. Reproduction

From the repository root run:

    python 05-knowledge/experiments/verify_seeded_positive_elliptic_paths_20261019.py
    python -O 05-knowledge/experiments/verify_seeded_positive_elliptic_paths_20261019.py
    python 05-knowledge/experiments/verify_seeded_positive_elliptic_paths_20261019.py --length 64

The default runs produce identical JSON with4703 checks. Besides the first32 seed multiples and sampled family members, the checker verifies the displayed polynomial identities exactly in Z[a,b] and checks372 ordered primitive Pythagorean triples. The default receipt is [seeded_positive_elliptic_paths_controls.json](seeded_positive_elliptic_paths_controls.json); a longer receipt records4991 checks. These controls support the displayed identities. The all-parameter proofs are the algebraic identities, positivity coefficients, non-torsion argument, and height estimate above.
