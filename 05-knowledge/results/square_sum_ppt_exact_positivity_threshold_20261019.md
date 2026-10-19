# Exact positivity threshold for primitive-Pythagorean elliptic seeds

**2026-10-19. PROVED, conditional only on the established primitive-PPT non-torsion theorem and standard real elliptic group topology.** This strengthens the explicit positive parameter family in [the seed interpolation paper](square_sum_ppt_seed_elliptic_interpolation_20261019.md) to an exact classification of every positive primitive-Pythagorean AP seed. It concerns this specified elliptic recurrence, not ordinary square-sum existence, which Gerbicz settled in2018 for all sufficiently large consecutive intervals.

## 1. Statement

Let (a,b,c) be a primitive Pythagorean triple with 0<a<b and c²>4ab. Use the curve and point

    A=2(a²+b²−4ab), Y=16a²b,
    U=4Y(b−a)−A², V=Y²−AU,
    E: y²=x³−Ax²+Ux+V, P=(A,Y).

The inherited theorem proves E smooth and P non-torsion. Its abscissas begin with the AP triangle A,2c²,2(c²+4ab), and every consecutive sum is a square.

Define

    R(t)=t⁶+4t⁵+3t⁴−8t³+35t²−12t+1.

Let beta be the larger of its two roots in (0,2−sqrt(3)). Numerically beta≈0.228228400091583; this decimal is descriptive, not its definition.

**Theorem.** Every x(nP), n≥1, is strictly positive if and only if

    a/b < beta.

If beta<a/b<2−sqrt(3), infinitely many positive multiples nP have negative abscissa. All multiples remain defined and have pairwise distinct abscissas in both regimes.

The threshold is specific to the chosen relative signs of the two initial square-root witnesses. It does not rule out a different recurrence or an integer extension of the same seed.

## 2. Normalization and the two real regimes

Scale x=b²X,y=b³Z and set t=a/b. The normalized coefficients are

    A₀=2(t²−4t+1), Y₀=16t²,
    U₀=64t²(1−t)−A₀², V₀=8R(t).

The cubic discriminant has sign opposite to

    G(t)=t³+6t²+14t−2.

Since G'(t)=3t²+12t+14>0 for t>0, G has one positive root rho. It lies between13488/100000 and13489/100000. Neither rho nor beta is rational: for G the rational-root candidates ±1,±2 fail, and for the monic integer polynomial R with constant1 the only rational-root candidates ±1 fail.

For t<rho the cubic has three distinct real roots; for t>rho it has only one. The primitive-PPT input has rational t, so the equality case cannot occur.

## 3. Why the entire three-root regime stays positive

Write F(X)=X³−A₀X²+U₀X+V₀. We have A₀>0 by the positive-AP condition, and

    F'(A₀)=A₀²+U₀=4Y₀(1−t)>0.

When F has three roots, its larger critical point is

    r_c=(A₀+sqrt(A₀²−3U₀))/3.

The inequality U₀>−A₀² gives r_c<A₀. Since F(A₀)=Y₀²>0, A₀ lies strictly above the largest root r₃. Also r₃≥A₀/3>0 because the sum of the three roots is A₀.

The point P therefore lies in the identity component of E(R), namely its unbounded real component. That component is a subgroup, so every multiple remains there, with x≥b²r₃>0. This does not require the whole real curve to have positive abscissas: its other component may contain negative points, but this orbit never enters it.

## 4. The sign threshold in the connected regime

For 0<t<alpha=2−sqrt(3),

    R''(t)=30t⁴+80t³+36t²−48t+70>0,

since alpha<1. Thus R is strictly convex on the AP interval. We have R(0)=1 and R(alpha)=32alpha⁴>0, the latter following from A₀=0 and V₀=Y₀² at the endpoint. Direct exact evaluations give R(1/5)=−904/15625<0. Hence R has precisely two zeros in this interval.

Both endpoints of the rational bracket for rho have negative R-values. Convexity therefore gives R(rho)<0. Consequently the smaller R-root lies below rho and the larger root beta lies above it. For rho<t<beta, V₀<0; for beta<t<alpha, V₀>0.

In the connected regime F has one real root r. Its sign is opposite the sign of V₀: F(0)=V₀ and the leading coefficient is positive. Thus r>0 for rho<t<beta, and every finite real curve point has positive x. Combined with Section3, this proves positivity throughout0<t<beta.

For beta<t<alpha, the real root is negative. The real curve contains a nonempty open set with r<X<0. A connected real elliptic curve is a circle group, and the positive multiples of a non-torsion point are dense in that circle. Therefore infinitely many nP enter this negative-abscissa open set. This proves necessity as well as sufficiency.

The standard real-group facts are stated explicitly in Noam Elkies's [Harvard Math223 elliptic-curve notes](https://people.math.harvard.edu/~elkies/M223.24/index.html): the real group has one or two circle components, and a non-torsion point generates a dense subgroup of its circle. Density here is used only to prove eventual sign failure, not to establish any integer coverage claim.

## 5. Golden five-term positivity is weaker than orbit positivity

The inherited fourth abscissa is

    x(4P)=2(b²−a²−4ab),

so the first five labels are positive exactly when t<sqrt(5)−2. The exact all-orbit threshold beta is smaller. For example beta<23/100<sqrt(5)−2 follows from R(23/100)>0 and the preceding root classification, with23/100>rho.

An explicit primitive triple comes from Euclidean parameters(m,n)=(34,27):

    (a,b,c)=(427,1836,1885), t=427/1836.

It lies between beta and sqrt(5)−2. Its first eight abscissas are positive, but

    x(9P)=−126100899082622355694/917204544598225.

Thus the golden condition ensuring two extra positive integer labels cannot be promoted to an infinite positive orbit. The counterexample is to that promotion, not to square-sum extendibility of the seed.

## 6. Consequences for the Pythagorean tree

The positive family is much broader than one parameter ray. In Euclidean coordinates a=m²−n²,b=2mn with m>n and a<b, the exact condition becomes

    1 < m/n < beta+sqrt(1+beta²).

Every coprime, opposite-parity pair in this interval gives an infinite positive rational square-sum recurrence through its AP seed. Such pairs are infinite; for instance m=n+1 yields t→0 through valid primitive pairs. The previously constructed outer Berggren ray remains one explicit subfamily.

No single finite integer scaling covers this infinite rational orbit: density in its identity circle gives infinitely many distinct abscissas inside a bounded real interval, whereas a fixed-denominator lattice contains only finitely many numbers there. Finite prefixes retain the inherited exact square-scaling and height bounds. Ordinary interval Hamiltonicity, prescribed seed-edge extension, and common-center AP-petal compatibility are separate properties.

## 7. Reproduction

Run `python 05-knowledge/experiments/verify_ppt_positivity_threshold_20261019.py`. It uses exact rational arithmetic for the discriminant bracket, polynomial sign witnesses, the counterexample's first nine multiples, and its PPT/interval conditions. The receipt is [ppt_positivity_threshold_controls.json](ppt_positivity_threshold_controls.json). The real-topology argument and convexity proof above establish the infinite conclusions.
