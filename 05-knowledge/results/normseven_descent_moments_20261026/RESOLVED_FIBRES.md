# Two arithmetic gates: the same root torsor and good-prime counts, but lift degree two or four

2026-10-26. **PROVED**, R1-R4 below. This continues the incoming
`normseven_arithmetic_f1_20261026/boundary/RESULT.md` C2-C3, rather than
repeating its rational-node example. The new example has **nonrational
nodes** and compares their minimal resolution-lifting degrees. Its residue
torsor, exceptional conic and every good-prime local zeta function are explicit.

## R1. A simultaneous root-torsor and tangent-conic family

Over Q take nonzero D,delta and the same seven branch forms, with twists

    d=(2/delta, 1, 1, 1, -D, 1, -1),
    s_i^2=d_i L_i(E,x,z).

Use either marked norm-seven complementary pair. Above the rational triple
P=[0:0:1], the ordinary quotient product has exactly one geometric point y.
It is Q-rational: the triple 012 maps to a basis in each quotient deck,
so each quotient fibre has just one geometric point. This assertion is
about points, not the nonreduced fibre's scheme length.

The normal root surface X has eight geometric nodes above y:

    [s_0:...:s_6]=[0:0:0:1:epsilon sqrt(D):eta:theta],
    epsilon,eta,theta in {+1,-1}.                         (R1)

For nonsquare D these are four closed degree-two points, all with residue
field K=Q(sqrt(D)). The finite reduced fibre is four copies of Spec(K).
For square D it is eight rational points. In either case its entire
arithmetic structure is independent of delta.

Let pi: X_tilde -> X be the resolution obtained by blowing up these nodes.
The inverse image of each residue-field node is the conic

    Q_delta: A^2+B^2=delta C^2                           (R2)

over that node's residue field. To see the actual map, use s_3=1. Then
E+z=1, E=(delta/2)s_0^2, and the equations in the three vanishing roots give
s_1^2+s_2^2=delta s_0^2. The other three roots are units at (R1), with
nonzero derivatives, so adjoining them is etale locally. The projectivized
tangent cone is (R2), with [A:B:C]=[s_1:s_2:s_0]. Blowing up the nondegenerate
quadratic cone resolves it. This is the same local mechanism as incoming C2,
now over K instead of Q.

Thus for nonsquare D the **reduced complete resolved fibre over y** is four
copies of Q_delta considered as a Q-scheme through Spec(K)->Spec(Q).
This is restriction of the structure morphism, NOT Weil restriction.
For any extension field F/Q, this fibre has an F-point iff

    D is a square in F, and A^2+B^2=delta C^2
    has a nonzero solution over F.                       (R3)

Both directions retain y and the original twist. An F-point has to choose
an embedding K->F and then a point on its exceptional conic; conversely
these choices give a point upstairs. In particular adjoining the root-torsor
splitting field need not finish the resolution lift.

## R2. Same residue torsor, different minimum degree

Take D=13 and compare delta=2 with delta=3. The central divisor E=0,
the eight geometric nodes, their common field K and the original quotient
point y are the same. Only the transverse coefficient d_0 changes.

For delta=2, [A:B:C]=[1:1:1] is rational. Hence the minimum degree of a
number field F/Q admitting a point of the resolved fibre is exactly two:
K suffices and degree one cannot contain sqrt(13).

For delta=3 the conic has no K-point. Here is an elementary local proof
which uses no global local-to-global theorem. Since 13=1 mod3, the roots
of X^2-13 modulo3 are simple and lift digit by digit to Q_3. Explicitly,
from a root mod3^e choose a correction t3^e, t in {0,1,2}; its nonzero
derivative 2x determines the unique next digit. Thus K embeds in Q_3.
But A^2+B^2=3C^2 has no nonzero Q_3 solution: scale to integral coordinates
with at least one unit. Modulo3 forces 3|A,B; then the equation forces
3|C, a contradiction. A K-point would give a Q_3-point, impossible.

Any field supporting a resolved lift must contain K by R3. A degree-two
field containing K is K itself, which fails. The biquadratic extension

    F=Q(sqrt(13),sqrt(2)),  [A:B:C]=[1:sqrt(2):1]

works and has degree four by independent prime-valuation squareclasses.
All extension degrees containing K are even, so four is the minimum.
We have therefore proved

    min degree(delta=2,D=13)=2,
    min degree(delta=3,D=13)=4.                          (R4)

This is an actual obstruction to lifting a specified rational point through
normalization **and then resolution**. It is not a claim about all rational
points on either surface or about integral hourglasses.

The two gates can also help each other. With delta=3,D=2 the root-torsor
field already gives [1:sqrt(2):1], so the minimum is two. With D=5 use
[(1+sqrt(5))/2:(1-sqrt(5))/2:1]. Thus the extra resolution degree is not
always two: it depends on the conic over the actual residual field.
More generally any nonsquare D which is square in Q_3 gives the same
degree-four conclusion for delta=3; adjoining sqrt(2) always supplies the
upper bound, and its squareclass cannot equal that of D in this case.

## R3. Every good finite-field local zeta factor agrees

Now keep D=13 and reduce at any prime p outside {2,3,13}. For every
q=p^f and both delta=2,3, the complete resolved fibre has

    8(q+1) points if 13 is square in F_q,
    0 points otherwise.                                (R5)

Indeed the root-torsor gate supplies eight rational nodes or none. Every
nondegenerate conic A^2+B^2=delta C^2 over an odd finite field has q+1
points, independent of the nonzero delta. A direct proof is short:
count (x,y) with y^2=x^2-delta by (x-y)(x+y)=delta, obtaining q-1 pairs,
so sum_x chi(x^2-delta)=-1. The affine norm conic therefore has
q-chi(-1) points; its infinity has 1+chi(-1), giving q+1.

For a fixed base F_q and variable extension degree m, let T denote the
zeta variable. The local zeta function of this entire reduced fibre is

    [(1-T)(1-qT)]^-8                   if 13 is square in F_q,
    [(1-T^2)(1-q^2 T^2)]^-4            otherwise.        (R6)

This follows by substituting R5 into exp(sum_m N_m T^m/m). In the nonsplit
case odd counts vanish and every even extension has 8(q^m+1) points.
Consequently delta=2 and delta=3 have identical zeta factors for this
specified fibre at **every** good prime, including all extension counts,
although R4 distinguishes their minimum lift degrees over Q.

This does not contradict a local-global principle: the obstruction sits
at the excluded prime 3 after the K-embedding. “All good finite reductions”
omits p-adic solubility at bad places. Nor is it a counterexample asserting
equal global surface zeta functions; only this closed resolved fibre is
being counted. Geometric conics and their rational-coefficient cohomological
counts can forget arithmetic norm obstructions.

## R4. Relative F1 obstruction with the residue field retained

Over K the conic Q_2 is P^1 and has a split-torus decomposition. Q_3 has
no K-point by R2, and cannot be a finite disjoint union of nonempty split
K-tori: every such torus has its identity K-point. Thus their identical
good-prime counts do not certify a split-torus decomposition over K.

Over Q even the delta=2 complete fibre has no Q-point, because its nodes
have residue K. It too admits no decomposition into nonempty split Q-tori.
Do not move the relative torification statement from K to Q without that
residue-field sidecar. This is a local closed-fibre statement, not an F1
model for X, X_tilde or a general obstruction to all F1 frameworks.

## Contract, inheritance, and audit boundary

Source: the specified rational coarse point y and fixed twists (D,delta).
Target: four K-conics, their actual extension points, and their finite-field
point counts. Map: normalize to the eight nodes, then projectivize the
tangent cone; finite reduction preserves both named parameters.
Preserved: common source, residual field, transverse norm equation and
minimum splitting-field degree. Lost by good-prime counts: the norm
obstruction over K at its places above 3. Necessary sidecar: the tangent
quadratic form over K, not just its dimension, motive or smooth point count.
Decisive test: D=13,delta=2/3; positive interaction controls D=2/5,delta=3.

The incoming C2/C3 packet proves the rational-node mechanism and Gaussian
norm criterion over Q; it did not claim this comparison over nonsplit node
fields. R2's obstruction over K is proved directly and does not apply the
Q prime-factor criterion to an arbitrary number field. No Hasse principle,
Chebotarev density or Brauer-group classification is needed.
