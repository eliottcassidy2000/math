---
id: THM-4454
title: "Sharp global signed-root duplication stability"
status: >
  PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED. The exact
  best dimension-independent stability constant is K3; every eligible
  finite signed-root list has strict inequality and an explicit actual
  finite family approaches the constant. This does not prove actual
  Laurent two-rung separation or LRC(14).
source: complementary secant and regional tail packing, September 6, 2026
depends_on: []
primary_script: 04-computation/open_frontier_sep06_stability_complement.py
primary_output: 05-knowledge/results/open_frontier_sep06_stability_complement.out
primary_script_sha256: bc1c612f840d5e9f1e8135fa17086f3244442d13962e8e37b633080c77c2d5a1
primary_output_sha256: c43d2e8da2947154415299d155766777a3384754221d4526dff92b00ed3dda9d
report: 05-knowledge/results/open_frontier_sep06_stability_complement.md
audit: 05-knowledge/results/open_frontier_sep06_stability_global_audit.md
hash_basis: raw LF repository bytes
---

# THM-4454 — sharp global signed-root duplication stability

## 1. The sharp theorem

Let r be any finite real list, padded with zeros if necessary, with

    p1=sum r_i=1, p2=sum r_i^2=1, E=(1-p4)/2>0.

Put G(s)=product_i(1+r_i s) and

    D=-[s^4]G(s)^2=(5-8p3+3p4)/6, J=D/E,
    c_*=(13-8sqrt(2))/3.

Let a>=b>=0 be the largest two positive entries, using a zero if the
second is absent. The squared distance to the permutations of
(1/sqrt(2),1/sqrt(2),0,...) is

    d^2=2-sqrt(2)(a+b).

Then, for every such finite list,

    J-c_* > K3 d^2,
    K3=4sqrt(3)/[3(1+sqrt(2))(1+sqrt(3))].                (1)

This constant is best possible over all finite list lengths. There is
no finite equality case, but an explicit finite family below approaches
the quotient K3. In particular the exact global infimum is approximately
0.3501345012.

The formulas for E and D follow from Newton identities: e2=0,
e3=(p3-1)/3 and e4=(4p3-1-3p4)/12. The coefficient of s^4 in G^2 is
2e4+2e3. The distance formula follows by choosing the two largest
positive entries for the two nonzero target coordinates. It is positive:
equality in a+b<=sqrt(2) forces only two entries 1/sqrt(2), contradicting
p1=1. Arbitrary negative entries and zero padding are allowed.

## 2. One objective, two complementary regions

Write

    u=sqrt(2), v=sqrt(3), z=1/v, h=1/u,
    A=2-u, B=u-1, gamma=3K3/(4u)=(2-u)(3-v)/4,
    C0=A+u gamma=(1+u+v-uv)/2,
    C=C0-gamma(a+b)=A+(3K3/8)d^2.

Direct substitution gives the exact objective identity

    F_actual=1-C-p3+C p4
            =(3E/4)(J-c_*-K3 d^2).                     (2)

We prove F_actual>0, using a signed secant when b<=z and a convex
tail-square comparison when b>=z. Neither surrogate is asserted to
preserve p1; that normalization remains attached to the actual list.

## 3. The signed secant for b<=z

We have C>=A>0 and C<=C0-2gamma b. The function
b(C0-2gamma b) is increasing on [0,z]. Indeed
C0-4gamma z>0. Therefore

    1-2Cb >= 1-2z(C0-2gamma z)=2-v>0.

The function f_C(x)=x-Cx^2 is strictly increasing on [0,b]. Each positive
tail entry after a is bounded above by b, and a negative entry contributes
strictly less than zero to p3-Cp4. It follows that

    F_actual >= F(a,b),
    F(a,b)=1-C-a^3-(1-a^2)b+C[a^4+(1-a^2)b^2].           (3)

At fixed a,

    F_b=(1-a^2)T,
    T=gamma(1+a^2-2ab-3b^2)+2C0 b-1,
    U(b)=gamma(2-6b^2)+2C0 b-1,
    U-T=gamma[(1-a^2-b^2)+2b(a-b)]>=0.

On [0,z], U'>=2C0-12gamma z>0 and U(z)=2C0/v-1<0.
For exact sign certificates,

    2C0-12gamma z=7-2u-(5-2u)v>0,

since both compared sides are positive and their squared difference is
32u-42>0. Also U(z)<0 is equivalent to uv>1+u, whose successive
squared comparison reduces to9>8. The derivative bound above implies
C0-4gamma z>0 as used in the secant.

Thus F decreases strictly with b for a<1. It suffices to check
b=min(a,z,sqrt(1-a^2)). There are three cases.

For 0<=a<=z, b=a, the exact factorization is

    F(a,a)=2gamma(1-a)(a-z)(a-h)>=0.                     (4)

For z<=a<=u/v, b=z, put

    P(a)=a^3-(1+u)a^2+2a/3+2/3.

Then

    F(a,z)=gamma(1-a)(a-z)P(a)>=0,                      (5)

because P'<=8/3-2(1+u)/v<0 and
P(u/v)=2u(2v-3)/9>0.

Finally, for u/v<=a<1 and b=sqrt(1-a^2), put t=a+b and E2=a^2b^2.
Eliminating ab=(t^2-1)/2 gives

    g2=B-a^3-b^3+A(a^4+b^4),
    g2/[E2(2-ut)]=u(At+1)/(t+1)^2.

Here 1<t<=u/v+1/v<u. The right side implies

    4g2/[3E2(2-ut)]
       >=4u(3-u)/[3(1+u)^2]
        =4(13u-18)/3 >1/2>K3.                          (6)

The first strict comparison reduces to21632>21609. The second follows
from v/(1+v)<2/3 and u>1, which give K3<4/9. Thus F>0 on this last
boundary, except F(1,0)=0 at the excluded zero-energy endpoint.

The only envelope zeros in this region are (z,z) and (1,0). At (z,z),
a+b>1 forces a negative tail entry, making (3) strict. At (1,0), E=0.
Every eligible actual list with b<=z therefore has F_actual>0.

## 4. Convex tail packing for b>=z

Let c^2=1-a^2-b^2 be the total tail square mass. Then
0<=c<=z<=b, a+b>=2z and

    0<C<=C3=(3-v)/2<3v/8.

For 0<x<=c^2<=1/3, H_C(x)=x^(3/2)-Cx^2 satisfies

    H_C''(x)>=7v/4-3>0.

Strict convexity and H_C(0)=0 give
sum H_C(x_i)<=H_C(sum x_i). Replacing a negative tail entry by its
absolute value increases the cubic term strictly. Since a+b>1 forces
such a negative entry, we obtain the strict objective comparison

    p3-Cp4 < a^3-Ca^4+b^3-Cb^4+c^3-Cc^4.               (7)

It remains to prove nonnegativity for the three-root surrogate, which
has a>=b>=z>=c>=0 and square norm one. Write t=a+b. Its moments are

    p3=(3t-t^3-3tc^2)/2+c^3,
    p4=(3c^4-2c^2t^2-2c^2-t^4+2t^2+1)/2.

At fixed c, the exact allowed interval is

    z+sqrt(2/3-c^2)<=t<=sqrt(2(1-c^2)).                 (8)

For F(t,c)=1-C-p3+C p4, direct differentiation yields

    F_tt=3t-2C0(3t^2+c^2-1)+2gamma t(5t^2+3c^2-3).

This is strictly negative on (8). To check the entire interval, use
0<gamma<3/16, C0>27/32 and 8/7<2/v<=t<=u<10/7. The coefficient of
c^2 is negative, so F_tt<=Q(t), where

    Q(t)=10gamma t^3-6C0t^2+(3-6gamma)t+2C0,
    Q'(t)<(3/8)(15t^2-27t+5)<0.

The convex quadratic in this bound is negative at 8/7 and 10/7, with
values -307/49 and -145/49. Also

    Q(2/v)<35/(4v)-81/16<0,

whose last comparison is 140^2<81^2*3. The parameter bounds follow
from 24/17<u<17/12 and 19/11<v<26/15: gamma<35/187<3/16 and
C0=1-(u-1)(v-1)/2>61/72>27/32. These are exact rational-square checks.

Concavity reduces the minimum to the two endpoints of (8). At b=z,
the monotonicity of x-Cx^2 on [0,z] bounds the surrogate below by (5),
which is nonnegative. At a=b=x in [z,h], c=sqrt(1-2x^2), the exact
tail correction to the diagonal envelope is

    F(x,x,c)=2gamma(1-x)(x-z)(x-h)
              +c^2(x-c)[1-C(x+c)].

Using x-c=3(x-z)(x+z)/(x+c) and h-x=c^2/[2(h+x)], this becomes

    F(x,x,c)=c^2(x-z){
       3(x+z)/(x+c)[1-C(x+c)]-gamma(1-x)/(h+x)}.          (9)

The braces exceed 1/2. Indeed x+c<=2z, C<=C3,
1-C(x+c)>=2-v>1/4, (x+z)/(x+c)>=1,
(1-x)/(h+x)<=1 and gamma<1/4. Thus (9) is nonnegative, including
both endpoints. This and concavity prove surrogate nonnegativity.
The strict signed comparison (7) proves F_actual>0 in this region.

The two regions exhaust the domain and prove (1).

## 5. An actual finite family proves global sharpness

For each integer n>=3, set

    x_n=1+1/n, y_n=1-2/n, L=n^2, Q_n=3+6/n^2,
    q_n=(9-Q_n)/[L(3+sqrt(Q_n+(9-Q_n)/L))].

Take raw roots x_n,x_n,y_n and exactly L roots -q_n. Divide all roots
by S_n=3-Lq_n>0. The exact identity

    9-Q_n-6Lq_n+L(L-1)q_n^2=0

shows that p1=p2=1 after normalization. All roots are nonzero and E>0.
Since 0<q_n<2/L,

    S_n^2=Q_n+Lq_n^2<3+10/n^2<=3(1+1/n)^2.

Thus the two largest normalized positive roots exceed z for every n>=3.
As n tends to infinity, S_n tends to v, all three positive roots tend
to z, and the negative dust square mass tends to zero. The limiting
energy is 1/3 and the limiting squared distance is 2-2u/v>0. Consequently

    (J-c_*)/d^2 -> (3-4v/3-c_*)/(2-2u/v)=K3.

The normalization retains the nonvanishing negative dust first moment.
This is an actual finite-list sharpness sequence, not an infeasible
positive-only equality row or a fractional-multiplicity relaxation.

## 6. Scope, verification and lineage

The theorem closes the best dimension-independent constant for this
specified signed-root duplication quotient. It neither classifies all
near-minimizing sequences nor supplies the actual binomial coefficient
transport needed for general Laurent two-rung separation. No external
priority or proof-assistant claim is made.

It improves the earlier explicit constants in
[the first quantitative bound](../../05-knowledge/results/quantitative_stability_empty_core_morning_sep06.md)
and [the energy-adjusted bound](../../05-knowledge/results/open_frontier_sep06_stability.md).
Those inequalities remain valid; their open optimal-constant question is
settled by this theorem. The proof above is elementary and self-contained.

The [complement source](../../04-computation/open_frontier_sep06_stability_complement.py)
has 90 exact gates: formal identities over Q(sqrt2,sqrt3), interval-wide
radical signs, and 30 eligible actual rows from 49 declared rational prefixes.
The [packing source](../../04-computation/open_frontier_sep06_stability_packing.py)
has 72 exact gates: full curvature/boundary identities, 24 declared prefix
controls, four actual sharpness-family rows and surrogate boundary controls.
All exclusions are printed. The banks check implementations; the displayed
analytic inequalities and limits prove the universal theorem.

Root independently audited the regional proof/source and replayed its
normal/optimized output. The regional author independently audited the
complementary proof/source and its normal/optimized output. A third reader
audited the entire global argument, strictness and actual sharpness.
See the [complementary audit](../../05-knowledge/results/open_frontier_sep06_stability_complement_audit.md),
[global audit](../../05-knowledge/results/open_frontier_sep06_stability_global_audit.md)
and [root audit](../../05-knowledge/results/open_frontier_sep06_root_audit.md).
Reproduction commands and all frozen hashes are in those reports.
