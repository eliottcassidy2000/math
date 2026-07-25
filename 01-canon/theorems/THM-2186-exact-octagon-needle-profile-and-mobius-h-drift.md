---
id: THM-2186
title: "Exact octagon needle profile and Mobius H-drift"
status: >
  PROVED + VERIFIED-EXACT. For
  V_c=(1,c,c+1,2c,2c+1,...,6c,6c+1), the maximum lonely
  margin is an exact residue-eight rational function for every c>=1.
  The previously proposed value 7c/(8(7c+1)) holds exactly when 8|c.
  The proof identifies V_c as the difference set of seven arithmetic
  circle points and one needle point, classifies the four regular-octagon
  equality states, computes their four polyhedral gap gauges, and reduces
  the upper bound to nine independently audited finite cases. Along the
  stable dyadic branch the projective height H=M/(1/8-M) is exactly 7c,
  so scale insertion is a Mobius map in M. A general toric-needle theorem
  gives weak-* safe-mass convergence, an effective eventual exit, and an
  eventual residue-wise rational law at isolated torsion maxima.
source: codex-2026-07-24-octagon-needle
depends_on: []
related:
  - THM-2083
  - THM-2171
  - THM-2174
  - THM-2184
script: 04-computation/lrc14_octagon_needle_profile_thm2186.py
output: 05-knowledge/results/lrc14_octagon_needle_profile_thm2186.out
script_sha256: efb5ebe7527196472a49f8636a3dd86c9674d589f90d800446511660c96183ba
output_sha256: 82747a29ee0f8ef5a7d93e7cbfc9b7934a62418831b3ab47ad007e3f69cb97e3
hash_basis: working-tree bytes (LF)
---

# THM-2186 -- exact octagon needles and projective drift

For a finite integer row `V`, put

```text
M(V)=max_(t in R/Z) min_(v in V) ||vt||.              (1)
```

For an integer `c>=1`, define

```text
V_c=(1,c,c+1,2c,2c+1,...,6c,6c+1).                  (2)
```

At `c=1` this is a multiset profile because the first two entries agree.
For every `c>=2` it is a row of thirteen distinct positive speeds.

The original motivation is the zero-layer decomposition

```text
V_c=R+cU,
R=(1,0,1,0,1,...,0,1),
U=(0,1,1,2,2,...,6,6).                              (3)
```

The binary THM-2171 ladder uses `c=2^h`. The exact law below explains both
why the proposed closed form was seen on that ladder and why it is false on
an arbitrary affine slope.

## 1. Exact residue-eight law

Let `r` be the residue of `c` modulo eight. Then

```text
            7c
r=0:   ------------
        8(7c+1)

r=1:       1/8

            3c-2
r=2:   ------------
        8(3c+1)

            7c-5
r=3:   ------------
        8(7c+1)

            7c-4
r=4:   ------------
        8(7c+1)

            7c-3
r=5:   ------------
        8(7c+1)

            7c-2
r=6:   ------------
        8(7c+1)

            7c-1
r=7:   ------------ .                                (4)
        8(7c+1)
```

Equivalently,

```text
M(V_c)=1/8-m_r(c),                                   (5)

m_0=1/[8(7c+1)],       m_1=0,
m_2=3/[8(3c+1)],       m_3=6/[8(7c+1)],
m_4=5/[8(7c+1)],       m_5=4/[8(7c+1)],
m_6=3/[8(7c+1)],       m_7=2/[8(7c+1)].             (6)
```

In particular,

```text
M(V_c)=7c/[8(7c+1)]       iff       8 divides c.     (7)
```

This is an all-positive-integer statement, not an eventual asymptotic.

## 2. The hidden object is an eight-point circular code

For `(x,y) in (R/Z)^2`, set

```text
P_i=-ix,                    0<=i<=6,
Y=y,
Q(x,y)={P_0,...,P_6,Y},                              (8)

F(x,y)=min(
  ||y||,
  min_(1<=k<=6)||kx||,
  min_(1<=k<=6)||kx+y||
).                                                    (9)
```

Every difference `P_j-P_i` is `(i-j)x`, and every difference `Y-P_i` is
`y+ix`. Therefore

```text
F(x,y)=the minimum circular distance between two points of Q(x,y),
M(V_c)=max_t F(ct,t).                                (10)
```

This is the missing geometric object behind (2): seven points in one
arithmetic orbit and one additional needle point.

Eight cyclic gaps sum to one, so

```text
F(x,y)<=1/8.                                         (11)
```

Equality holds only when all eight gaps are `1/8`. Since
`P_i=iP_1`, regular-octagon equality forces

```text
P_1=a/8,        a odd,
Y=7a/8.
```

Writing `s=-a mod 8`, this is exactly

```text
E={ (s/8,s/8):s in {1,3,5,7} }.                     (12)
```

Conversely every point in (12) gives a regular octagon. The needle
`x=cy` meets (12) exactly when

```text
(c-1)s=0 mod 8.
```

Because `s` is odd, this is equivalent to `c=1 mod 8`. This proves the
`r=1` line of (4), including its equality classification.

## 3. The four exact octagon gap gauges

Near an equality point write

```text
x=s/8+u,              y=s/8+v.                       (13)
```

As long as the equality cyclic order is retained, every gap is `1/8`
minus one of the following linear deficits:

```text
s=1, order P0,Y,P6,P5,P4,P3,P2,P1:
  -v, 6u+v, -u, -u, -u, -u, -u, -u;

s=3, order P0,P5,P2,Y,P4,P1,P6,P3:
  5u, -3u, -2u-v, 4u+v, -3u, 5u, -3u, -3u;

s=5, order P0,P3,P6,P1,P4,Y,P2,P5:
  3u, 3u, -5u, 3u, -4u-v, 2u+v, 3u, -5u;

s=7, order P0,P1,P2,P3,P4,P5,P6,Y:
  u, u, u, u, u, u, -6u-v, v.                      (14)
```

Thus define the four positive homogeneous polyhedral gauges

```text
L_1(u,v)=max(-v,6u+v,-u),
L_3(u,v)=max(5u,-3u,-2u-v,4u+v),
L_5(u,v)=max(3u,-5u,-4u-v,2u+v),
L_7(u,v)=max(u,-6u-v,v).                             (15)
```

Inside the corresponding chamber,

```text
F(s/8+u,s/8+v)=1/8-L_s(u,v).                         (16)
```

The needle constraint becomes

```text
u-cv in (c-1)s/8+Z.                                  (17)
```

These gauges, rather than a residue label alone, are the critical
archimedean sidecar.

## 4. Exact support minima

For a real lift `rho`, minimize `L_s(u,v)` subject to

```text
u-cv=rho.                                             (18)
```

The unit balls of the four gauges have the following vertices:

```text
L_1<=1:
  (-1,-1), (-1,7), (1/3,-1);

L_3<=1:
  (-1/3,-1/3), (-1/3,7/3), (1/5,-7/5), (1/5,1/5);

L_5<=1:
  the negatives of the L_3 vertices;

L_7<=1:
  (1,1), (1,-7), (-1/3,1).                           (19)
```

Taking the maximum and minimum of `u-cv` over these polygons and scaling
gives

```text
             3rho/(3c+1),          rho>=0,
h_1(rho)={
             -rho/(7c+1),          rho<=0;

             5rho/(7c+1),          rho>=0,
h_3(rho)={
             -3rho/(7c+1),         rho<=0;

             3rho/(7c+1),          rho>=0,
h_5(rho)={
             -5rho/(7c+1),         rho<=0;

             rho/(7c+1),           rho>=0,
h_7(rho)={
             -3rho/(3c+1),         rho<=0.            (20)
```

Here `h_s(rho)` is the exact constrained minimum in (18). For example,
`{L_7<=1}` has support extrema

```text
max(u-cv)=7c+1,
min(u-cv)=-(3c+1)/3,                                 (21)
```

which gives the final two lines of (20). The other three follow directly
from (19).

For a residue `r`, put

```text
a_s=((r-1)s mod 8),          0<=a_s<=7.              (22)
```

If `a_s!=0`, only the two lifts

```text
rho=a_s/8,              rho=-(8-a_s)/8               (23)
```

can minimize (20); every other lift has the same sign and larger magnitude.
If `a_s=0`, the minimum is zero. Substitution into (20), followed by the
minimum over `s=1,3,5,7`, gives exactly (6). One set of winning lifts and
unit-ball directions is

```text
r   s    rho       (u/m_r,v/m_r)
0   7     1/8          (1,-7)
1   1       0           (0,0)
2   7    -1/8        (-1/3,1)
3   1    -3/4          (-1,7)
4   3     1/8       (1/5,-7/5)
5   1    -1/2          (-1,7)
6   5     1/8       (1/3,-7/3)
7   7     1/4          (1,-7).                       (24)
```

For every row of (24), set

```text
x=s/8+u,             y=s/8+v,             t=y.       (25)
```

Equation (17) holds, so `x=ct mod 1`. Directly from (14), all eight gaps
are positive and their minimum is `1/8-m_r`. This proves every lower bound
in (4), including the small values of `c`.

## 5. Localization proves all but nine finite cases

It remains to show that no other cyclic order does better. Suppose

```text
F(ct,t)>1/8-m,                  m=m_r(c)>0.           (26)
```

Sort the eight distinct points of `Q(ct,t)`, beginning at `P_0=0`, as

```text
0=z_0<z_1<...<z_7<1,
e_j=z_j-j/8.                                           (27)
```

Every cyclic gap is greater than `1/8-m`. Summing the first `j` gaps and
then their complementary `8-j` gaps gives

```text
-jm<e_j<(8-j)m,
|e_j|<7m.                                             (28)
```

Let `a_i` be the octagon rank of `P_i`. Since

```text
P_(i+1)=P_i+P_1 mod 1,
```

the multiple

```text
(a_(i+1)-a_i-a_1)/8
```

has a representative of absolute value less than `21m`. Therefore, if

```text
21m<1/8,                                              (29)
```

that multiple is zero modulo one. It follows inductively that

```text
a_i=ia_1 mod 8.                                      (30)
```

The seven `P_i` are distinct, so `a_1` is odd. The point `Y` occupies the
missing rank `7a_1=-a_1`; hence the configuration lies in one of the four
orders in (14). Equations (16)--(20) then force

```text
L_s(u,v)>=m_r,
```

contradicting (26).

Condition (29) covers

```text
r=0: every c>=8,       r=1: equality already handled,
r=2: every c>=26,      r=3: every c>=19,
r=4: every c>=20,      r=5: every c>=13,
r=6: every c>=14,      r=7: every c>=7.              (31)
```

The complete exceptional set is therefore

```text
X={2,3,4,5,6,10,11,12,18}.                           (32)
```

Exact affine-cell evaluation gives

```text
c:        2     3     4      5      6      10    11    12    18
M(V_c):  1/14  1/11  3/29   1/9   5/43   7/62  3/26  2/17  13/110.
                                                               (33)
```

These are precisely the values in (4). This completes the all-`c` proof.
QED.

As a quick independent obstruction to the discarded all-slope formula,
any positive maximum of a one-dimensional lower envelope of the functions
`||vt||` occurs at a kink or at the crossing of two affine branches. Its
reduced denominator is at most

```text
2 max(V_c)=12c+2.                                    (34)
```

The reduced denominator of `7c/[8(7c+1)]` is

```text
8(7c+1)/gcd(c,8).                                    (35)
```

When `8` does not divide `c`, (35) is at least `14c+2`, already larger
than (34). Thus the false extension cannot even be a lower-envelope
vertex value.

## 6. Independent exact verification

The companion has two separate exact evaluators.

1. The affine-cell evaluator partitions `[0,1]` at every kink
   `k/(2v)`. On each cell all `||vt||` are affine, so the lower envelope is
   concave and its maximum occurs at a cell endpoint or an in-cell
   pair-line intersection. It checks every case in (32).
2. The integer breakpoint evaluator independently enumerates denominators
   `2v`, `v+w`, and `|v-w|`, evaluates distances using integer residues,
   and compares rational values by cross multiplication. It checks every
   `1<=c<=512` and hostile slopes through `10002`.

The checker also evaluates the witnesses in (24), the projective identities
below, and the normal/optimized-Python transcripts. No floating-point
comparison enters the audit.

## 7. Toric-needle dichotomy

The octagon profile illustrates a general affine-ray theorem. Let

```text
R,U in Z^n,
F_(R,U)(x,y)=min_i ||U_i x+R_i y||,
M_c=max_t F_(R,U)(ct,t),                              (36)

Omega_alpha={(x,y):F_(R,U)(x,y)>alpha},
0<alpha<1/2.                                         (37)
```

Let `nu_c` be the pushforward of Haar measure on `R/Z` by

```text
t |->(ct,t).                                         (38)
```

At a character `(m,n) in Z^2`,

```text
hat(nu_c)(m,n)
 =integral exp(2pi i(mc+n)t)dt
 ={1, mc+n=0;
   0, otherwise}.                                    (39)
```

For every fixed nonzero `(m,n)`, the exceptional equality in (39) occurs
for at most one positive integer `c`. Trigonometric-polynomial density
therefore gives

```text
nu_c -> Haar measure on (R/Z)^2 weak-*.               (40)
```

The boundary of `Omega_alpha` lies in the finite union

```text
||U_i x+R_i y||=alpha.                               (41)
```

Every nondegenerate member of (41) is a finite union of affine torus
circles and has Haar measure zero. If some `(U_i,R_i)=(0,0)`, then
`Omega_alpha` is empty. Hence in all cases

```text
measure{t:min_i||(R_i+cU_i)t||>alpha}
   -> area(Omega_alpha).                              (42)
```

Because a nonempty torus-open set has positive area, there is a sharp
strict dichotomy:

```text
Omega_alpha empty:
  no affine-ray member is strictly alpha-safe;

Omega_alpha nonempty:
  every sufficiently large c is strictly alpha-safe. (43)
```

There is also an elementary effective form. If

```text
F_(R,U)(x_0,y_0)>=alpha+eta,
R_max=max_i |R_i|,
```

choose an integer `k` nearest to `cy_0-x_0` and put

```text
t=(x_0+k)/c.
```

Then `ct=x_0 mod 1`, `||t-y_0||<=1/(2c)`, and the circle norm is
one-Lipschitz. Therefore

```text
F_(R,U)(ct,t)>=alpha+eta-R_max/(2c).                 (44)
```

Every

```text
c>R_max/(2eta)                                       (45)
```

is strictly safe.

Writing `A=max F_(R,U)`, equations (43)--(45) separate three cases:

```text
A>alpha: eventual strict exit, effectively;
A<alpha: no member reaches alpha;
A=alpha: no strict exit, and weak safety is exactly incidence
         of the needle with {F_(R,U)=alpha}.          (46)
```

Thus an unbounded repeated zero-layer family `R+c_jU`, including
`c_j=q^j`, can remain strictly zero-safe only if its complete two-clock
profile is everywhere at most the target. At equality, the equality
skeleton is a necessary arithmetic sidecar.

## 8. Eventual torsion-tangent law

The critical case of (46) has a general exact refinement. Let `F` be the
minimum of finitely many torus-distance affine forms, let

```text
alpha=max F,
E={p:F(p)=alpha},                                    (47)
```

and assume:

1. `E` is a finite set of points on the `1/q` torus grid;
2. every `p in E` is an isolated polyhedral maximum.

Near `p`, there is then a positive homogeneous polyhedral gauge `L_p` such
that

```text
F(p+(u,v))=alpha-L_p(u,v).                           (48)
```

Let

```text
P_p={(u,v):L_p(u,v)<=1},
H_p^+(c)=max_(P_p)(u-cv),
H_p^-(c)=min_(P_p)(u-cv).                            (49)
```

For a fixed residue `c=r mod q`, put

```text
theta_(p,r)=r p_y-p_x mod 1,       0<=theta<1.       (50)
```

If `theta=0`, the needle passes through `p`. Otherwise the exact local
deficit is

```text
D_(p,r)(c)=min(
 theta/H_p^+(c),
 (1-theta)/(-H_p^-(c))
).                                                    (51)
```

Indeed, the local needle equation is

```text
u-cv in theta+Z.                                     (52)
```

Positive homogeneity and the two support extrema in (49) prove (51);
among positive lifts only `theta` can win, and among negative lifts only
`theta-1` can win.

Apply (43) at every level below `alpha`. It gives

```text
max_t F(ct,t) -> alpha.                              (53)
```

Hence every maximizing needle point eventually lies in one of the
neighborhoods (48), and

```text
alpha-max_t F(ct,t)=min_(p in E) D_(p,r)(c)           (54)
```

for all sufficiently large `c=r mod q`.

The polygon `P_p` has finitely many vertices. Consequently each support
extremum in (49) is the maximum or minimum of finitely many affine
functions of `c`, and one affine branch wins for all large `c`. Comparing
the finitely many expressions in (54) yields:

> **Torsion-tangent conclusion.** On each residue class modulo `q`, the
> critical needle deficit is eventually either zero or exactly
>
> ```text
> A/(Bc+D)                                           (55)
> ```
>
> for rational constants determined by one equality point, one sign of
> the phase lift, and one vertex of its tangent polygon.

The octagon calculation (14)--(33) is a case where an exact finite audit
removes the word "eventually".

## 9. Dyadic toothpick self-similarity and the H-coordinate

On the stable residue-zero branch,

```text
M_c=7c/[8(7c+1)],
d_c=1/8-M_c=1/[8(7c+1)].                            (56)
```

Define the projective height

```text
H_c=M_c/d_c=8M_c/(1-8M_c).                           (57)
```

Then

```text
H_c=7c.                                               (58)
```

Thus scale insertion `c |-> qc` is linear in `H` and Mobius in `M`:

```text
H_(qc)=qH_c,

M_(qc)=qM_c/[1+8(q-1)M_c].                           (59)
```

Equivalently, for `z_c=8d_c=1/(7c+1)`,

```text
z_(qc)=z_c/[q-(q-1)z_c].                             (60)
```

For the binary ladder,

```text
M(V_2)=1/14,
M(V_4)=3/29,
M(V_(2^h))=7*2^h/[8(7*2^h+1)]       for h>=3,       (61)

M_(2c)=2M_c/(1+8M_c)                 when 8|c.       (62)
```

This is the exact toothpick self-similarity. At the optimizer there are
seven equal short gaps and one long seam; binary insertion contracts the
seam excess by the projective law (60).

The projective coordinate is affine on every noncritical residue:

```text
r=0: H_c=7c,                 r=1: H_c=infinity,
r=2: H_c=c-2/3,             r=3: H_c=(7c-5)/6,
r=4: H_c=(7c-4)/5,          r=5: H_c=(7c-3)/4,
r=6: H_c=(7c-2)/3,          r=7: H_c=(7c-1)/2.      (63)
```

Thus the apparent nonlinear drift is projectivized affine scale.

For an odd radix `q`, its powers remain units modulo eight. If
`q=1 mod 8`, every power lies on exact octagon equality. If
`q=3,5,7 mod 8`, the powers alternate between `q mod 8` and `1`, hitting
equality every even exponent. Every even radix reaches residue zero by
the third power. This is the finite torsion clock coupled to the
archimedean drift.

## 10. Relation to THM-2184 and the LRC frontier

THM-2184 controls the **safe measure** of an endpoint-grid core plus an
affine tail `NLc_i+r_i`, with an effective `O(1/N)` error to a fixed
two-torus continuation profile. Equations (36)--(45) are its qualitative
empty-core/affine-ray counterpart for arbitrary threshold and arbitrary
tail length. The statements are complementary:

```text
THM-2184:
  quantitative safe-measure continuation from a grid core;

THM-2186:
  exact maximum profile, critical equality skeleton, and tangent drift. (64)
```

The row (2) can be placed under THM-2184 by taking core `{1}` and fixing
`c` modulo fourteen, but the present octagon argument is much sharper:

```text
M(V_2)=1/14,
M(V_c)>1/14                    for every c>=3.        (65)
```

Thus this complete zero-layer ray is discharged for LRC(14), with one
threshold equality at the consecutive row `{1,...,13}`. More generally,
(43) says that a repeated pump state with a nonempty strict two-clock body
cannot remain a counterexample at unbounded scale. What remains for a
general THM-2171 state is either a toric profile everywhere below target
or a critical equality skeleton requiring its own phase/tangent sidecar.

This theorem does not prove that every LRC(14) counterexample lies on one
affine ray and does not close the finite rank-eleven or rank-twelve
residuals.

## 11. Cyclic orders, tournaments, and an exploratory Fano probe

The intrinsic finite combinatorial datum in this proof is a cyclic order
of eight labelled phase points. A tournament orientation can name an
order chamber after choosing a cut, but it discards the gap lengths and
therefore does not preserve the target. The faithful typed carrier is

```text
vertices:             P_0,...,P_6,Y;
observable:           cyclic order and adjacent gap vector;
target:               minimum gap;
orientation gauge:    choice of circle cut and reversal;
lost by order alone:  all archimedean gap magnitudes;
needed sidecar:       one of the gauges L_s and its phase lift.          (66)
```

The `s=1,7` chambers are the ordinary and reversed octagon orders. The
`s=3,5` chambers are star-polygon orders. At the winning tangent vertices,
the `s=3` gap word has three short and five long gaps, while `s=5` has
five short and three long gaps.

> **EXPLORATORY ONLY -- OPEN.** The resulting `3+5` / `5+3` star-order
> split is a cheap object to compare with the still-open rank-one scalar
> `5+3` Fano carrier. No incidence-preserving map to that carrier
> has been proved. This observation is not a dependency, theorem,
> reduction, or claimed explanation of the Fano residual.

The proved lesson is narrower: at a critical toric maximum, retain the
equality skeleton, the cyclic chamber, and the polar gap gauge. A finite
phase label without its tangent magnitude repeats THM-2174's loss of the
archimedean scale coordinate. QED.
