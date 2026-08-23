---
id: THM-3872
title: "Three-cusp polarization branches and the minimal affine square-residual gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In
  the THM-3864 seminormalization, simultaneous square/cube descent is exactly
  the product condition H(a)H'(a)=0 at each of the three cusp addresses.
  Hence the descent locus is the union of eight linear cusp-polarization
  branches.  Modulo the branch-ring ideal vanishing at all three cusps, every
  defect class has an explicit affine representative classification.  In the
  resulting seven noncanonical minimal slices, the only nontrivial square
  residual is the THM-3869 representative h_1-4(x+1), whose residual is
  (9x+4)^4.  Its extra Cardano line remains genuinely ramified.  The three
  natural generators of the omitted cusp-value-zero ideal produce no new
  square on their entire constant linear span.  Polynomial-coefficient and
  higher additions from that ideal, alternative Delta-lifts, a Keller atlas,
  and JC(2) remain OPEN.
source: jc_zero_debt_lift / post-THM-3864 noncanonical representative lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_quartic_c3_construct, 2026-08-23).
  The audit independently rederived the local square/cube iff, the selector
  and derivative matrices, the exact R/J fibres and all eight strata, every
  scaling boundary, and the necessity and sufficiency of all square
  recurrences.  It also reconstructed J as the reduced three-point ideal and
  checked the full constant span: the odd-degree splits force the sole seam,
  while an independent sextic recurrence has x-adic obstructions at u=+/-6
  and remainder gcd u^3 elsewhere.  This confirms the saturated unit ideal
  without relying on the primary Groebner path.  The companion verifies the cusp-value
  selector and defect-derivative matrices, every mixed product descent, exact
  Delta divisibility for all searched families, the unique h_1-4(x+1)
  square, two unit sextic-square ideals on the d_0 branch, the complete
  two-parameter square recurrence and unit ideal on the d_+ branch, both
  two-cusp odd-degree ladders, the final plus/minus exceptional quartic, and
  the normalization involution.  It also verifies the exact three-generator
  presentation of the cusp-value-zero ideal, all three mixed descents, and
  the three generator-ray obstructions, the full constant-span leading seam,
  and its saturated sextic-square unit ideal.  Normal and optimized runs
  byte-match the frozen transcript; the polynomial-coefficient scope remains
  explicitly open.
depends_on:
  - THM-3864-integrated-three-cusp-conductor-seminormal-three-direction-gate
  - THM-3869-three-cusp-square-residual-cardano-line-ramification
related:
  - THM-3854-integrated-three-cusp-quintic-s5-natural-completion-obstruction
  - THM-3874-three-cusp-quadratic-k3-affine-class-group
script: 04-computation/jc2_three_cusp_polarization_branches_thm3872.py
output: 05-knowledge/results/jc2_three_cusp_polarization_branches_thm3872.out
script_sha256: ccc01d00f16f021f808bff0b6201a305c22fd3a96a1f638494e68600e7877a53
output_sha256: 7a94c8b43527961378c7f03bab377c584df3262deafb428196eb56052636239c
semantic_sha256: f3a13ca8c79e4e782fee246a3231531ea22caadbaf2924048348e63beb4dae7a
hash_basis: raw LF bytes
---

# THM-3872 -- eight cusp polarizations, one minimal square

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero.  Retain
the THM-3854 normalization and the THM-3864 rings

```text
x=t^4-2t^2,                   y=3t^5-5t^3,
R=k[x(t),y(t)] subset S subset k[t],                         (1)
```

where `S` is the seminormalization of `R`.  The cusp addresses are

```text
A={0,1,-1}.                                                   (2)
```

Let `h_1,h_2,h_3` be the defect basis from THM-3864.  This theorem has two
parts.

1. The full square/cube-descent locus in `S` is the union of eight explicit
   linear cusp-polarization branches.
2. After quotienting branch-ring additions by the ideal which vanishes at
   all three cusps, the smallest affine section of all seven noncanonical
   branches has exactly one nontrivial square residual.  It is

```text
h_*=h_1-4(x+1),                 residual=(9x+4)^4.            (3)
```

The second assertion is deliberately about the displayed affine section and
the displayed polynomial descents.  It does not classify arbitrary additions
from the cusp-value-zero ideal or arbitrary changes of the two lifts.

## 1. The local product law is an exact iff

THM-3864 proves the global descriptions

```text
S={polynomials satisfying the three node-pair value equalities},
R={F in S : F'(a)=0 for every a in A}.                       (4)
```

Let `H in S`.  Its powers automatically retain every node equality.  At a
cusp `a`,

```text
(H^2)'(a)=2H(a)H'(a),
(H^3)'(a)=3H(a)^2H'(a).                                     (5)
```

Consequently

```text
H^2,H^3 in R  iff  H(a)H'(a)=0 for every a in A.             (6)
```

In fact `H^2 in R` already implies the right side, and the right side implies
both power conditions.  There is no quotient loss in `(6)`.

For each cusp independently, choose either the value equation `H(a)=0` or
the derivative equation `H'(a)=0`.  Each of the `2^3=8` choices cuts out a
linear subspace of `S`, and `(6)` says that their union is exactly the descent
locus.  Intersections record points for which both local coordinates vanish;
they do not create extra branches.

## 2. Exact affine fibers modulo the cusp-value-zero ideal

The three branch-ring selectors

```text
e_0=x+1,
e_+=(-2x-y)/4,
e_-=(-2x+y)/4                               in R             (7)
```

have value matrix

```text
             t=0   t=1   t=-1
e_0            1     0      0
e_+            0     1      0
e_-            0     0      1.                              (8)
```

Thus cusp evaluation `R -> k^3` is onto.  Put

```text
J={r in R : r(0)=r(1)=r(-1)=0}.                             (9)
```

Then `R/J=k e_0 direct_sum k e_+ direct_sum k e_-`.

For `C=(C_1,C_2,C_3)`, let `h_C=sum C_ih_i`.  It vanishes at all three
cusps, and its derivative vector is

```text
d_0=-20C_2,
d_+=10(-C_1-2C_2+2C_3),
d_-=10(C_1-2C_2+2C_3).                                    (10)
```

Every representative of its defect class is `h_C+r`, `r in R`.  Modulo `J`
it has the unique form

```text
H=h_C+v_0e_0+v_+e_++v_-e_-.                                (11)
```

Because every `r in R` has zero cusp derivatives, `(6)` becomes exactly

```text
v_0d_0=v_+d_+=v_-d_-=0.                                    (12)
```

This classifies the affine fibers: a cusp value is forced to zero precisely
when the corresponding defect derivative is nonzero, and it is free when
that derivative vanishes.  The eight derivative strata and their minimal
representative slices are therefore:

```text
no zero d_i:       h_C;
d_0=0:             C_2=0,                 h_C+a e_0;
d_+=0:             C_1=-2C_2+2C_3,       h_C+a e_+;
d_-=0:             C_1= 2C_2-2C_3,       h_C+a e_-;
d_0=d_+=0:         C proportional ( 2,0,1), h_C+a e_0+b e_+;
d_0=d_-=0:         C proportional (-2,0,1), h_C+a e_0+b e_-;
d_+=d_-=0:         C proportional ( 0,1,1), h_C+a e_++b e_-;
all d_i=0:         C=0,                    H in R.           (13)
```

The last equality uses the nonzero determinant `-8000` of the derivative
matrix.  It gives only the trivial descents `P=H^2,Q=H^3` and zero residual,
not a cubic extension.

## 3. The bounded polynomial-descents universe

For each `h_C`, use the canonical square/cube descents `P_C,Q_C` of
THM-3864.  For a selector addition `r` in `(11)`, the product `rh_C` belongs
to `R` by `(12)`; use the minimum-degree mixed descent `A`.  Then the searched
lifts are exactly

```text
P=P_C+2A+r^2,
Q=Q_C+3rP_C+3rA+r^3,                                      (14)
```

and pull back to `(h_C+r)^2,(h_C+r)^3`.  Exact division defines

```text
R_(C,r)=(P^3-Q^2)/Delta in k[x,y].                          (15)
```

The companion verifies every mixed descent and every division in `(14)-(15)`.
Since `k[x,y]` is a UFD and `k` is algebraically closed, a polynomial which is
a square in `k(x,y)` is already a polynomial square up to a square constant.
All nonzero specializations below are therefore valid necessary square tests.

## 4. The `d_0=0` branch has the unique positive point

Write

```text
H=C_1h_1+C_3h_3+a e_0.                                    (16)
```

The two mixed descents used in `(14)` are

```text
e_0h_1=-11x^2-7x+y^2,
e_0h_3=27x^2y+69xy+y^3+38y.                                (17)
```

First take `C_3=0` and scale `C_1=1`.  The `y=0` residual is the quartic

```text
G_a(x)=
 6561x^4+(3888-486a^2-3888a)x^3
 +(-3a^4+72a^3+234a^2-2352a)x^2
 +(-6a^4+24a^3+336a^2)x
 -3a^4-16a^3.                                               (18)
```

Comparison with `(A_2x^2+A_1x+A_0)^2` has reduced lexicographic basis

```text
256A_2-243A_0a^2-1944A_0a-5184A_0,
32A_1-21A_0a^2-168A_0a-480A_0,
A_0^2+96a^2+768a+1280,
(a+4)^3.                                                     (19)
```

Thus `a=-4`, and direct substitution gives

```text
R_(1,0,-4)=(9x+4)^4.                                       (20)
```

It remains to prove that adding an `h_3` component creates no further
square.  Scale `C_3=1` and write `c=C_1`.  If `c!=0`, the specialized
quartic is `c^6G_(a/c)(x)`, so `(19)` forces `a=-4c`.  On that seam,

```text
R_(c,1,-4c)=(c+y)^2K_c,

K_c(0,y)=-16(-16c^4-64c^3y+10488c^2y^2-127072cy^3
              +445835y^4+37044cy^5+37044y^6).              (21)
```

The coefficient ideal for `K_c(0,y)` to be a cubic square is the unit ideal
in `k[c]`.  If `c=0`, the original residual is `y^2L_a`, where

```text
L_a(0,y)=
 -3a^4-472a^3y-27816a^2y^2-(8a^3+727776a)y^3
 - (1008a^2+7133360)y^4-42336ay^5-592704y^6.               (22)
```

Its cubic-square coefficient ideal is again the unit ideal.  Equations
`(18)-(22)` prove that, up to overall scalar, `(3)` is the unique nontrivial
square residual on the whole minimal `d_0=0` slice.

## 5. The one-cusp `d_+=0,d_-=0` branches are empty

On the open `d_+=0` stratum one has `C_2!=0`; scale `C_2=1` and put

```text
(C_1,C_2,C_3)=(-2+2s,1,s),                 r=a e_+.          (23)
```

Let `p(x)=R_(s,a)(x,0)=sum_(i=0)^7 p_ix^i`.  Its constant term is
`p_0=-320000`.  If it were the square of a cubic, put `u_i=A_0A_i` for the
square-root coefficients.  The first three rows force

```text
u_1=p_1/2,
u_2=(p_2-u_1^2/p_0)/2,
u_3=(p_3-2u_1u_2/p_0)/2.                                   (24)
```

The remaining necessary and sufficient rows are

```text
p_7=0,
p_4=(u_2^2+2u_1u_3)/p_0,
p_5=2u_2u_3/p_0,
p_6=u_3^2/p_0.                                              (25)
```

After clearing only the nonzero constant denominator `p_0`, the exact ideal
of `(25)` in `k[s,a]` has Groebner basis `[1]`.  Hence this whole affine slice
is empty.  The involution `t -> -t`, `(x,y)->(x,-y)` swaps `e_+` with `e_-`
and sends `h_2` to `-h_2`; it proves the identical statement for `d_-=0`.

## 6. The two-cusp intersections are empty

For `d_0=d_+=0`, scale the defect class to `2h_1+h_3` and put

```text
H=2h_1+h_3+a e_0+b e_+.                                    (26)
```

The `y=0` residual has degree seven with leading coefficient

```text
[x^7]R=-6561b^3/8.                                          (27)
```

Thus `b!=0` is excluded by odd degree.  At `b=0`, `(26)` is the `c=2`
instance of Section 4 and is already nonsquare.  The normalization involution
gives the `d_0=d_-=0` branch.

For `d_+=d_-=0`, scale the class to `h_2+h_3` and write

```text
H=h_2+h_3+a e_++b e_-.                                     (28)
```

The first two leading rows are

```text
[x^7]R=-6561(a-b-152)^3/8,
[x^6]R=-729(a-b-152)^2(5a-5b-652)/2.                        (29)
```

On the seam `a-b=152`, the next odd row is

```text
[x^5]R=-(243*277248/64)(a-76)^2.                            (30)
```

Hence odd degree excludes every point except `(a,b)=(76,-76)`.  There

```text
R(x,0)=-512(9x+5)^4,                                       (31)
```

but at `x=0`, with `q=y^2`, one gets

```text
-592704q^4+4946089q^3-11376906q^2+210000q-320000.           (32)
```

The coefficient ideal for `(32)` to be a quadratic square is `[1]`.  This
closes the final nontrivial two-cusp branch.

## 7. Exact conclusion and the surviving payment

Sections 4-6 exhaust the seven noncanonical branches of the minimal affine
section `(11)-(14)`.  The only nontrivial square residual is `(20)`, namely
the representative

```text
h_*=(t^2-1)(9t^4-18t^2+4)=h_1-4(x+1).                      (33)
```

THM-3869 constructs its natural nonmonogenic cubic maximal order and proves

```text
Disc_div=V(Delta)+2V(9x+4).                                 (34)
```

Thus the unique minimal square does not remove the extra Cardano
ramification: it merely converts the order discriminant factor `(9x+4)^4`
to the genuine field-discriminant factor `(9x+4)^2`.

The proved and independently audited THM-3874 K3/class-group theorem gives the
global version of this warning: `Cl(Q)=Z` rules out every connected cyclic
cubic cover of the quadratic resolvent which is unramified in codimension one.
Thus **any** future representative from the omitted `J`-directions must pay
some extra divisor, not merely that the minimal affine section pays `9x+4`.
This consequence does not close the omitted mixed or polynomial
`J`-directions themselves, and it is not used in the proof above.

## 8. The three first `J`-generator rays do not deform the square

The ideal omitted by the finite quotient has the exact presentation

```text
J=(j_1,j_2,j_3)R,
j_1=x(x+1),              j_2=y(x+1),              j_3=y^2+4x.       (35)
```

Indeed these polynomials have Groebner basis

```text
4x+y^2,                    xy+y,                    x^2+x,           (36)
```

so their quotient has basis `1,x,y`, dimension three, and the evaluation
matrix `(8)` identifies it with the three reduced cusp points.  Thus `(35)`
is the full cusp-value-zero ideal, not merely an annihilating subideal.

The minimum mixed descents of `j_ih_*` are

```text
A_1=-15x^3-15x^2+xy^2-4x,
A_2=-15x^2y-15xy+y^3-4y,
A_3=81x^4+9x^3-44x^2+15xy^2-16x+4y^2.                     (37)
```

For each `i`, use the natural one-parameter lifts

```text
P_lambda=P_*+2lambda A_i+lambda^2j_i^2,
Q_lambda=Q_*+3lambda j_iP_*+3lambda^2j_iA_i+lambda^3j_i^3. (38)
```

They pull back to `(h_*+lambda j_i)^2,(h_*+lambda j_i)^3`.  Exact division
gives three immediate obstructions:

```text
i=1: deg_y R=2, [y]R=0, [y^2]R=-8lambda^3x^3, R(0,0)=256;
i=2: [y^5]R(0,y)=-8lambda^3;
i=3: [x^7]R(x,0)=52488lambda^3.                             (39)
```

For `i=2,3`, nonzero `lambda` gives an odd-degree nonzero specialization.
For `i=1`, a square root would have the form `U(x)+yV(x)`.  Its missing
linear `y` term forces `UV=0`; `U=0` contradicts the nonzero constant part,
while `V=0` forces `lambda=0` through the displayed `y^2` coefficient.
Thus every generator ray returns to the original square only at
`lambda=0`.

### 8.1 The entire constant generator span is empty

In fact the three rays interact without creating a hidden square.  Put

```text
r=uj_1+vj_2+wj_3,                  A=uA_1+vA_2+wA_3,        (40)
```

and use `(14)` with `h_*`, `r`, and `A`.  If `w=0`, the specialization at
`x=0` has degree at most five, and its degree-five coefficient is

```text
[y^5]R(0,y)=-8v^3.                                          (41)
```

Thus `v!=0` gives odd degree, while `v=0` is exactly the already closed
`j_1` ray.

Suppose `w!=0`.  The `y=0` specialization has degree seven with

```text
[x^7]R(x,0)=243w^2(u^2+216w).                               (42)
```

Off the cancellation seam it is nonsquare by odd degree.  On the seam

```text
w=-u^2/216,
```

one necessarily has `u!=0`.  The parameter `v` disappears from the
specialization, which becomes a sextic `F_u(x)`.  Compare it with
`(B_3x^3+B_2x^2+B_1x+B_0)^2` and saturate by introducing `z` with
`uz-1=0`.  The exact lexicographic Groebner basis of the seven coefficient
rows together with `uz-1` is

```text
[1].                                                         (43)
```

Therefore no nonzero constant triple `(u,v,w)` deforms `(20)` to another
square residual.

The theorem does **not** search polynomial-coefficient combinations
`f_1(x,y)j_1+f_2(x,y)j_2+f_3(x,y)j_3`, higher `J`-adic additions,
alternative lifts differing by multiples of `Delta`, or a polynomial-plane
Keller atlas.  Those are the precise remaining noncanonical directions.  No
Jacobian counterexample is claimed.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_three_cusp_polarization_branches_thm3872.py
python3 -O 04-computation/jc2_three_cusp_polarization_branches_thm3872.py
```

Both runs must byte-match
`05-knowledge/results/jc2_three_cusp_polarization_branches_thm3872.out`.
