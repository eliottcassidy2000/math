---
id: THM-2855
title: "Shifted positive-cone transverse asymptotic family"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  integer n>=1 there is exactly one positive
  pair (x_n,y_n) for which the quadratic factorial moment on
  span{d_n+x_n d_(n+2), d_(n+1)+y_n d_(n+2)} divides the cubic.
  The intersection is transverse.  The branch has an exact real-analytic
  expansion at infinity, and its fourth factorial moment escapes for all
  sufficiently large n.  All-depth fourth-moment escape is not claimed.
source: root/gmc-shifted-transverse-family-2026-07-28
depends_on: []
related:
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
  - THM-2848-whitened-moving-plane-multipole-and-pearson-boundary
  - THM-2854-four-slot-shifted-window-macaulay-box
script: 04-computation/gmc_shifted_positive_cone_transverse_family_thm2855.py
output: 05-knowledge/results/gmc_shifted_positive_cone_transverse_family_thm2855.out
script_sha256: d2db5c0fd1051cfcd11eb23e504accc2a3d5602f5e38e9cda468c659aa8396a6
output_sha256: 3028e4e777f36d5c74fe6581f2e736268a94bffee404611ce491bb6423bd2180
hash_basis: LF-normalized bytes
---

# THM-2855 -- shifted positive-cone transverse asymptotic family

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Put

```text
L(s^m)=m!,                f_m=s^m/m!,
d_m=f_(m+1)-f_m.                                      (1)
```

For every integer `n>=1`, there is exactly one pair

```text
x_n>0,                  y_n>0                         (2)
```

such that, for

```text
U_n=d_n+x_n d_(n+2),       V_n=d_(n+1)+y_n d_(n+2),   (3)
```

the binary quadratic

```text
Q_n(alpha,beta)=L((alpha U_n+beta V_n)^2)             (4)
```

divides the corresponding binary cubic

```text
C_n(alpha,beta)=L((alpha U_n+beta V_n)^3).             (5)
```

The common intersection is transverse.  Consequently either projective
root `(alpha:beta)` of `Q_n` gives a nonzero polynomial

```text
H_n=alpha U_n+beta V_n
```

with

```text
L(H_n)=L(H_n^2)=L(H_n^3)=0.                            (6)
```

This extends the single `n=1` hostile in THM-2846 to an infinite shifted
family.  It remains a finite-moment hostile, not a counterexample to GMC2.

## 1. Normalize the exact tensors

For offsets `0<=a,b,c<=2`, define

```text
B_ab(n)=L(d_(n+a)d_(n+b))/L(d_(n+2)^2),
T_abc(n)=L(d_(n+a)d_(n+b)d_(n+c))/L(d_(n+2)^3).       (7)
```

The adjacent-difference identities give

```text
L(d_i d_j)=binom(i+j,i),                               (8)

L(d_i d_j d_k)
 =(i+j+k)!/(i!j!k!)
  *[2(i+j+k+1)^2+(i+j+k)(ij+ik+jk)-ijk]
   /[(i+1)(j+1)(k+1)].                                (9)
```

Thus `(7)` consists of rational functions of `n` whose denominators are
positive for `n>=1`.  If `t=1/n`, each extends analytically to `t=0`, and

```text
B_ab(infinity)=2^(a+b-4),
T_abc(infinity)=3^(a+b+c-6).                          (10)
```

Use `(7)` to form

```text
g11=L(U^2),     g12=L(UV),      g22=L(V^2),
t111=L(U^3),    t112=L(U^2V),
t122=L(UV^2),   t222=L(V^3),                            (11)
```

with the common positive normalizing factors suppressed.  As in THM-2824,
quadratic divisibility is equivalent to

```text
I1=3t112 g11 g22-t222 g11^2-2t111 g12 g22=0,
I2=3t122 g11 g22-2t222 g12 g11-t111 g22^2=0.          (12)
```

## 2. Exact all-integer elimination

Clear the positive denominators in `(12)`, obtaining

```text
P1(n,x,y),                  P2(n,x,y) in Z[n,x,y],
deg_x(P1,P2)=(4,3).                                    (13)
```

Their exact subresultant degree profile in `x` is

```text
4,3,2,1,0.                                             (14)
```

The last member factors as

```text
Res_x(P1,P2)
 =R(n) G_n(y)^2 p_n(y),                                (15)
```

up to a nonzero rational integer, where every factor of `R(n)` is
strictly positive for `n>=1`,

```text
G_n(y)
 =(n+2)+2(2n+3)y+2(2n+3)y^2>0                         (16)
```

for `y>0`, and `p_n` has degree fifteen.

Write the coefficients of `p_n`, in descending powers of `y`, as
`c_0(n),...,c_15(n)`.  The first three have only positive coefficients as
polynomials in `n`.  Each later `c_j` has exactly one sign variation as a
polynomial in `n`, and its unique positive zero lies in

```text
j              3   4   5   6   7   8  9 10 11 12 13 14 15
root interval  (101,102),(42,43),(24,25),(15,16),(10,11),
               (7,8),(5,6),(4,5),(3,4),(2,3),
               (1,2),(1,2),(0,1).                    (17)
```

The right endpoints in `(17)` are nonincreasing.  Therefore, for every
integer `n>=1`, the coefficient list

```text
c_0(n),...,c_15(n)
```

has exactly one sign change.  Its leading coefficient is positive and its
constant coefficient is negative.  Descartes' rule and the intermediate
value theorem give exactly one positive root `y_n`, counted with
multiplicity.  In particular that root is simple.

After its common content is removed, the linear subresultant is

```text
A_n(y)x+C_n(y),                                        (18)
```

where all `251` coefficients of `A` as a polynomial in `(n,y)` are
positive and all `253` coefficients of `C` are negative.  Hence

```text
x_n=-C_n(y_n)/A_n(y_n)>0.                              (19)
```

The nonzero linear subresultant proves that `(19)` is the unique common
`x`-root.  The factor `p_n` occurs to exponent one in `(15)`, its positive
root is simple, and the `x`-leading coefficients in `(13)` never vanish
on `n>=1,x>0,y>0`.  The intersection is therefore transverse.  This also
rules out a second positive branch at another scale.

## 3. Blow up the singular limiting line

At `t=0`, the normalized invariants are

```text
I1(0,x,y)
 =-(4x+1)(6x-5y-1)^2(108xy+48x+17y+7)/186624,

I2(0,x,y)
 =-(2y+1)(6x-5y-1)^2(54xy+21x+11y+4)/46656.           (20)
```

Thus the naive limiting Jacobian vanishes on the whole line

```text
6x-5y-1=0.                                             (21)
```

The coefficient of `t` after imposing `(21)` has the common factor

```text
(4x+1)^2(9x+1)(18x^2-36x-7).                          (22)
```

The unique positive choice is

```text
x_0=1+5sqrt(2)/6,              y_0=1+sqrt(2).          (23)
```

Now put

```text
x=x_0+a t,                     y=y_0+b t,
D=6a-5b.                                               (24)
```

There are nonzero constants

```text
p1=25(-99sqrt(2)-140)/5184,        p2=(6/5)p1
```

such that

```text
t^-2 I1=p1 F1+O(t),             t^-2 I2=p2 F2+O(t),    (25)
```

where, in coordinates `(A,D)=(a,6a-5b)`,

```text
36F1=
 D^2+(12+38sqrt(2)/5)D-(2sqrt(2)/5)A
 -829sqrt(2)/45-13867/540,

36F2=
 D^2+(12+23sqrt(2)/3)D-(2sqrt(2)/5)A
 -827sqrt(2)/45-13871/540.                            (26)
```

The crucial triangular difference is

```text
F1-F2=-(sqrt(2)/540)(D+2/3-sqrt(2)/18).               (27)
```

Thus the unique blown-up zero is

```text
D_0=-2/3+sqrt(2)/18,

a_0=-(74184+52463sqrt(2))/1296,
b_0=-(14808+10495sqrt(2))/216.                        (28)
```

Direct differentiation gives

```text
det d(F1,F2)/d(a,b)=1/4860.                            (29)
```

Equivalently the triangular `(A,D)` determinant has magnitude `1/24300`;
the factor five is the determinant of the coordinate change in `(24)`.
The real-analytic implicit-function theorem now supplies a unique
real-analytic branch for all sufficiently small positive `t`.  It agrees
with the all-integer solution of Section 2, and

```text
x_n=
 1+5sqrt(2)/6
 -(74184+52463sqrt(2))/(1296n)+O(n^-2),

y_n=
 1+sqrt(2)
 -(14808+10495sqrt(2))/(216n)+O(n^-2).                (30)
```

This is the load-bearing blow-up.  Applying the implicit-function theorem
directly to `(20)` would be invalid because of the squared limiting line.

## 4. Eventual fourth-moment exit

Let

```text
F_n(alpha,beta)=L((alpha U_n+beta V_n)^4).
```

Reduce `F_n(1,z)` modulo `Q_n(1,z)`.  If the coefficient of `z` in the
remainder is cleared by `g22^3` and the positive fourth-multinomial scale,
its limit along `(30)` is

```text
3(208885+147704sqrt(2))/262144>0.                      (31)
```

The normalized Gram determinant simultaneously satisfies

```text
g11g22-g12^2
 =(17+12sqrt(2))/(1152n)+O(n^-2)>0.                   (32)
```

At the upper-half-plane root of `Q_n`, the imaginary part of the fourth
moment is the product of the positive coefficient in `(31)` and the
positive imaginary part supplied by `(32)`, divided by a positive power
of `g22`.  Hence

```text
L(H_n^4)!=0
```

for all sufficiently large `n`; the conjugate root gives the conjugate
nonzero value.  Since `n>=1`, write `H_n=s h_n` and set

```text
P_n(Z,W)=W+Z h_n(ZW),                 W=conj(Z).        (33)
```

Charge balance gives

```text
E[P_n^(2k)]=binom(2k,k)L(H_n^k),
E[P_n^(2k-1)]=0.                                      (34)
```

Thus the associated two-charge Gaussian envelope has moments through six
zero and its eighth moment nonzero at every sufficiently large depth.

This section is deliberately **eventual**.  The exact elimination in
Section 2 proves the cubic-null branch for every integer `n>=1`, but
`(31)` alone does not prove fourth-moment escape at every finite depth.
THM-2846 supplies the `n=1` control, and THM-2849 covers the first-window
supports through its stated finite range.  A uniform all-depth quartic
claim would still need either an exact resultant/nonvanishing argument or
an explicit asymptotic threshold plus a finite certified bridge.

## 5. Scope and evidence

The exact companion:

1. reconstructs the bilinear and cubic tensors independently from the
   `2^k` adjacent-difference expansion at literal depths `1,2,3`;
2. verifies every limiting, blow-up, triangular, and Jacobian identity;
3. computes the full symbolic subresultant chain, factors the resultant,
   and certifies the coefficient-threshold staircase `(17)`;
4. proves coefficientwise positivity/negativity of the linear selector
   `(18)`; and
5. derives the exact Gram and quartic limits `(31)--(32)`.

All truth-bearing checks use exact arithmetic and explicit exceptions, so
optimized execution retains the proof gates.  An independent hostile audit
rederived the adjacent-difference cubic tensor, the Sturm sign staircase,
the subresultant and selector signs, transversality, every blow-up constant
and Jacobian, and the eventual quartic/Gaussian exit.  Independent normal
and optimized runs byte-match the stored transcript and both declared LF
hashes.  No arbitrary-cone, all-support, all-depth quartic, SFC(4), or GMC2
conclusion is claimed.

**QED.**
