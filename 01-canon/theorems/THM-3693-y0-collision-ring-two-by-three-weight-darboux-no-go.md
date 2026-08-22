---
id: THM-3693
title: "Two-by-three weight Darboux no-go in the y=0 collision ring"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the y=0
  collision ring R of THM-3686, no Darboux pair {P,Q}=1 can have at most two
  active grading weights in one output and at most three in the other.  For
  nonarithmetic supports, a sumset classification plus Poisson-centralizer
  transitivity kills the unique duplicated bucket.  Arithmetic supports
  have two possible middle constant buckets.  A complete all-integer case
  split handles nonpositive and weight-zero edges, endpoint valuations force
  one negative weight to be -2, and four explicit first-order ODEs factor
  the alleged unit through a nonunit endpoint Wronskian.  The two arithmetic
  buckets are proved separately; no false h/u reflection is used.  Hence
  every counterexample survivor in R needs at least three active weights in
  each output.  Three-by-three noncentered supports, JC(2), and arbitrary
  quartic C3 data remain OPEN.
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  PASS -- root independently rederived the nonarithmetic centralizer case,
  both arithmetic support parameterizations, all nonpositive/weight-zero
  edges, the endpoint min(a,b)=2 reduction, and all four ODE/factor families,
  including both gcd-three cube-root branches.
depends_on:
  - THM-3686-y0-collision-normalization-and-bracket-anatomy
  - THM-3691-y0-collision-ring-two-weight-darboux-no-go
related:
  - THM-3692-ordinary-two-by-three-shift-root-peeling
script: 04-computation/jacobian_y0_two_by_three_weight_no_go_thm3693.py
output: 05-knowledge/results/jacobian_y0_two_by_three_weight_no_go_thm3693.out
script_sha256: 9727857098803f3a615f8ab2ff7aad1aead11d724f8932332ca83898eeb08072
output_sha256: 95a9fa1b59ef5c4496dc00767c4650b1df38794dbb8ed5f4f8e76a0176d52758
hash_basis: LF-normalized bytes
---

# THM-3693 -- two-by-three grading support cannot carry the collision

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Retain the collision ring and grading of THM-3686/3691:

```text
R=C[A,B,C] subset C[x,z],
u=1-x^2z,                    h=1-u^2,
A=3x^-2h,                    B=2x^-1h,        C=xu.    (1)
```

Every weight-`r` element is `x^r p(u)`, and

```text
{x^r p(u),x^s q(u)}
 =x^(r+s+1)[s p'q-r p q'].                             (2)
```

Moreover,

```text
r=-R<0 => h^ceil(R/2)|p,          r=R>0 => u^R|p.      (3)
```

We prove that `{P,Q}` is not a nonzero constant whenever

```text
|supp_wt(P)|<=2,                  |supp_wt(Q)|<=3,      (4)
```

or with the two outputs reversed.  Scalar components are bracket-invisible
and are not active weights.  THM-3691 already handles every proper subcase,
so assume the support sizes are exactly two and three.

## 1. A Poisson-centralizer transitivity lemma

Let `K=C(x,u)`.  If nonconstant `F,G,T in K` satisfy

```text
{F,G}=0={F,T},                                           (5)
```

then `{G,T}=0`.  Indeed, nondegeneracy of the two-variable Poisson bracket
makes `(5)` equivalent to

```text
dF wedge dG=0=dF wedge dT.                              (6)
```

In the two-dimensional `K`-space of rational differentials, `dG` and `dT`
are both multiples of the nonzero vector `dF`, so their wedge is zero.

We also retain the homogeneous commutation consequences of THM-3691.  Two
strictly opposite nonzero weights cannot commute.  A commuting weight-zero
component paired with a nonzero weight is scalar.  Same-sign commuting
components are powers of one common coefficient polynomial, with exponents
fixed by the gcd of their weights.

## 2. The two-by-three sumset grammar

Write the supports by their gaps:

```text
supp(P)={r,r+d},
supp(Q)={s,s+e,s+e+f},                  d,e,f>0.         (7)
```

A constant needs at least two bracket edges in one output-weight bucket;
a unique edge would be a forbidden homogeneous unit bracket.  Equality of
two edge sums requires

```text
d in {e,f,e+f}.                                         (8)
```

If the supports are not arithmetic with `d=e=f`, exactly one bucket is
duplicated.  Every other edge is unique and therefore commutes.  This already
kills the duplicate:

1. if `d=e`, the top component of `Q` commutes with both components of `P`;
2. if `d=f`, the bottom component of `Q` does; and
3. if `d=e+f`, the middle component of `Q` does.

Centralizer transitivity then makes the two `P` components commute.  Each
remaining `Q` component already commutes with one of them, so a second use of
transitivity makes both edges in the duplicated bucket commute as well.
They cannot sum to one.

Thus only the arithmetic supports remain:

```text
supp(P)={r,r+d},               supp(Q)={s,s+d,s+2d}.   (9)
```

There are two duplicated middle buckets, and the constant can occur in
either.  We prove them separately because the coefficient modules in `(3)`
are asymmetric at `h=0` and `u=0`.

## 3. First arithmetic middle bucket

If the first middle bucket is constant, parameterize the weights as

```text
supp(P)={-a,b-1},
supp(Q)={-b,a-1,2a+b-2},               a+b>1.          (10)
```

The extreme brackets vanish, the first middle bucket equals one, and the
second middle bucket cancels:

```text
{P_-a,Q_-b}=0,
{P_-a,Q_(a-1)}+{P_(b-1),Q_-b}=1,
{P_-a,Q_(2a+b-2)}+{P_(b-1),Q_(a-1)}=0,
{P_(b-1),Q_(2a+b-2)}=0.                                (11)
```

### 3.1 Integer and endpoint reductions

If `a<=0` or `b<=0`, one extreme contains strictly opposite nonzero weights,
or a weight-zero component is forced scalar.  The result reduces to
THM-3691.  If `a=1`, lower-extreme commutation gives coefficient cores
`H,H^b`; both constant edges vanish at `h=0`.  The case `b=1` is symmetric.
Hence `a,b>=2`.

The constant bucket in `(11)` contains the negative components of weights
`-a` and `-b`.  If both `a,b>=3`, `(3)` makes their coefficient polynomials
and their derivatives vanish at `u=+1,-1`; the alleged unit vanishes there.
Thus

```text
min(a,b)=2.                                              (12)
```

### 3.2 The branch `a=2`

If `b` is odd, lower-extreme commutation gives cores `H^2,H^b`.  Since the
weight-`-2` core is divisible by the squarefree `h`, one has `h|H`, and both
constant edges vanish at `h=0`.  Hence put `b=2m`.

The lower cores are now `H,H^m`.  For the upper extreme set

```text
epsilon=gcd(b-1,b+2) in {1,3},
k=(b-1)/epsilon,                   n=(b+2)/epsilon.     (13)
```

Its cores are `K^k,K^n`, with `u^epsilon|K`.  Write the weight-one middle
coefficient as `L`.  The cancelling bucket in `(11)` is a first-order ODE.
After absorbing the nonzero extreme scalars, its polynomial solutions are

```text
epsilon=1:  L=kappa K+rho H K^3,

epsilon=3:  L=rho H K,
            plus kappa T only when K=T^3.               (14)
```

Here `rho` is the forced nonzero multiple of `n/k`; its exact value is
irrelevant to the factor below.  Substitution into the constant bucket gives

```text
epsilon=1:   (H'K+2HK') Phi(H,K)=1,

epsilon=3, no cube root:
             (3H'K+2HK') Phi(H,K)=1,

epsilon=3, K=T^3:
             (H'T+2HT') Psi(H,T)=1.                    (15)
```

All displayed `Phi,Psi` are polynomials.  The first and third leading
Wronskians are nonconstant because `h|H` and `u|K` or `u|T`.  In the middle
line `u^3|K`, so `3H'K+2HK'` is divisible by `u^2`.  None is a unit.

### 3.3 The branch `b=2`

If `a` is odd, lower commutation again upgrades the endpoint order and the
constant bucket vanishes at `h=0`.  Put `a=2m`.  Normalize the four extreme
cores as

```text
P_-a: H^m,       Q_-2: H,
P_1: K,          Q_(2a): K^(2a).                       (16)
```

If `S` is the coefficient of `Q_(a-1)`, the cancelling-bucket ODE has

```text
S=kappa K^(a-1)+rho H^m K^(2a-1).                      (17)
```

The forced particular coefficient is a nonzero multiple of `2a`.  Exact
substitution factors the constant equation as

```text
(H'K+2HK') Phi(H,K)=1,                                 (18)
```

again impossible.  This closes the first arithmetic bucket.

## 4. Second arithmetic middle bucket

If the second middle bucket is constant, the correct parameterization is

```text
supp(P)={-a,b-1},
supp(Q)={1-a-2b,-b,a-1}.                               (19)
```

The equations are now

```text
{P_-a,Q_(1-a-2b)}=0,
{P_-a,Q_-b}+{P_(b-1),Q_(1-a-2b)}=0,
{P_-a,Q_(a-1)}+{P_(b-1),Q_-b}=1,
{P_(b-1),Q_(a-1)}=0.                                   (20)
```

Nonpositive parameters reduce by opposite-weight or weight-zero extremes.
If `a=1<b`, the upper weight-zero component is scalar; if `b=1<a`, the
weight-zero component of `P` is scalar; and if `a=b=1`, the constant bucket
is divisible by `h`.  Thus `a,b>=2`, and the same `h=0` evaluation of the
constant bucket forces `(12)`.

### 4.1 The branch `a=2`

The four extreme cores and the unknown weight-`-b` coefficient `L` are

```text
P_-2: H^2,             Q_-(2b+1): H^(2b+1),
P_(b-1): K^(b-1),      Q_1: K,
Q_-b: L.                                                 (21)
```

The first cancelling bucket in `(20)` integrates exactly to

```text
L=kappa H^b+rho H^(2b-1)K^(b-1),                       (22)
```

where the normalized particular coefficient is `(2b+1)/2`.  The unit
bucket then has the stronger factorization

```text
H(H'K+HK') Phi(H,K)=1.                                 (23)
```

The visible factor `H` is already a nonunit.

### 4.2 The branch `b=2`

Put

```text
epsilon=gcd(a,a+3) in {1,3},
k=a/epsilon,                       n=(a+3)/epsilon.     (24)
```

The extreme cores are `H^k,H^n,K,K^(a-1)`.  If `L` is the unknown
weight-`-2` coefficient, the cancelling ODE gives

```text
L=rho H^(3/epsilon)K,

plus kappa H^2 when epsilon=1,
or plus kappa T^2 when epsilon=3 and H=T^3.             (25)
```

The normalized particular coefficient is `n/k`.  Substitution gives

```text
epsilon=1:   (H'K+HK') Phi(H,K)=1,

epsilon=3, no cube root:
             (H'K+3HK') Phi(H,K)=1,

epsilon=3, H=T^3:
             (T'K+TK') Psi(T,K)=1.                    (26)
```

For `epsilon=1`, divisibility gives `h|H,u|K`, so the leading factor is
nonconstant.  For `epsilon=3`, the weight-`-a` condition on `H^k` forces
`h^2|H`; hence `H'K+3HK'` is divisible by `h`.  On the cube branch it forces
`h|T`, and `T'K+TK'` is again nonconstant.  This closes the second bucket.

## 5. Consequence and exact frontier

Sections 2--4 prove, in every target degree,

```text
min(|supp_wt(P)|,|supp_wt(Q)|)<=2 and
max(|supp_wt(P)|,|supp_wt(Q)|)<=3
              => {P,Q} notin C*.                       (27)
```

Therefore any Darboux pair in the collision ring must satisfy

```text
|supp_wt(P)|>=3,                    |supp_wt(Q)|>=3.    (28)
```

Such a pair would still identify the THM-3686 collision and refute `JC(2)`.
The theorem eliminates the entire two-by-three grading frontier, not all
three-by-three supports.  In particular noncentered `3x3`, larger supports,
general quartic `C3` incidence data, and `JC(2)` remain open.

## 6. Exact reproduction

Run

```bash
python3 -B 04-computation/jacobian_y0_two_by_three_weight_no_go_thm3693.py
python3 -B -O 04-computation/jacobian_y0_two_by_three_weight_no_go_thm3693.py
```

Both modes byte-match the stored transcript.  The companion checks the
sumset classification on a hostile integer window and symbolically verifies
all four ODE solutions/factorizations, including both gcd-three cube-root
branches.  The support and endpoint arguments above are all-degree proofs.

**QED.**
