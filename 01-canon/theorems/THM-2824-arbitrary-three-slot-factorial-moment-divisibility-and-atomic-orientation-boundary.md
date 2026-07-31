---
id: THM-2824
title: "Arbitrary three-slot factorial moment-three detection by tilted atomic orientation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For any three
  factorial slots,
  mean-zero polynomials form a binary plane.  Their second and third
  moments have a common complex projective zero exactly when two explicit
  real divisibility invariants vanish.  One invariant factors through a
  strictly positive cubic orientation and a telescoping sum of atomic
  determinants.  A positive hockey-stick expansion and a strict
  likelihood-ratio tilt prove every atom nonnegative, with equality
  exactly for one consecutive atom.  Consequently the first three
  factorial moments detect every arbitrary three-slot polynomial, and the
  associated two-charge Gaussian envelope with lowest slot a>=1 is
  detected by moment at most six.  Exact arithmetic through top exponent
  50 audits every formula, sign, and equality boundary.
source: root/audit-2809-2026-07-28
depends_on: []
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-2812-consecutive-three-slot-factorial-moment-six-detection
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - HYP-8765-gmc2-radial-channel-return-tower
script: 04-computation/gmc_arbitrary_three_slot_atomic_orientation_thm2824.py
output: 05-knowledge/results/gmc_arbitrary_three_slot_atomic_orientation_thm2824.out
script_sha256: ee30ff28881556cb00da6d80a085d5e104a02eaae55c14258a7f9ad47df2750d
output_sha256: 5046f21acf1030c725a10b53db4044107969645bf51e86262f290d64dfe35fb3
hash_basis: LF-normalized bytes
---

# THM-2824 -- arbitrary three-slot factorial moment-three detection

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2812 proves moment-three detection for three consecutive factorial
slots.  The first apparent extension fails if one chases signs in a
Euclidean remainder: on nonconsecutive supports the two real remainder
coefficients need not have a uniform termwise presentation.

The invariant replacement is cleaner.  The second moment is a positive
definite real binary quadratic.  A real binary cubic meets it over
`C P^1` exactly when the quadratic divides the cubic.  The two
division-free divisibility invariants expose a strictly oriented factor,
and one of them telescopes to atomic determinants.

The missing atomic inequality has a positive discrete model.  After a
hockey-stick expansion, its numerator and denominator are transforms of
two positive finite measures under the kernel
`H_n(j)=binom(n+j,j)`.  The numerator measure has weighted mean strictly
above `c-1`; the denominator is supported at or below `c-1`.  Successive
rows of `H` are increasing likelihood-ratio tilts, so the relevant
increment ratio is strictly increasing.  This proves the universal atomic
inequality and its exact equality case.  The proof and its independent
hostile audit are complete here.

## 1. Normalize the three slots

Let

```text
L:Q[s] -> Q,                           L(s^n)=n!,        (1)

f_n=s^n/n!,                            L(f_n)=1.         (2)
```

Fix arbitrary integers

```text
0<=a<b<c.                                               (3)
```

Every polynomial on these slots has a unique expression

```text
H=x f_a+y f_b+z f_c.                                   (4)
```

If `L(H)=0`, then `x+y+z=0`.  Put

```text
U=f_b-f_a,                         V=f_c-f_b.           (5)
```

Then there are unique `alpha,beta in C` with

```text
H=alpha U+beta V.                                      (6)
```

Define the real Gram and cubic tensors

```text
g11=L(U^2),       g12=L(UV),        g22=L(V^2),

t111=L(U^3),      t112=L(U^2 V),
t122=L(U V^2),    t222=L(V^3).                         (7)
```

Thus

```text
Q(alpha,beta)=L(H^2)
 =g11 alpha^2+2g12 alpha beta+g22 beta^2,              (8)

C(alpha,beta)=L(H^3)
 =t111 alpha^3+3t112 alpha^2 beta
  +3t122 alpha beta^2+t222 beta^3.                     (9)
```

For real `(alpha,beta)!=(0,0)`, the integral model gives

```text
Q(alpha,beta)
 =integral_0^infinity (alpha U(s)+beta V(s))^2 e^(-s) ds
 >0.                                                   (10)
```

Hence

```text
g11>0,        g22>0,        g11 g22-g12^2>0,           (11)
```

and `Q` is irreducible and squarefree over `R`.

## 2. Common nullity is binary quadratic divisibility

The two projective roots of `Q` are nonreal conjugates.  If the real cubic
`C` vanishes at one of them, it vanishes at the other.  Therefore

```text
Q and C have a common point in C P^1
 iff
Q divides C in R[alpha,beta].                          (12)
```

Write a possible quotient as `r alpha+s beta`.  Comparing the four
coefficients of

```text
C=(r alpha+s beta)Q                                   (13)
```

first forces

```text
r=t111/g11,                     s=t222/g22.            (14)
```

The two middle coefficients then vanish exactly when

```text
I1=
 3t112 g11 g22-t222 g11^2-2t111 g12 g22=0,            (15)

I2=
 3t122 g11 g22-2t222 g12 g11-t111 g22^2=0.            (16)
```

Consequently

```text
L(H)=L(H^2)=L(H^3)=0 for some H!=0 on {a,b,c}
 iff
I1=I2=0.                                               (17)
```

This is the first load-bearing simplification.  No sign comparison between
the two roots or between raw Euclidean-remainder coefficients is needed.  A
nonzero real linear remainder cannot meet the nonreal quadratic roots; the
division-free pair `(15)--(16)` says exactly when that remainder is zero.

## 3. The first cubic orientation is always strict

Put

```text
r=b-a,             X=f_b/f_a=s^r/(b!/a!),             (18)

F_k=L(f_a^3 X^k),                 k=0,1,2,3.
```

Then

```text
t111=F_3-3F_2+3F_1-F_0.                               (19)
```

The consecutive ratios are

```text
rho_k=F_(k+1)/F_k
     =product_(j=1)^r (3a+kr+j)/(a+j).                (20)
```

If `a=0`, then

```text
rho_0=1,
rho_1=binom(2r,r)>=2,
rho_2=binom(3r,r)>=3.
```

Therefore

```text
t111/F_0
 =rho_1(rho_2-3)+2>0.                                 (21)
```

If `a>0`, every factor in `rho_0` is greater than `1`, every factor in
`rho_1` is greater than `2`, and every factor in `rho_2` is at least `3`.
Hence

```text
t111/F_0
 =rho_0(rho_1(rho_2-3)+3)-1>0.                        (22)
```

Thus, for every `0<=a<b`,

```text
t111=L((f_b-f_a)^3)>0.                                (23)
```

This strict sign is independent of the third slot `c`.

## 4. The atomic orientation is universally nonnegative

Define

```text
D(a,b,c)=2t222 g12-3t122 g22.                         (24)
```

Equation `(16)` factors exactly as

```text
I2=-g11 D(a,b,c)-t111 g22^2.                          (25)
```

By `(11)` and `(23)`, the second term is strictly negative.  In particular,

```text
D(a,b,c)>=0                       implies I2<0,         (26)
```

which excludes common moment nullity without inspecting `I1`.

Now put

```text
d_i=f_(i+1)-f_i,                    0<=i<b.             (27)
```

Because

```text
U=sum_(i=a)^(b-1) d_i,                                (28)
```

and `(24)` is linear in its `U` occurrence,

```text
D(a,b,c)=sum_(i=a)^(b-1) D_i(b,c),                    (29)

D_i(b,c)
 =2L(V^3)L(d_i V)-3L(d_i V^2)L(V^2).                 (30)
```

Equivalently, each atom is the exact determinant

```text
D_i(b,c)
 =det [[L(d_i V),   L(d_i V^2)],
       [3L(V^2),    2L(V^3)]].                        (31)
```

We now prove the sharp sign of `(31)`.

### 4.1. Turn the global moments into one long secant

For `n>=0`, put

```text
h1(n)=L(f_n V),                 h2(n)=L(f_n V^2),

A_n=h1(n+1)-h1(n)=L(d_n V),
B_n=h2(n+1)-h2(n)=L(d_n V^2),
R_n=B_n/A_n.                                             (32)
```

Since `b>=1`, one has `V(0)=0`.  Also

```text
V'=f_(c-1)-f_(b-1).                                     (33)
```

Integration by parts against `e^(-s) ds` gives

```text
L(V^2)=2L(VV')
      =2 sum_(n=b-1)^(c-2) A_n,

L(V^3)=3L(V^2 V')
      =3 sum_(n=b-1)^(c-2) B_n.                         (34)
```

Thus the global ratio in `(31)` is a positive weighted average of the
increment ratios `R_n`.  It remains to prove that `R_n` is strictly
increasing.

### 4.2. A positive hockey-stick model

Set

```text
K_n(m)=L(d_n f_m)=binom(n+m,m-1),
H_n(j)=binom(n+j,j).                                    (35)
```

The hockey-stick identity says

```text
K_n(m)=sum_(j=0)^(m-1) H_n(j).                          (36)
```

Consequently

```text
A_n=K_n(c)-K_n(b)
   =sum_(j=b)^(c-1) H_n(j)>0.                           (37)
```

Write

```text
alpha=binom(2b,b),       beta=binom(b+c,b),
gamma=binom(2c,c),       r=c-b.                         (38)
```

The divided-power product law gives

```text
B_n=alpha K_n(2b)-2beta K_n(b+c)+gamma K_n(2c)
   =sum_(j=0)^(2c-1) q_j H_n(j),                        (39)
```

where

```text
q_j=
 alpha-2beta+gamma,                 0<=j<2b,
 gamma-2beta,                      2b<=j<b+c,
 gamma,                           b+c<=j<2c.             (40)
```

These are nonnegative weights.  Indeed,

```text
gamma/beta
 =product_(k=1)^r (2b+r+k)/(b+k)
 >=2^r>=2.                                               (41)
```

The middle weight in `(40)` is therefore nonnegative, and the first is
strictly positive because it is the middle weight plus `alpha`.  In
particular `B_n>0` and `R_n` is well-defined.

### 4.3. Strict separation of the tilted means

Let `mu_q(n)` be the `j`-mean of the positive weights
`q_j H_n(j)`, and let `mu_A(n)` be the `j`-mean of
`H_n(j)` on `b<=j<c`.  The latter support immediately gives

```text
mu_A(n)<=c-1.                                            (42)
```

At `n=0`, `H_0(j)=1`.  Summing each constant block in `(40)` gives the
exact centered first moment

```text
sum_(j=0)^(2c-1) (j-(c-1))q_j
 =c gamma+(r-1)(b+c)beta-(2r-1)b alpha.                 (43)
```

This is strictly positive.  Successive central-binomial ratios give

```text
gamma/alpha
 =product_(t=1)^r [4-2/(b+t)]
 >=3^r>2r-1,                                             (44)
```

because `b+t>=2`; the middle term in `(43)` is nonnegative and `c>b`.
Therefore

```text
mu_q(0)>c-1.                                             (45)
```

The kernel rows satisfy the exact likelihood-ratio tilt

```text
H_(n+1)(j)/H_n(j)=(n+j+1)/(n+1).                        (46)
```

If `Var_q(n)` is the variance under the normalized weights
`q_jH_n(j)`, then

```text
mu_q(n+1)-mu_q(n)
 =Var_q(n)/(n+1+mu_q(n))>0.                             (47)
```

Strictness follows already from the positive weights at `j=0,1`.
Combining `(42)`, `(45)`, and `(47)` yields

```text
mu_q(n)>c-1>=mu_A(n)                  for every n>=0.    (48)
```

Applying `(46)` once to the total masses in `(37)` and `(39)` now gives

```text
R_(n+1)/R_n
 =(n+1+mu_q(n))/(n+1+mu_A(n))>1.                        (49)
```

Thus `R_n` is strictly increasing for every `n>=0`, a statement stronger
than the range `n<b` needed here.

### 4.4. Atomic sign and equality

Let

```text
Rbar=
 [sum_(n=b-1)^(c-2) A_n R_n]/
 [sum_(n=b-1)^(c-2) A_n].                               (50)
```

Equations `(30)`, `(32)`, and `(34)` give the exact factorization

```text
D_i(b,c)
 =6 A_i [sum_(n=b-1)^(c-2) A_n] [Rbar-R_i].             (51)
```

For `i<b-1`, strict monotonicity gives
`R_i<R_(b-1)<=Rbar`.  For `i=b-1`, the average is equal to
`R_(b-1)` exactly when it has a single term, namely when `c=b+1`.
Therefore

```text
D_i(b,c)>=0                         for every 0<=i<b<c,

D_i(b,c)=0
 iff
(i,b,c)=(j,j+1,j+2)                for some j>=0.        (52)
```

Now `(29)` gives `D(a,b,c)>=0`.  Equations `(25)--(26)` give `I2<0`,
and `(17)` proves the universal three-slot theorem

```text
L(H)=L(H^2)=L(H^3)=0                 implies H=0         (53)
```

for every arbitrary three-slot factorial support.  The atomic statement
is a sufficient certificate for `(53)`; no converse is asserted.

## 5. A continuous Chebyshev-orientation sidecar

The adjacent difference pair already has a strict pointwise orientation.
Put

```text
p=b-a,       r=c-b,
A=b!/a!,     B=c!/b!,

X=s^p/A,     Y=s^r/B.                                  (54)
```

Then

```text
U=f_a(X-1),                  V=f_a X(Y-1),
```

and direct differentiation gives

```text
U V'-U' V
 =f_a^2 X/s [p+rXY-(p+r)Y].                            (55)
```

Every factor in `B` is strictly larger than every factor in `A`, so

```text
B^(1/r)>A^(1/p),                    X^r>Y^p.            (56)
```

Weighted AM--GM now gives

```text
p+rXY
 >=(p+r)(XY)^(r/(p+r))
 >(p+r)Y.                                             (57)
```

Thus

```text
U(s)V'(s)-U'(s)V(s)>0                  for every s>0.   (58)
```

This is a genuine extended-Chebyshev orientation, but it is logically
independent of the discrete proof above.  Pointwise orientation alone does
not transparently control `(31)`, because the atomic determinant mixes the
signed functions `d_i V` and `d_i V^2` against different global moments.
The hockey-stick likelihood-ratio coordinate is the sidecar that performs
that transfer.

## 6. Finite exact audit through exponent 50

The exact companion exhausts all

```text
0<=i<b<c<=50.                                          (59)
```

There are

```text
binom(51,3)=20,825
```

such atomic triples.  It finds

```text
D_i(b,c)>=0                         in all 20,825 cases, (60)
```

with exactly `49` zeros:

```text
(i,b,c)=(j,j+1,j+2),                  0<=j<=48.         (61)
```

The smallest positive value in the normalized factorial basis is

```text
D_0(2,3)=288.                                           (62)
```

For every `0<=a<b<c<=50`, the companion separately verifies the Gram
determinant, `(19)--(25)`, the telescoping identity `(29)`, and

```text
I1<0,                 I2<0.                            (63)
```

The companion additionally checks `(35)--(51)` in exact arithmetic:
the hockey-stick decomposition, positivity of all three `q` blocks,
the centered-mean identity, strict mean separation, the row-tilt
recursion, the ratio increment, and the final atomic factorization.
This is a bounded audit of the universal proof, not the source of its
quantifiers.  It uses integer multinomial moments only; there is no
floating-point sign decision.

## 7. Two-charge Gaussian consequence

Assume `a>=1` and write

```text
h=H/s,
P=W+Z h(ZW),                         W=conj(Z),         (64)
```

for a standard complex Gaussian `Z`.  Charge balance gives, for
`j=1,2,3`,

```text
E[P^(2j)]=binom(2j,j)L(H^j),          E[P^(2j-1)]=0.    (65)
```

Hence `(53)` proves Gaussian detection by moment at most six for every
arbitrary three-slot two-charge envelope `(64)`.

For exact three-slot `H`, the polynomial `P` has four monomials and primitive
return `R=2`, so six is the HYP-8765 cutoff `(k-1)R`.  This statement does
not separate unrelated Wick channels in a more general polynomial.

## 8. Information and failure ledger

| item | exact content |
|---|---|
| source | arbitrary normalized factorial slots `f_a,f_b,f_c` |
| mean-zero chart | `H=alpha(f_b-f_a)+beta(f_c-f_b)` |
| exact obstruction | simultaneous vanishing of `I1,I2` |
| strict invariant | `t111>0` for every `a<b` |
| atomic map | `D=sum D_i`, with `I2=-g11D-t111g22^2` |
| positive coordinate | hockey-stick weights `q_j>=0` |
| strict mechanism | numerator mean `>c-1`, denominator support `<=c-1` |
| equality | only `D_j(j+1,j+2)=0` |
| universal conclusion | arbitrary three-slot detection by moments `1,2,3` |
| continuous sidecar | strict Wronskian `(58)` |
| finite audit | all `c<=50`; `20,825` exact atomic checks |
| audit status | independently rederived; exact hostile controls through `c=1000` |

This theorem proves the arbitrary three-slot Strong Factorial statement and
the corresponding two-charge moment-six bound.  It does not prove general
HYP-8765, separate unrelated Wick channels, or give a new proof of all GMC2.

## 9. Independent hostile audit

The independent audit rederived, without importing the companion:

1. the three hockey-stick weight blocks and their nonnegativity;
2. the centered identity `(43)` and the product bounds `(41),(44)`;
3. the strict variance and affine likelihood-ratio recursions
   `(46)--(49)`;
4. the integration-by-parts long secant `(34)`;
5. the factorization `(51)` and its unique equality boundary; and
6. the binary quadratic/cubic divisibility conclusion.

A separate exact reconstruction checked `7,228` support pairs and
`310,934` atoms, including deterministic and random hostile edges through
`c=1000`.  It agreed with every sign, factorization, and equality claim.
As an orthogonal boundary control, the adjacent family `c=b+1` also admits
a direct descending-ratio proof; its exact formulas were independently
replayed for all `b<=100`.

## 10. Exact companion

Run

```text
python 04-computation/gmc_arbitrary_three_slot_atomic_orientation_thm2824.py
python -O 04-computation/gmc_arbitrary_three_slot_atomic_orientation_thm2824.py
```

Both executions byte-match

```text
05-knowledge/results/gmc_arbitrary_three_slot_atomic_orientation_thm2824.out.
```

The dependency-free companion uses exact integers and fractions.  It checks
the Gram, divisibility, strict-cubic, Wronskian, hockey-stick, positive-weight,
tilted-mean, ratio-monotonicity, telescoping, and atomic-factorization
identities on all `20,825` triples through `c=50`.  It has explicit
exception gates, no truth-bearing Python assertions, no floating point, and
no scratch dependency.

**QED.**
