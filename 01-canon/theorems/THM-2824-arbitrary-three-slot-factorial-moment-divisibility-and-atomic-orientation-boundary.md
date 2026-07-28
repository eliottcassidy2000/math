---
id: THM-2824
title: "Arbitrary three-slot factorial divisibility and atomic orientation boundary"
status: >
  RESERVED / PROOF-COMPLETE REDUCTION CANDIDATE + VERIFIED-EXACT /
  AWAITING INDEPENDENT HOSTILE AUDIT.  For any three factorial slots,
  mean-zero polynomials form a binary plane.  Their second and third
  moments have a common complex projective zero exactly when two explicit
  real divisibility invariants vanish.  One invariant factors through a
  strictly positive cubic orientation and a telescoping sum of atomic
  determinants.  Universal nonnegativity of those atoms would prove
  arbitrary three-slot Strong Factorial detection and the associated
  two-charge Gaussian moment-six bound.  The atomic inequality is OPEN.
  Exact exhaustive arithmetic proves it for every 0<=i<b<c<=50, with
  equality only on consecutive atoms, and therefore closes every
  three-slot support whose top exponent is at most 50.
source: root/audit-2809-2026-07-28
depends_on: []
related:
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-2812-consecutive-three-slot-factorial-moment-six-detection
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - HYP-8765-gmc2-radial-channel-return-tower
script: 04-computation/gmc_arbitrary_three_slot_atomic_orientation_thm2824.py
output: 05-knowledge/results/gmc_arbitrary_three_slot_atomic_orientation_thm2824.out
script_sha256: ebc659e273db7b956cdd536001ab91b3f02445cf7fea050186df76b2c9c762ef
output_sha256: ba0be7465d0842d398bedf06217d47404154ac1454e9bc440b6657f3948c8c16
hash_basis: LF-normalized bytes
---

# THM-2824 -- arbitrary three-slot nullity reaches one atomic inequality

**RESERVED / PROOF-COMPLETE REDUCTION CANDIDATE + VERIFIED-EXACT /
AWAITING INDEPENDENT HOSTILE AUDIT.**

THM-2812 proves moment-three detection for three consecutive factorial
slots.  The first apparent extension fails if one chases signs in a
Euclidean remainder: on nonconsecutive supports the two real remainder
coefficients need not have a uniform termwise presentation.

The invariant replacement is cleaner.  The second moment is a positive
definite real binary quadratic.  A real binary cubic meets it over
`C P^1` exactly when the quadratic divides the cubic.  The two
division-free divisibility invariants expose a strictly oriented factor,
and one of them telescopes to a single atomic inequality.

The reduction and all orientation statements below are proved.  The
universal atomic inequality `(24)` is explicitly **OPEN**.  Its exact
verification through top exponent `50` is a finite theorem, not evidence
silently promoted to the universal claim.

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

## 4. The remaining obstruction telescopes atomically

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

The sharp remaining universal question is:

> **OPEN atomic orientation inequality.**  Is
>
> ```text
> D_i(b,c)>=0                         for every 0<=i<b<c? (32)
> ```

If `(32)` holds, then `(29)` gives `D(a,b,c)>=0` for every three-slot
support.  Equations `(25)--(26)` then give `I2<0`, and `(17)` proves

```text
L(H)=L(H^2)=L(H^3)=0                 implies H=0       (33)
```

for arbitrary three-slot factorial support.

The implication is one-way: `(32)` is a sufficient atomic
total-positivity statement, not claimed necessary for `(33)`.  Direct
factorial substitution gives the exact equality family

```text
D_(b-1)(b,b+1)=0.                                      (34)
```

The open issue is proving nonnegativity, and classifying equality, beyond
the finite range in Section 6.

## 5. A continuous Chebyshev-orientation sidecar

The adjacent difference pair already has a strict pointwise orientation.
Put

```text
p=b-a,       r=c-b,
A=b!/a!,     B=c!/b!,

X=s^p/A,     Y=s^r/B.                                  (35)
```

Then

```text
U=f_a(X-1),                  V=f_a X(Y-1),
```

and direct differentiation gives

```text
U V'-U' V
 =f_a^2 X/s [p+rXY-(p+r)Y].                            (36)
```

Every factor in `B` is strictly larger than every factor in `A`, so

```text
B^(1/r)>A^(1/p),                    X^r>Y^p.            (37)
```

Weighted AM--GM now gives

```text
p+rXY
 >=(p+r)(XY)^(r/(p+r))
 >(p+r)Y.                                             (38)
```

Thus

```text
U(s)V'(s)-U'(s)V(s)>0                  for every s>0.   (39)
```

This is a genuine extended-Chebyshev orientation and a plausible source for
`(32)`.  It does **not** by itself prove `(32)`: the atomic determinant
mixes the signed functions `d_i V` and `d_i V^2` against different global
moments.  The missing lemma must transfer the pointwise orientation through
that signed factorial integral.

## 6. Finite exact theorem through exponent 50

The exact companion exhausts all

```text
0<=i<b<c<=50.                                          (40)
```

There are

```text
binom(51,3)=20,825
```

such atomic triples.  It finds

```text
D_i(b,c)>=0                         in all 20,825 cases, (41)
```

with exactly `49` zeros:

```text
(i,b,c)=(j,j+1,j+2),                  0<=j<=48.         (42)
```

The smallest positive value in the normalized factorial basis is

```text
D_0(2,3)=288.                                           (43)
```

For every `0<=a<b<c<=50`, the companion separately verifies the Gram
determinant, `(19)--(25)`, the telescoping identity `(29)`, and

```text
I1<0,                 I2<0.                            (44)
```

Therefore `(33)` is **FINITE-EXACT** for every three-slot support with
top exponent at most `50`.  This conclusion uses integer multinomial
moments only; there is no floating-point sign decision.

## 7. Conditional two-charge Gaussian consequence

Assume `a>=1` and write

```text
h=H/s,
P=W+Z h(ZW),                         W=conj(Z),         (45)
```

for a standard complex Gaussian `Z`.  Charge balance gives, for
`j=1,2,3`,

```text
E[P^(2j)]=binom(2j,j)L(H^j),          E[P^(2j-1)]=0.    (46)
```

Hence the universal atomic inequality `(32)` would imply Gaussian detection
by moment six for every arbitrary three-slot two-charge envelope `(45)`.
The exact finite theorem proves this whenever `c<=50`.

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
| continuous sidecar | strict Wronskian `(39)` |
| finite theorem | all `c<=50`; `20,825` exact atomic checks |
| sharp observed equality | consecutive atom `(j,j+1,j+2)` |
| first missing datum | universal proof of `(32)` |

This theorem does not prove the arbitrary three-slot Strong Factorial
Conjecture, general HYP-8765, or a new proof of all GMC2.  It replaces that
three-slot problem by the explicit signed total-positivity inequality `(32)`
and proves the entire bounded exponent-`50` bank.

## 9. Exact companion

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
the Gram, divisibility, strict-ratio, Wronskian, telescoping, and atomic
identities on all `20,825` triples through `c=50`.  It has explicit
exception gates, no truth-bearing Python assertions, no floating point, and
no scratch dependency.

**Awaiting independent hostile audit; not QED for the universal atomic
inequality.**
