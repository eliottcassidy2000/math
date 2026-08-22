---
id: THM-3667
title: "LRC optimal positive three-twist frame and rigidity tradeoff"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  Among all normalized positive three-site averaging masks on F13^2, the
  worst nontrivial Fourier multiplier is maximized by explicit real
  cyclotomic weights.  The exact lower and upper frame constants and all
  equality frequencies are determined.  At the unique optimizer (up to
  swapping the positive sites), exactly thirteen eigenvalue pairs collide,
  enlarging the full linear centralizer from dimension 169 to 195.  Equal
  weights retain simple spectrum but have a smaller gap.  This optimizes an
  abstract target-twist observer; it does not prove covering-row
  nonconstancy or LRC(14).
source: kps-s192 / weighted continuation of THM-3665, 2026-08-21
depends_on:
  - THM-3665-lrc-support-minimal-three-twist-target-detector
related:
  - THM-3666-lrc-owner-pivot-dual-pair-swap-twist-basis
script: 04-computation/lrc_optimal_positive_three_twist_frame_thm3667.py
output: 05-knowledge/results/lrc_optimal_positive_three_twist_frame_thm3667.out
script_sha256: f90565ddce543be2a8a37fe715e32a2eb1d660380b51315ba3847fd3b5889dce
output_sha256: 8c70caa5ad5fb990c4078e1af3c072471d6e2352cd9800ed65000f6a72642a88
semantic_sha256: 458ca3b1a351bb98be8968f79fc9e7c54b3bb571d654eb80cdd26d6842fca53c
hash_basis: raw LF bytes
---

# THM-3667 -- optimal conditioning forces thirteen spectral collisions

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  THM-3665 used equal positive weights because they give a simple
spectrum and an especially clean maximum principle.  They are not the
best-conditioned positive weights.  The exact optimizer exposes a sharp
conditioning-versus-rigidity boundary.

## 1. Positive three-site family

On `F13^2`, fix an ordered basis `(e_1,e_2)` and let

```text
m_(a,b)=a delta_0+b delta_(e_1)-delta_(e_2),
a>0, b>0, a+b=1.                                    (1)
```

Every mask in (1) has mean zero and Fourier multiplier

```text
lambda_(a,b)(u,v)=a+b zeta^(-u)-zeta^(-v),
zeta=exp(2*pi*i/13).                                 (2)
```

The triangle-equality argument of THM-3665 applies with unequal positive
weights: a zero in (2) would make a convex combination of two unit complex
numbers have modulus one.  Both phases must coincide with the third, so
`u=v=0`.  Thus every member of (1) is a support-minimal detector.

Normalize its worst multiplier by

```text
mu(a)=min_((u,v)!=(0,0)) |lambda_(a,1-a)(u,v)|.     (3)
```

Swapping the two positive support sites by an affine rebasing gives
`mu(a)=mu(1-a)`, so it suffices to optimize on `a>=1/2`.

## 2. Two hostile modes force the optimizer

Put

```text
s=sin(pi/13), c=cos(pi/13).                         (4)
```

The modes `(1,0)` and `(2,1)` give

```text
f_1(a)=2(1-a)s,

f_2(a)=2s sqrt(s^2+(2a-1)^2 c^2).                  (5)
```

On `a>=1/2`, `f_1` is strictly decreasing and `f_2` is strictly increasing.
Therefore

```text
mu(a)<=min(f_1(a),f_2(a))                           (6)
```

is globally maximized only where the two values meet.  Solving gives

```text
a_*=(4c^2-2)/(4c^2-1),
b_*=1/(4c^2-1).                                     (7)
```

Equivalently, with

```text
y=2cos(2*pi/13),                                    (8)
```

the weights are

```text
a_*=y/(y+1),       b_*=1/(y+1).                     (9)
```

Both are positive.  Numerically they are approximately

```text
(a_*,b_*)=(0.6391079971,0.3608920029).              (10)
```

## 3. Exact all-frequency lower certificate

Let

```text
x=2cos(pi/13),
P(x)=x^6-x^5-5x^4+4x^3+6x^2-3x-1=0,
D=x^2-1, A=x^2-2.                                   (11)
```

Thus `(a_*,b_*)=(A/D,1/D)`.  Define

```text
C_n(x)=2cos(n*pi/13),
R_k=2-C_(2k)(x)=4sin(k*pi/13)^2.                    (12)
```

An exact expansion of (2) gives

```text
D^2 |lambda_(a_*,b_*)(u,v)|^2
 =A D R_v+D R_(u-v)-A R_u.                          (13)
```

Reducing (13) modulo `P(x)`, then evaluating at the isolated root

```text
1941883/10^6 < x < 1941884/10^6,                   (14)
```

proves

```text
D^2 |lambda(u,v)|^2 >= R_1                         (15)
```

for every nontrivial frequency.  Equality occurs exactly at

```text
(1,0),(12,0),(2,1),(11,12).                         (16)
```

The companion proves all 164 strict signs by rational interval arithmetic;
the four equalities reduce to the zero polynomial in `Q[x]/(P)`.  Consequently
the upper bound forced in Section 2 is attained, and

```text
max_(0<a<1) mu(a)^2
 =(4-x^2)/(x^2-1)^2.                                (17)
```

The optimizers are exactly `a_*` and `b_*`, interchanged by swapping the two
positive sites.

For comparison, equal normalized weights have gap

```text
mu(1/2)=2sin(pi/13)^2 approximately 0.114544,
mu(a_*) approximately 0.172729.                     (18)
```

Thus optimal weighting improves the amplitude gap by about 51 percent.

## 4. Exact upper frame constant

The same exact scan gives

```text
max_(u,v)|lambda_(a_*,b_*)(u,v)|
 =2cos(pi/26),                                      (19)
```

with equality exactly at `(0,6),(0,7)`.  Hence, for every twist profile `H`,

```text
[(4-x^2)/(x^2-1)^2] ||H-H_bar||_2^2
 <=||m_(a_*,b_*)*H||_2^2
 <=4cos(pi/26)^2 ||H-H_bar||_2^2.                  (20)
```

Both constants are sharp character profiles.

## 5. The optimizer loses simple spectrum

Multiplying (2) by `D=y+1` gives the cyclotomic integer

```text
D lambda(u,v)=y+zeta^(-u)-(y+1)zeta^(-v).           (21)
```

There are exactly 156 distinct values in the 169-frequency spectrum:

```text
143 singleton eigenvalues,
13 double eigenvalues.                              (22)
```

The double classes have the uniform description

```text
{(u,u+1),(u+3,u+2)},             u in F13,          (23)
```

with coordinates reduced modulo thirteen.  Their identity follows already
from

```text
b_*(1-zeta^(-3))=zeta^(-1)(1-zeta^(-1)),            (24)
```

because `b_*=1/(1+zeta+zeta^(-1))`.  Exact reduction in the basis
`1,zeta,...,zeta^11` proves that (23) is the complete collision list.  The
zero eigenvalue remains a singleton at `(0,0)`, so detector invertibility on
mean-zero profiles is unaffected.

For a diagonalizable operator, the full linear centralizer dimension is the
sum of squared eigenvalue multiplicities.  Thus

```text
dim Cent_End(C_m)=143+13*2^2=195.                   (25)
```

By contrast, THM-3665's equal-weight mask has 169 distinct eigenvalues and
centralizer dimension 169.  Exact maximin conditioning therefore forces a
strict loss of linear rigidity.

The loss is confined to the exact optimizer.  Pairwise collisions impose
only finitely many exceptional algebraic weight ratios, so generic weights
arbitrarily close to `(a_*,b_*)` recover simple spectrum with arbitrarily
near-optimal gap.

## 6. Interpretation for the LRC target plane

In THM-3666's pair-swap coordinates, the optimal local test is

```text
a_* H(x,y)+b_* H(x-1,y)-H(x,y-1).                  (26)
```

Its vanishing is a positive weighted harmonic recurrence.  The finite
maximum-modulus principle still makes global vanishing equivalent to a
constant twist profile, while (20) gives the best possible uniform stability
among normalized positive three-site tests.

There are now two useful observables rather than one universally superior
choice:

```text
equal rational weights:
  simpler arithmetic, simple spectrum, strongest centralizer rigidity;

optimal cyclotomic weights:
  largest frame gap, but thirteen spectral collisions.             (27)
```

Neither observable proves that the canonical marked current is nonconstant
on a hypothetical covering row.  The theorem concerns the unrestricted
mod-thirteen target aggregate, not the all-`91`-unit projector, and does not
close LRC(14).

## 7. Exact companion

Reproduce with

```bash
python3 -B 04-computation/lrc_optimal_positive_three_twist_frame_thm3667.py
python3 -B -O 04-computation/lrc_optimal_positive_three_twist_frame_thm3667.py
```

Both streams match the stored transcript.  The assertion-free companion
proves the isolating interval from rational pi bounds and alternating cosine
bounds, performs every sign comparison in the exact degree-six real
cyclotomic field, classifies both equality sets, and reduces all 169
eigenvalues in `Q(zeta_13)` to obtain the complete collision census.  **QED.**
