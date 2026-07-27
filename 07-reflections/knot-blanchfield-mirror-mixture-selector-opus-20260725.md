---
source: codex-2026-07-25-knot-blanchfield-selector
status: >
  PROVED + CITED APPLICATION + INDEPENDENTLY HOSTILE-AUDITED;
  SUPERSEDED / GENERALIZED BY
  `knot-torus-mirror-selector-stable-plane-opus-20260725`. A
  cross-frequency inertia functional is
  an explicit additive Gordian-1-Lipschitz invariant. For K=T(2,7), a
  mirror-switched pair has coefficient matrix [[3,-1],[-1,3]], so every
  P,Q mirror mixture and every common-context displacement is at least
  P+Q+2|P-Q|. This equals the real Blanchfield dimension. The stable
  mixture defect is at most 4min(P,Q), sharpening THM-2308's projective
  profile to delta*z<=h(z)<=min(4z,delta(1+z)/2). If delta=4 the entire
  stable two-plane is forced. This does not compute the balanced defect
  delta or an ordinary finite connected-sum unknotting number.
depends_on:
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2308-mirror-double-nakanishi-floor-and-sharp-stable-mixture-profile
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2292-common-catalytic-section-and-helly-calibration-nerve
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
  - knot-torus-mirror-selector-stable-plane-opus-20260725
  - knot-torus-mirror-curvature-extreme-atlas-opus-20260725
  - knot-curvature-mixed-cohabitation-kernel-opus-20260725
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088."
  - "Maciej Borodzik and Stefan Friedl, The unknotting number and classical invariants I, arXiv:1203.3225."
  - "Maciej Borodzik and Stefan Friedl, The unknotting number and classical invariants II, arXiv:1207.2413."
  - "Charles Livingston, Signature invariants related to the unknotting number, arXiv:1710.10530."
---

# A Blanchfield selector sees mirror mixtures that concordance erases

> **GENERALIZED FAMILY THEOREM.** The later torus--mirror packet proves
> the complete two-frequency selector polygon and both stable endpoint
> rigidities for every `T(2,2g+1)`, `g>=2`.  This file remains the
> independently audited `T(2,7)` seed and its Brittenham--Hermiller
> family-specific calibration.

## 1. The primary theorem

For a knot `J`, let `sigma_J(z)` denote its Levine--Tristram signature on
`S^1`, and let `eta_J(z)` denote the associated nullity on `C*`, with the
Borodzik--Friedl convention

```text
eta_J(1)=0.
```

Define

```text
mu_BF(J)
 =(1/2)[max_(z in S^1)(eta_J(z)+sigma_J(z))
       +max_(z in S^1)(eta_J(z)-sigma_J(z))],

eta_BF(J)=max_(z in C*) eta_J(z).                   (1)
```

Borodzik--Friedl prove

```text
n_R(J)=max(mu_BF(J),eta_BF(J)),                     (2)
```

where `n_R` is the minimal real Blanchfield matrix size; their earlier
Blanchfield theorem gives

```text
n_R(J)<=u(J).                                       (3)
```

The two separate maxima in (1) are load-bearing. Replacing them by one
pointwise absolute-signature maximum loses the improvement below.

## 2. Exact `T(2,7)` signature table

Let

```text
K=T(2,7).
```

Use the standard `A_6` Seifert matrix

```text
V_(i,i)=1,

V_(i,i+1)=-1,

all other entries zero.                             (4)
```

The tridiagonal determinant recurrence gives

```text
det(tV-V^T)
 =t^6-t^5+t^4-t^3+t^2-t+1
 =Phi_14(t).                                        (5)
```

For `z=exp(i theta)` and `s=sin(theta/2)`, a diagonal unitary
conjugation turns the Levine--Tristram Hermitian matrix into the real
tridiagonal matrix with diagonal `4s^2` and off-diagonal `-2s`. Its six
eigenvalues are

```text
4s[s-cos(j*pi/7)],                 1<=j<=6.         (6)
```

At the six roots

```text
z=exp(i*pi*r/7),            r in {1,3,5,9,11,13},
```

equation (6) gives, up to the global chirality sign,

```text
sigma_K=(1,3,5,5,3,1),

eta_K=1.                                           (7)
```

On the intervening open arcs,

```text
sigma_K=0,2,4,6,4,2,0,

eta_K=0.                                           (8)
```

All roots of `Phi_14` are simple, so

```text
eta_BF(K)=1.                                       (9)
```

The mirror reverses every signature and preserves every nullity.

## 3. Exact mirror-mixture lower bound

For integers `P,Q>=0`, put

```text
J_(P,Q)=(#^P K)#(#^Q mirror(K)),

N=P+Q,                         d=|P-Q|.             (10)
```

Connected-sum block additivity gives, at every root in (7),

```text
eta_(J_(P,Q))=N,

sigma_(J_(P,Q))=(P-Q)sigma_K.                      (11)
```

With the chirality in (4) fixed, swap `K` and `mirror(K)` if necessary so
that `P>=Q`. Equations (7)--(8) give

```text
max(eta+sigma)=N+5d,

max(eta-sigma)=N-d.                                (12)
```

The off-root value `6d` never exceeds `N+5d`, since `N>=d`.
Therefore

```text
mu_BF(J_(P,Q))=N+2d,

eta_BF(J_(P,Q))=N,

n_R(J_(P,Q))=N+2d.                                 (13)
```

Combining (3) and (13) gives the exact classical lower certificate

```text
u(J_(P,Q))
 >=P+Q+2|P-Q|
 =3max(P,Q)-min(P,Q).                              (14)
```

This is a lower bound, not an equality claim for ordinary unknotting
number.

## 4. Cross-frequency inertia is an explicit Gordian gauge

For arbitrary fixed `alpha,beta in S^1\{1}`, define

```text
Phi_(alpha,beta)(J)
 =(1/2)[eta_J(alpha)+sigma_J(alpha)
       +eta_J(beta)-sigma_J(beta)].                 (15)
```

This is integer-valued and additive. More importantly, it is
**Gordian-1-Lipschitz**. Indeed, two knots differing by one crossing
change admit common even-size Seifert matrices for which every
Levine--Tristram Hermitian form changes by a rank-one semidefinite
matrix, with the same sign at every frequency. For a positive update,

```text
(eta+sigma)/2 = n/2-n_-       changes by 0 or +1,

(eta-sigma)/2 = n/2-n_+       changes by 0 or -1.   (16)
```

The two changes in (16) therefore sum to `-1`, `0`, or `1`; a negative
update reverses the signs. Consequently

```text
|Phi_(alpha,beta)(J)-Phi_(alpha,beta)(J')|
 <=d_G(J,J').                                       (17)
```

This direct inertia proof is stronger than merely observing that one
particular value lies below the Borodzik--Friedl maximum.

## 5. The faithful binary object is a switched gauge pair

Let

```text
zeta_r=exp(i*pi*r/7).
```

The two additive functionals

```text
L(J)
 =(1/2)[eta_J(zeta_5)+sigma_J(zeta_5)
       +eta_J(zeta_1)-sigma_J(zeta_1)],

Lbar(J)=L(mirror(J))                                (18)
```

obey

```text
|L(J)|<=u(J),              |Lbar(J)|<=u(J),

L(K)=3,                  L(mirror(K))=-1,

Lbar(K)=-1,              Lbar(mirror(K))=3.         (19)
```

Hence

```text
max(|L|,|Lbar|)(J_(P,Q))
 =3max(P,Q)-min(P,Q).                               (20)
```

Because the gauges are additive, (17) also gives, for **every** context
`C`,

```text
d_G(J_(P,Q)#C,C)
 >=max(|3P-Q|,|-P+3Q|)
 =3max(P,Q)-min(P,Q).                               (21)
```

Thus the same bound survives catalytic localization and connected-sum
homogenization; a common context cannot erase it. It is exactly (14).

This is not naturally a tournament. The faithful carrier is a
mirror-switched ordered pair of signed spectral gauges, with root nullity
as the mirror-even coordinate and signature as the mirror-odd coordinate.
Since `mirror(K)` represents `-[K]` in concordance, any concordance
homomorphism kills `K#mirror(K)` and sees only `P-Q`; it cannot recover
the `P+Q` nullity term in (13). The Blanchfield sidecar is therefore
genuinely nonconcordance information.

After homogenization, (15) says that the stable dual ball on
`span{K,mirror(K)}` contains the two mixed supports

```text
(3,-1),                  (-1,3).                    (22)
```

Their coefficient matrix has determinant `8`. They are the
operation-ready replacement for a cosmetic binary relation and already
prove stable rank two on this mirror lattice without a common-fibre or
Yang-type sidecar.

## 6. Stable defect and curvature budget

Use THM-2308's notation

```text
a=[K],                  b=[mirror(K)],

p=u_hash,

D(P,Q)=3(P+Q)-p(Pa+Qb),

h(z)=D(1,z),            delta=h(1),       1<=delta<=4.     (23)
```

The lower endpoint is the explicit Brittenham--Hermiller construction

```text
u(K#mirror(K))<=5,
```

iterated under connected sum, while the upper endpoint is the new
Blanchfield floor `p(a+b)>=2` from (13). Thus the interval
`1<=delta<=4` records the exact present gap between a diagrammatic upper
route and the strongest classical lower certificate; neither endpoint is
being assumed.

Apply (14) to every common multiple and divide. For `P>=Q`,

```text
p(Pa+Qb)>=3P-Q,

D(P,Q)<=4Q,

h(z)<=4z.                                          (24)
```

Together with THM-2308's existing calibration and mirror envelope,

```text
delta z
 <=h(z)
 <=min(4z,delta(1+z)/2).                            (25)
```

Let the canonical curvature measure on `(0,1)` be

```text
kappa=-D^2h,

A=integral_(0,1)(1-t)dkappa(t),

B=integral_(0,1)t dkappa(t).                        (26)
```

The representation

```text
h(z)=delta z+integral G(z,t)dkappa(t)
```

turns (25) into the exact moment budgets

```text
A<=4-delta,

B<=delta/2.                                        (27)
```

The old THM-2308-only budget `A<=6-delta` remains an abstract statement
about the information available before (14), but it is no longer sharp
for the actual knot.

## 7. Sharp new stable-norm envelope

Put

```text
e=a+b,                  o=a-b,

X=|s|,                  Y=|t|.
```

Equations (24)--(25) give the global envelope

```text
max((6-delta)X, 2X+4Y, 6Y)
 <=p(se+to)
 <=(6-delta)max(X,Y)+delta Y.                       (28)
```

For `1<=delta<4`, both sides are attained by admissible abstract stable
norms with all current knot data. In curvature coordinates:

```text
lower-defect / upper-norm extremal:
  kappa=0,
  h(z)=delta z;

upper-defect / lower-norm extremal:
  t_0=delta/(8-delta),
  m_0=(8-delta)/2,
  kappa=m_0 delta_(t_0),
  h(z)=min(4z,delta(1+z)/2).                         (29)
```

The second measure saturates both budgets:

```text
m_0(1-t_0)=4-delta,

m_0t_0=delta/2.                                    (30)
```

Thus the curvature atlas remains infinite-dimensional for every
`delta<4`, but with a strictly smaller feasible moment polytope.

## 8. Rigidity at the endpoint `delta=4`

At `delta=4`, equation (27) gives

```text
A<=0.
```

Since `kappa` is nonnegative and `1-t>0` on the canonical open interval,

```text
kappa=0,

h(z)=4z.                                           (31)
```

The two norm envelopes in (28) coincide:

```text
p(se+to)
 =2max(|s|,|t|)+4|t|.                              (32)
```

The formal limit of the second atom in (29) is `2delta_1`, but endpoint
mass is invisible to the Green kernel and is not part of the canonical
measure. Equation (31), not an endpoint atom, is the faithful statement.

Therefore the earlier THM-2308-only claim of continuum ambiguity at every
`delta in [1,4]` must be restricted, after the Blanchfield input, to

```text
1<=delta<4.                                        (33)
```

This does not prove `delta=4`; it proves that **if** the balanced stable
defect has that endpoint value, the entire stable plane is already
determined.

## 9. General torus-mirror family

The same calculation works for

```text
K_g=T(2,2g+1).
```

The standard `A_(2g)` Seifert matrix has simple unit-root signatures

```text
1,3,...,2g-1,...,3,1
```

and root nullity one. Hence

```text
n_R((#^P K_g)#(#^Q mirror(K_g)))
 =P+Q+(g-1)|P-Q|,                                  (34)

u>=P+Q+(g-1)|P-Q|.
```

Relative to the separate cost `g(P+Q)`, the stable mixture defect is at
most

```text
2(g-1)min(P,Q).                                    (35)
```

The `g=3` case is exactly (24). This family identifies the mechanism:
mirror-even root nullity and mirror-odd signature interact through the two
separate maxima in (1).

## 10. Scope

The theorem supplies a new raw lower certificate, a sharper stable profile
envelope, and conditional rigidity at `delta=4`. It does not:

```text
compute delta;

prove ordinary or stable unknotting additivity;

produce a positive Gordian catalyst;

classify knots;

or determine an interior profile when delta<4.                      (36)
```

## 11. Exact audit

```text
C:\Users\Eliott\.cache\codex-runtimes\codex-primary-runtime\
dependencies\python\python.exe
  04-computation/knot_t27_blanchfield_gordian_selector_probe.py
```

The audit reconstructs `Phi_14` by the tridiagonal determinant
recurrence, verifies the six root and seven off-root inertia cells by the
exact monotonic comparison `sign(2j-|7-r|)`, checks the
Borodzik--Friedl mixture formula and the determinant-eight gauge on all
`0<=P,Q<=8`, and checks the curvature extremal budgets exactly.

Stored transcript:

```text
05-knowledge/results/knot_t27_blanchfield_gordian_selector_probe.out
```

Normal, optimized (`-O`), and stored transcripts match. LF-byte SHA-256:

```text
script  87df1cf64b7c3772565f0a99f2dcbf97d279c88d95e2084025043c92aadbcd28
output  de36da1ccbeb21d697273f452217d9b7ada596f81d3bdd2b20128cd1d5d2cb68
```
