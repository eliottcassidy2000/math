---
id: THM-2573
title: "Logarithmic Abel normal and common-endpoint jump pairing"
status: >
  PROVED + VERIFIED-EXACT.  Let L_s,R_s be finite complex step layers on
  the circle with L_s conjugate(R_s)=0 almost everywhere, and fix a physical
  frequency offset M.  Their ordinary full-X endpoint sum is zero.  If the
  left layer alone is Poisson weighted by rho^|X|, its first potentially nonzero boundary
  scale is logarithmic: after division by
  (1-rho)log(1/(1-rho)), the limit is -1/(2 pi^2) times the Mth Fourier
  moment of the aggregate common-jump pairing
  sum_x Delta L_s(x) conjugate(Delta R_s(x)) delta_x.  Weighting both
  endpoints with the same rho doubles the normal; the general factor is the
  total Abel speed.  For nonnegative disjoint layers the negative jump
  pairing is a positive handoff measure.  A finite common-target DFT commutes
  with the normal extraction, and some nonzero target colour fires for some
  integer M exactly when these boundary-pair measures are not all equal.
  But the LRC deep bank samples only M=m c_3 with gcd(m,91)=1, which has a
  sharp alternating-grid blind kernel.  On the THM-2569 packet, freezing the
  target-informed selector still makes any resulting colour auxiliary;
  lawfulness requires a covariant whole-layer orbit before smoothing.
  Coincident factor endpoints must be merged into total-layer one-sided
  jumps.  The normal survives only as a singular Abel boundary coefficient,
  not as the ordinary full-X pushforward or one fixed-X current.  No
  THM-2334 relation current, row exclusion, or LRC(14) conclusion follows.
source: common-endpoint-seam-2026-07-27-abel-normal
depends_on: []
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2448-right-endpoint-cospan-transition-atlas
  - THM-2568-full-x-transition-annihilation-and-refined-pair-drift-boundary
  - THM-2569-stationary-diagonal-conditioned-paired-corner-and-frozen-future-role-boundary
  - THM-2574-oriented-tooth-component-holonomy-and-fixed-frequency-descent
script: 04-computation/lrc14_logarithmic_abel_normal_thm2573.py
output: 05-knowledge/results/lrc14_logarithmic_abel_normal_thm2573.out
script_sha256: 27406e4b4490f0a799ea08910d924b95919055c4a42c18aa37ec5e8f321bb749
output_sha256: d61b1a50d001d5f23536e32554aa77872c92b322eb99fd6fd6488b752d69f75f
hash_basis: working-tree bytes (LF)
---

# THM-2573 -- the first boundary term after full-X annihilation

**PROVED + VERIFIED-EXACT.**

THM-2568 proves that a complete danger-to-safe endpoint transition vanishes
after the ordinary physical-frequency sum.  That is the value of the
endpoint pairing at the boundary of the Poisson disk.  It does not say that
the approach to the boundary is flat.

For finite step layers the first potentially nonzero universal scale is not an ordinary derivative:
the derivative diverges logarithmically.  The correctly renormalized normal
is an explicit codimension-one object:

```text
ordinary full-X value                    = pointwise overlap = 0;

logarithmic Abel normal                  = common-boundary jump pairing.
                                                               (1)
```

This is the physical normal/jet option left open by THM-2568.  It gives a
precise lawful target observable once a lawful whole-layer orbit is supplied,
and it identifies exactly why the present THM-2569 packet does not yet supply
one.

## 1. Whole-layer Abel normal

Use the Fourier convention

```text
f_hat(X)=integral_T f(x)exp(-2 pi i Xx)dx.                    (2)
```

Fix an integer physical offset `M`.  For each target label `s`, let
`L_s,R_s` be complex-valued step functions on the circle with finitely many
jumps and

```text
L_s(x)conjugate(R_s(x))=0                 almost everywhere. (3)
```

At a discontinuity `x`, write the cyclic one-sided **total-layer** jump

```text
Delta L_s(x)=L_s(x+)-L_s(x-),

Delta R_s(x)=R_s(x+)-R_s(x-).                         (4)
```

Changing open to closed interval conventions changes point values only and
does not change (2) or (4).  If several factor boundaries coincide, first
form the whole products `L_s,R_s`, take their two one-sided traces, and only
then use (4).

For `0<rho<1`, smooth the left whole layer and define

```text
J_rho(s;M)
 =sum_(X in Z)rho^|X|
    L_s_hat(X)conjugate(R_s_hat(X+M)).                       (5)
```

The unweighted sum is absolutely convergent by Cauchy--Schwarz.  Endpoint
Parseval and (3) give

```text
J_1(s;M)
 :=sum_X L_s_hat(X)conjugate(R_s_hat(X+M))

 =integral_T L_s(x)conjugate(R_s(x))exp(2 pi i Mx)dx
 =0.                                                         (6)
```

Define the common-jump moment

```text
B_s(M)
 =sum_(x in Jump(L_s) intersect Jump(R_s))
    Delta L_s(x)conjugate(Delta R_s(x))exp(2 pi i Mx).       (7)
```

Then the logarithmic Abel normal exists and is

```text
N_s(M)
 :=lim_(rho->1-)
     J_rho(s;M)/[(1-rho)log(1/(1-rho))]

 =-B_s(M)/(2 pi^2).                                         (8)
```

More generally, for fixed `alpha,beta>=0`, `alpha+beta>0`, replacing the
weight in (5) by

```text
rho^(alpha|X|+beta|X+M|)                                   (9)
```

multiplies the right side of (8) by `alpha+beta`.  Thus smoothing the right
layer alone gives the same canonical normal; smoothing both layers with the
same `rho` doubles it.  Equivalently, smoothing both with `sqrt(rho)` gives
exactly (8).  The Abel-speed invoice is part of the theorem, not a convention
hidden in the constant.

## 2. Proof of the jump formula

For a step function `F`, distributional integration by parts gives, when
`X!=0`,

```text
F_hat(X)
 =1/(2 pi i X) sum_x Delta F(x)exp(-2 pi i Xx).              (10)
```

For `X` different from `0,-M`, substitute (10) into one cross term:

```text
L_hat(X)conjugate(R_hat(X+M))
 =1/[4 pi^2 X(X+M)]
   sum_(x,y)Delta L(x)conjugate(Delta R(y))
      exp(2 pi i X(y-x))exp(2 pi i My).                      (11)
```

The finitely many exceptional frequencies contribute `O(1-rho)` to
`J_rho-J_1`.  In (11), a coincident pair `x=y` has no oscillation.  Replacing
`X(X+M)` by `X^2` changes the Abel difference by `O(1-rho)`, and the positive
and negative tails each obey

```text
sum_(n>=1)(rho^n-1)/n^2
 =-(1-rho)log(1/(1-rho))+O(1-rho).                         (12)
```

One proof of (12) differentiates the dilogarithm series:

```text
sum_(n>=1)(rho^n-1)/n^2
 =integral_1^rho [-log(1-t)/t]dt.                            (13)
```

Hence the two tails contribute `-2`, and the factor `1/(4 pi^2)` in
(11) gives `-1/(2 pi^2)`.

When `x!=y` on the circle, the differentiated tail is a bounded Abel sum of

```text
exp(2 pi i n(y-x))/n.                                      (14)
```

Dirichlet summation therefore makes its contribution `O(1-rho)`.  The
`O(1/|X|^3)` denominator correction is absolutely `O(1-rho)`.  Only the
coincident pairs survive division by the logarithmic scale, proving (8).
The same argument with (9) replaces the tail exponent by
`(alpha+beta)|X|+O(1)` and proves the general factor.

This proof handles arbitrary complex jump amplitudes.  The conjugation,
the phase `exp(2 pi i Mx)`, the exceptional shift `X=-M`, and both Fourier
tails are all fixed by (10)--(11); no reality assumption is used.

## 3. Positive handoff measure for Boolean endpoint layers

Suppose now that `L_s,R_s` are real, nonnegative, and disjoint.  At one
common discontinuity put

```text
(a,b)=(L_s(x-),L_s(x+)),

(c,d)=(R_s(x-),R_s(x+)).                                  (15)
```

Disjointness gives `ac=bd=0`.  A direct four-case check gives

```text
(b-a)(d-c)<=0.                                             (16)
```

Define the finite positive handoff measure

```text
nu_s
 =-sum_x Delta L_s(x)Delta R_s(x) delta_x.                  (17)
```

Equation (8) becomes

```text
N_s(M)
 =1/(2 pi^2) integral_T exp(2 pi i Mx)dnu_s(x).             (18)
```

In particular `N_s(0)>=0`, with strict inequality exactly when at least one
total-layer boundary passes directly between positive left and positive
right traces.  Disjoint positive supports separated by a gap have
`nu_s=0`; positivity of both layers alone does not force a normal.

The adjective **total-layer** is essential.  A common filter may die at the
same point at which a danger gate meets its complement.  Then one of the
whole layers has zero jump and that apparent factorwise handoff contributes
nothing to (17).  This is the boundary analogue of forming Boolean products
before Poisson smoothing in THM-2452.

## 4. Common-target DFT and the exact drift dichotomy

Put `p=13`, `zeta=exp(2 pi i/p)`, and suppose the thirteen pairs
`(L_s,R_s)` are obtained from one endpoint cospan by a **common lawful target
action**.  Define

```text
Jhat_rho(q;M)
 =1/p sum_(s in F_p)J_rho(s;M)zeta^(q s),

Nhat(q;M)
 =1/p sum_(s in F_p)N_s(M)zeta^(q s).                       (19)
```

The finite target sum commutes with the limit in (8), so

```text
lim_(rho->1-)
 Jhat_rho(q;M)/[(1-rho)log(1/(1-rho))]
 =Nhat(q;M).                                                (20)
```

Thus the ordinary full-`X` pushforward still vanishes by (6), in every
target character, while its singular logarithmic normal can be nonzero.
This is a different boundary observable, not a contradiction to THM-2568.

For the complex measures

```text
beta_s
 =sum_x Delta L_s(x)conjugate(Delta R_s(x))delta_x,          (21)
```

the following are equivalent:

```text
Nhat(q;M)=0 for every q!=0 and every M in Z;

beta_s is independent of s.                                (22)
```

Indeed the first line says every integer Fourier moment of `beta_s` is
independent of `s`; Fourier uniqueness for finite atomic measures gives the
second.  Conversely the second makes every nonzero target DFT vanish.  In
the positive case one may replace `beta_s` by `nu_s`.

Equation (22) is the codimension-one functional form of target drift:

```text
moving boundary-pair measure
  or
target-invariant boundary-pair measure.                     (23)
```

It is strictly sharper than comparing total boundary masses.  Spatial
motion can retain `nu_s(T)` while changing its nonzero physical moments.

## 5. The deep `91`-unit bank sees only a restricted projection

The LRC triangle does not permit every physical offset.  Its offsets have
the form

```text
M=m c_3,                   gcd(m,91)=1.                       (24)
```

Let `pi_(c_3):x -> c_3x` on the circle.  Equation (18) samples only the
unit-frequency moments of the pushforward

```text
(pi_(c_3))_*nu_s.                                           (25)
```

Therefore (22) cannot be invoked from the allowed bank alone.  The exact
restricted dichotomy is only

```text
some allowed unit moment varies with s
  -> a nonzero target Abel normal;

all allowed unit moments are constant in s
  -> deep-unit-blind boundary branch.                        (26)
```

The second branch need not make the measures equal.  A sharp positive
hostile is obtained by pulling an alternating fourteen-boundary Boolean set
back through `x -> c_3x`.  The pushforward (25) is a positive multiple of the
uniform measure on a translated fourteen-point grid.  Its Fourier moments
vanish unless `14|m`; every `m` coprime to `91` is therefore blind, although
the translated measures are different.  No unrestricted Fourier-uniqueness
claim is available without an additional physical sidecar.

## 6. Exact typing on the THM-2569 common packet

Resolve the THM-2569 packet by its stationary root `h` and form complete
old-danger and repaired-safe endpoint layers.  They are rational step
functions and have the orthogonality (3), root stratum by root stratum.
Consequently (8) applies to their **actual total layers**, including the
stationary future multiplier.  It answers the analytic question exactly:

```text
ordinary full-X endpoint current                  zero;

logarithmic physical normal                       jump pairing (7). (27)
```

The future/root multiplier is not harmless for the normal.  It changes
one-sided total-layer traces and may kill, create, or weight common jumps.
The coincident-endpoint convention in Section 3 must be applied after this
multiplier is inserted.

There is, however, a separate target-action type gate.  In the current
THM-2569 composition, the target-informed old selector and its unshifted
danger facts remain frozen inside `w_(N,h)` while the moving endpoint factors
are co-shifted.  Formula (8) still computes the resulting analytic normal,
but a nonzero DFT of that frozen family is only an auxiliary colour.  Calling
it a THM-2334 target charge would repeat MISTAKE-266.

For (19)--(20) to be a lawful target character, one must first construct

```text
(L_(s,h),R_(s,h))_(s in F_13)                              (28)
```

by co-shifting **every** target-active factor, including the selector and
the later occurrence, before forming whole layers and before Abel smoothing.
The root equality `h=b` does not construct (28).  THM-2574's oriented tooth
component local system is a candidate reference/trivialization for this
missing lift, but it is complex and gauge-chosen rather than a positive
Boolean carrier.

Even after (28), the normal is not one fixed-`X` current: every fixed term is
lower order after division by the logarithmic scale, and the survivor comes
from the collective high-frequency tail.  Nor is whole-layer smoothing plus
a singular boundary renormalization automatically the factorwise Abel
relation current of THM-2334.  Word, owner, deep-unit admissibility, and the
relation-address pushforward still require a separate compatibility theorem.

## 7. Sharp controls

Four elementary examples separate the statements.

### 7.1 Moving complement: the normal can carry target charge

Let

```text
P_s=1_[s/13,s/13+1/4),

L_s=P_s,                  R_s=1-P_s.                         (29)
```

There are two unit handoff atoms, at `s/13` and `s/13+1/4`.  Hence at
`M=1`,

```text
N_s(1)
 =exp(2 pi i s/13)(1+i)/(2 pi^2).                           (30)
```

With the plus-sign DFT convention in (19), the unique nonzero target colour
is

```text
q=-1=12 mod 13.                                            (31)
```

The ordinary full-`X` current is nevertheless zero for every `s`.

### 7.2 Fixed complement: positive normal need not be charged

If `P_s=P` is independent of `s`, then `N_s(0)>0` whenever `P` has a
boundary, but every `q!=0` coefficient in (19) vanishes.  Boundary handoff is
not by itself target drift.

### 7.3 Separated supports: orthogonality need not have a normal

Take


```text
L=1_[0,1/4),                  R=1_[1/2,3/4).                 (32)
```

The jump sets are disjoint, so `B(M)=0` for every `M` and the logarithmic
normal vanishes.  The approach to the zero full-`X` value is only
`O(1-rho)`.

### 7.4 Coincident banks: factorwise boundary counting is false

Let

```text
P=1_[0,1/2),                  H=1_[0,3/4),

L=PH,                         R=(1-P)H.                       (33)
```

The naked pair `P,1-P` has handoffs at `0` and `1/2`.  The total layers in
(33) have a common jump only at `1/2`; the coincident death of `H` removes
the apparent handoff at `0`.  This is unchanged by strict/open endpoint
choices and proves that aggregate traces in (4) are mandatory.

## 8. Exact companion

Run

```bash
python3 04-computation/lrc14_logarithmic_abel_normal_thm2573.py
python3 -O 04-computation/lrc14_logarithmic_abel_normal_thm2573.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_logarithmic_abel_normal_thm2573.out.
```

The dependency-free exact referee checks:

- the two-tail `-2` normalization and `-1/(2 pi^2)` jump factor;
- complementary, separated-support, and coincident-total-layer jump
  measures on exact cyclic rational grids;
- all `26` moving-interval handoff atoms and the target sign `q=12`;
- the target-neutral positive-normal hostile;
- the uniform fourteen-boundary measure; and
- `1,584` signed `91`-unit samples through `|m|<=1000`, together with the
  symbolic implication `gcd(m,91)=1 => 14 does not divide m`.

The Abel asymptotic, off-diagonal Dirichlet bound, finite-measure uniqueness,
and lawfulness distinction are symbolic proofs above, not numerical
extrapolations. **QED.**
