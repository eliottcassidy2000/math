---
source: codex-2026-07-25-knot-curvature-cohabitation
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On each ordered
  multiplicity chamber of the torus--mirror stable plane, every mixed
  forward difference of the stable connected-sum defect is the integral
  of its curvature measure against one explicit nonnegative clipped-tent
  kernel.  All such chamberwise mixed differences vanish if and only if
  the interior curvature measure vanishes.  This is the stable
  continuum analogue of THM-2348's mixed-type cohabitation/ANOVA sector,
  with two load-bearing caveats: THM-2348 is a finite allocation theorem,
  and the unavoidable delta*min(P,Q) ridge lives on the diagonal and is
  invisible to any one interior chamber difference.  A complete row of
  chamber increments nevertheless telescopes to a clipped curvature
  moment, and that row together with its adjacent diagonal cell recovers
  delta exactly.  No finite unknotting equality, robust rectangularity,
  or knot realization is asserted.
depends_on:
  - knot-torus-mirror-selector-stable-plane-opus-20260725
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2346-global-allocation-anova-normal-form
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
  - knot-torus-mirror-curvature-extreme-atlas-opus-20260725
script: 04-computation/knot_curvature_mixed_cohabitation_kernel_probe.py
output: 05-knowledge/results/knot_curvature_mixed_cohabitation_kernel_probe.out
script_sha256: f0287fe41a65fb9828d67ac5dac1b35147d9cfe98c54f82ba4a0e3ec57da16ba
output_sha256: 8cc015ee0b3214cd0e53d40f94f8344670255096218aa8bdec75cf44006570ca
hash_basis: LF-normalized working-tree bytes
---

# Stable curvature is a mixed cohabitation kernel

## 1. From the projective profile back to multiplicities

Use the audited torus--mirror notation

```text
a=[K_g],                         b=[mirror(K_g)],

D(P,Q)=g(P+Q)-p(Pa+Qb),          P,Q>=0.             (1)
```

The defect is symmetric and homogeneous.  In the chamber `P>=Q`, write

```text
D(P,Q)=P h(Q/P),

h(z)=delta z+integral_(0,1)G(z,t)d kappa(t),         (2)

G(z,t)=min(z,t)-zt.
```

Multiplying the Green kernel by `P` exposes a much simpler lattice
formula:

```text
P G(Q/P,t)=min(Q,Pt)-Qt.                            (3)
```

Hence

```text
D(P,Q)
 =delta Q
  +integral_(0,1)[min(Q,Pt)-Qt]d kappa(t).          (4)
```

The affine terms in (4) will disappear under one mixed difference.  The
entire chamberwise interaction is therefore carried by `kappa`.

## 2. Exact clipped-tent formula

For integers

```text
P>=1,                       0<=Q<=P-1,              (5)
```

define the forward mixed cohabitation increment

```text
C_(P,Q)
 =D(P+1,Q+1)-D(P+1,Q)
  -D(P,Q+1)+D(P,Q).                                 (6)
```

Condition (5) keeps all four lattice points in the closed chamber
`P>=Q`.  Put

```text
clip(s)=
  0,        s<=0;
  s,        0<=s<=1;
  1,        s>=1.                                   (7)
```

Then

```text
C_(P,Q)=integral_(0,1) K_(P,Q)(t)d kappa(t),        (8)

K_(P,Q)(t)
 =clip((P+1)t-Q)-clip(Pt-Q).                        (9)
```

To prove this, observe for fixed `p,t` that

```text
min(Q+1,pt)-min(Q,pt)=clip(pt-Q).                   (10)
```

Apply (10) at `p=P+1` and `p=P` to the four `min` terms in
(6).  The terms `delta Q` and `-Qt` have zero mixed difference, giving
(8)--(9).

Since `t>0` and `clip` is nondecreasing,

```text
K_(P,Q)(t)>=0,

C_(P,Q)>=0.                                         (11)
```

This is a sign theorem for the stable defect which is not visible from
an arbitrary symmetric norm alone; it comes from the positive curvature
representation.

## 3. The lattice family detects every interior curvature atom

The kernels in (9) do not merely give necessary inequalities.  They
detect whether the curvature measure is zero.

Fix `0<t<1`.  For any positive integer `P`, choose

```text
Q=floor(Pt).                                        (12)
```

Then `0<=Q<=P-1` and, with `s=Pt-Q in [0,1)`,

```text
K_(P,Q)(t)
 =clip(s+t)-s
 =min(t,1-s)>0.                                     (13)
```

Thus the countable family of open positivity sets
`{K_(P,Q)>0}` covers `(0,1)`.  Since both `kappa` and every kernel are
nonnegative,

```text
C_(P,Q)=0 for every (P,Q) in (5)

iff

kappa=0.                                            (14)
```

The forward implication follows because every integral in (8) being
zero makes `kappa` vanish on its kernel's positivity set; countable
subadditivity and (13) then kill all of `(0,1)`.  The reverse implication
is immediate.

For an extreme one-kink profile `kappa=m delta_tau`,

```text
C_(P,Q)=m K_(P,Q)(tau),                             (15)
```

and for a two-kink extreme it is the corresponding positive sum of two
clipped tents.  The atomic extreme atlas therefore becomes an exact
finite-difference sampling model.

## 4. Complete-row telescoping and exact ridge recovery

The kernels have a second exact feature that a single chamber cell
conceals.  For `0<=t<=1`,

```text
sum_(Q=0)^(P-1) K_(P,Q)(t)
 =min(t,P(1-t)).                                   (16)
```

Indeed, for every `x>=0`,

```text
sum_(Q=0)^(n-1) clip(x-Q)=min(x,n).
```

Apply this once with `x=(P+1)t,n=P` and once with `x=Pt,n=P`.
Consequently the complete chamber row

```text
S_P=sum_(Q=0)^(P-1) C_(P,Q)
    =integral_(0,1) min(t,P(1-t))d kappa(t).       (17)
```

The row moments increase monotonically to the first curvature moment:

```text
S_P increases to integral_(0,1)t d kappa(t),       (18)

integral t d kappa-S_P
 =integral_[P/(P+1),1] ((P+1)t-P)d kappa(t).
```

Thus an atomic profile whose support lies below `P/(P+1)` has already
stabilized exactly at row `P`; no limiting argument is then needed.

Now retain the mixed cell crossing the symmetry diagonal,

```text
R_P
 =D(P+1,P+1)-D(P+1,P)-D(P,P+1)+D(P,P).             (19)
```

The ridge `delta min(P,Q)` contributes `delta` to `R_P`.  Each curvature
atom contributes `-2min(t,P(1-t))`, because the Green term vanishes at
the two diagonal corners and equals that clipped moment at each
off-diagonal corner.  Combining this with (17) gives the finite exact
Möbius identity

```text
delta=R_P+2S_P                                     (20)
```

for every `P>=1`.  Hence the diagonal coefficient is not lost by the
mixed-difference representation: it is lost only when one discards the
diagonal cell or the rest of its chamber row.

## 5. Relation to THM-2348 and the sign convention

THM-2348 defines, for finite knots,

```text
mu_K(P,Q)
 =d_G(K,P#Q)-d_G(K,P)-d_G(K,Q)+u(K),                (21)
```

and at the unknot proves

```text
mu_U(P,Q)=-sigma(P,Q),                              (22)
```

where `sigma` is the nonnegative connected-sum interaction cocycle of
THM-2176.  The function `D` in (1) is the stable torus--mirror
counterpart of this positive defect.  Consequently `C_(P,Q)` is the
positive mixed increment on the `sigma` side, or the negative mixed
increment in the corresponding stabilized `mu_U` sign convention.

This is a representation bridge, not an identification of THM-2348's
finite target-allocation tables with the stable norm.  It explains why
the faithful object in both settings is symmetric cohabitation data:

```text
finite allocation:
  mixed ANOVA tensors / conditional quotient rectangles;

stable multiplicity chamber:
  positive curvature measure / clipped-tent mixed increments.         (23)
```

Neither object is intrinsically a tournament.

## 6. The diagonal-ridge boundary

Equation (14) is deliberately chamberwise.  If `kappa=0`, then

```text
h(z)=delta z,

D(P,Q)=delta min(P,Q)                               (24)
```

on the full positive quadrant.  Every interior chamber difference (6)
vanishes, yet for `delta>0` the diagonal ridge in (19) is still a genuine
cross-type interaction.  Therefore

```text
all C_(P,Q)=0

does not imply THM-2348 robust prime-type rectangularity.              (25)
```

One must also inspect a diagonal crossing—or equivalently use (20)—to
kill the balanced ridge `delta`.  THM-2348 further requires vanishing
under arbitrary typewise perturbations and conditioning, not only
multiplicity cells in this stable two-generator plane.

If both

```text
delta=0,                         kappa=0,            (26)
```

then the entire relaxed stable defect vanishes.  This is the only point
at which the curvature and ridge obstructions simultaneously disappear.

## 7. New decisive measurements

The formula suggests a concrete hierarchy for future knot calculations.

```text
one positive C_(P,Q):
  certifies nonzero interior curvature and rules out the linear profile;

all chamber C_(P,Q) zero:
  forces kappa=0 but leaves the balanced ridge delta;

one complete row plus its diagonal cell:
  recovers delta exactly by (20), at every finite scale;

finite target-token conditioning:
  is still needed to compare with THM-2348's robust allocation
  factorization.                                      (27)
```

Thus a small exact bank of connected-sum distances can attack the
profile shape before the balanced defect itself is known.  The theorem
does not assert that the required finite ordinary distances equal their
stable limits, nor that every abstract atomic profile is knot-realizable.
