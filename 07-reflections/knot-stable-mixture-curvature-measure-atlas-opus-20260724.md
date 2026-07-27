---
source: codex-2026-07-24-knot-curvature-profile
status: >
  PROVED + INDEPENDENTLY AUDITED ABSTRACT CONVEX GEOMETRY, SUPERSEDED /
  SHARPENED FOR THE ACTUAL T(2,7) KNOT. Conditional only on the proved
  THM-2308 stable-plane data, the complete same-sign mixture profile is
  equivalent to one canonical positive curvature measure on (0,1) with
  two exact moment budgets. The later Blanchfield--Gordian selector
  sharpens A<=6-delta to A<=4-delta and forces rigidity at delta=4.
  This file retains the exact classification of the weaker THM-2308-only
  information surface; it is not the current sharp knot-specific atlas.
depends_on:
  - THM-2308-mirror-double-nakanishi-floor-and-sharp-stable-mixture-profile
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2292-common-catalytic-section-and-helly-calibration-nerve
---

# The stable knot mixture is a curvature measure, not a tournament

> **SHARPENED ACTUAL-KNOT BOUNDARY.** The abstract classification below
> is exact for the information inherited from THM-2308 alone. For the
> actual knot `T(2,7)`, see
> `07-reflections/knot-blanchfield-mirror-mixture-selector-opus-20260725.md`:
> its explicit Gordian gauge gives `A<=4-delta`, adds the norm support
> `2|s|+4|t|`, and eliminates continuum ambiguity at `delta=4`.

## 1. Inheritance and the corrected object

For

```text
K=T(2,7),                 a=[K],             b=[mirror(K)],
```

THM-2308 writes the stable same-sign defect as

```text
D(P,Q)=3(P+Q)-u_hash(Pa+Qb),
D(P,Q)=P h(Q/P),                    P>=Q>=0,
delta=h(1),                         1<=delta<=4.       (1)
```

The profile `h` is concave and obeys

```text
delta z<=h(z)<=min(6z,delta(1+z)/2).                 (2)
```

The binary symbiont relation remembers only whether a defect is positive.
Even finitely many exact mixture values do not retain the shape in (2).
The faithful projective coordinate is the distributional curvature of
`h`.

## 2. Canonical curvature representation

Concavity and `h(z)<=6z` imply that

```text
a=h'_+(0)=lim_(z downarrow 0) h(z)/z
```

exists and is finite. Define the canonical nonnegative Radon measure

```text
mu=-D^2 h
```

on the **open** interval `(0,1)`. Then

```text
h(z)=a z-integral_(0,z)(z-t)dmu(t).                 (3)
```

Let

```text
A=integral_(0,1)(1-t)dmu(t),
B=integral_(0,1)t dmu(t).                           (4)
```

Evaluating (3) at one and differentiating from the left give

```text
a=delta+A,
h'_-(1)=delta-B.                                    (5)
```

The two support-line bounds in (2) are therefore exactly the moment
budgets

```text
A<=6-delta,
B<=delta/2.                                         (6)
```

Their slacks have direct geometric meanings:

```text
6-delta-A=6-h'_+(0),
delta/2-B=h'_-(1)-delta/2.                          (7)
```

Substituting (5) in (3) yields the positive Green-kernel form

```text
h(z)=delta z+integral_(0,1) G(z,t)dmu(t),            (8)

G(z,t)=min(z,t)(1-max(z,t)).
```

Here `G(z,t)>0` for `z,t in (0,1)` and

```text
-D_z^2 G(z,t)=delta_t.                              (9)
```

Consequently the pair `(delta,mu)` reconstructs `h`, and the canonical
curvature of `h` recovers `mu`.

## 3. Converse and reconstructed norm

Conversely, fix `0<delta<6` and any finite nonnegative measure `mu` on
`(0,1)` satisfying (6). Define `h` by (8). Then `h` is continuous and
concave, has endpoints `0,delta`, and

```text
h(z)>=delta z,

G(z,t)<=z(1-t)       -> h(z)<=6z,

G(z,t)<=(1-z)t       -> h(z)<=delta(1+z)/2.          (10)
```

Thus (2) holds.

Extend the defect symmetrically and homogeneously:

```text
D(P,Q)
 =max(P,Q) h(min(P,Q)/max(P,Q))

 =a min(P,Q)
   -integral_(0,1)(min(P,Q)-t max(P,Q))_+ dmu(t).   (11)
```

In THM-2308's `(e,o)` coordinates put

```text
X=|s|,                   Y=|t|.
```

The associated norm has the exact formula

```text
p_mu(s,t)=6Y,                                      Y>=X,

p_mu(s,t)
 =(6-a)X+aY
  +integral_(0,1)((1-u)X-(1+u)Y)_+ dmu(u),         X>=Y.       (12)
```

Each central-chamber summand is convex. At `X=Y`, both formulas equal
`6X`; the outward slope jump is legal because `a<=6`. The quadrant
function is coordinatewise nondecreasing:

```text
partial_X p_mu>=6-a>=0,

min partial_Y p_mu
 =a-integral(1+u)dmu(u)
 =delta-2B>=0.                                      (13)
```

Quadrant convexity plus coordinatewise monotonicity proves that
`p_mu(|s|,|t|)` is an unconditional convex homogeneous function. Moreover

```text
p_mu(s,t)>=max((6-delta)|s|,6|t|),                  (14)
```

so it is a norm. It has exactly

```text
p_mu(e)=6-delta,       p_mu(o)=6,       p_mu(a)=p_mu(b)=3,

p_mu(s,t)=6|t|                         when |t|>=|s|.          (15)
```

For the knot range `delta<=4`, (2) also retains THM-2308's common
Alexander-fibre floor. Hence every measure satisfying (6) produces a
fully admissible abstract stable-plane norm with all currently proved knot
data.

## 4. The two sharp envelopes are curvature extremals

The lower defect envelope is

```text
mu=0,
h(z)=delta z.                                       (16)
```

It reconstructs THM-2308's upper norm

```text
p_high(s,t)
 =(6-delta)max(|s|,|t|)+delta|t|.                   (17)
```

At the other extreme put

```text
z_0=delta/(12-delta),
m_0=(12-delta)/2,
mu=m_0 delta_(z_0).                                 (18)
```

Both budgets are saturated:

```text
m_0(1-z_0)=6-delta,
m_0 z_0=delta/2.                                    (19)
```

The resulting profile is

```text
h(z)=min(6z,delta(1+z)/2),                           (20)
```

and the norm is THM-2308's lower envelope

```text
p_low(s,t)=max((6-delta)|s|,6|t|).                  (21)
```

The atom in (18) is the exact strategy-transition kink; its mass is the
slope drop `6-delta/2`.

## 5. Endpoint caveat and finite-relation no-go

The measure must be canonicalized on `(0,1)`. Atoms at zero and one are
invisible because

```text
G(z,0)=G(z,1)=0.                                    (22)
```

Thus a measure on the closed interval is a valid redundant
parameterization, but it is not faithful. Also `mu` alone is insufficient:
`delta` must be retained.

On the deliberately weaker **THM-2308-only** information surface, for
every fixed `delta in [1,4]`, continuum many feasible measures give
distinct profiles while preserving every axis and chamber-wall value.
This statement is superseded at `delta=4` for the actual knot by the
Blanchfield--Gordian selector. On the abstract surface, for example,

```text
mu_c=(1/2)delta_c,             1/4<=c<=3/4,          (23)
```

has strict slack in both budgets.

More strongly, fix any finite sample `z_1,...,z_n`. Choose `n+3` distinct
points `t_j in (0,1)`. The `(n+2) x (n+3)` matrix with rows

```text
G(z_i,t_j),          1-t_j,          t_j             (24)
```

has a nonzero kernel vector with both signs. A sufficiently small signed
perturbation of equal positive atomic weights gives two distinct feasible
measures `mu_+` and `mu_-` with:

```text
the same delta;
the same two curvature moments A,B;
the same values at every sampled ray z_i;
the same axes, walls, and outer law;
different profiles somewhere.                       (25)
```

The last line follows because `-D^2` recovers the canonical measure.
Therefore no fixed finite collection of sampled exact gains determines the
stable profile under the present axioms. Tournament edges, binary symbiont
signs, and any other finite comparison shadow are strictly coarser.

## 6. Information ledger

```text
source:
  THM-2308's unconditional stable norm and concave defect profile;

map:
  take the balanced defect delta and distributional curvature on (0,1);

preserved:
  every same-sign stable mixture value and, with the outer-law/sign
  sidecar, the complete norm on span{K,mirror(K)};

destroyed:
  raw unknotting numbers, integrality, catalysts, knot representatives,
  and any labelled sign relation not supplied separately;

canonical target:
  (outer-law/sign sidecar, delta, mu);

cheapest knot-specific test:
  compute one nontrivial curvature moment, one interior atom, or one
  interior rational slope by a recognizable additive lower bound or an
  explicit repeated-sum upper certificate.                           (26)
```

This is a complete classification of the **abstract stable-plane
uncertainty left by THM-2308**, not a classification of knots and not a
new value of `u` or `u_hash`.
