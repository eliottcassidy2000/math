---
source: codex-2026-07-25-knot-curvature-extremes
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  For every interior
  torus--mirror defect 0<delta<2(g-1), the curvature-measure relaxation
  inherited from the exact axes and the complete two-frequency
  Blanchfield selector bank has an explicit extreme-point atlas.  Its
  only extreme measures are zero; one atoms whose mass saturates at
  least one of the two moment budgets; and two atoms straddling the
  canonical barycentre delta/(4(g-1)-delta) and saturating both budgets.
  In particular the inherited ambiguity remains continuum-sized even
  among extreme profiles.  This classifies the abstract relaxation, not
  the actual stable Gordian norm and not realizability by knots.
depends_on:
  - knot-torus-mirror-selector-stable-plane-opus-20260725
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2346-global-allocation-anova-normal-form
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
  - knot-curvature-mixed-cohabitation-kernel-opus-20260725
---

# Extreme curvature profiles on the torus--mirror plane

## 1. The relaxation as a two-moment cone

Fix

```text
g>=2,                       L=2(g-1),

0<delta<L.                                           (1)
```

The audited torus--mirror plane theorem writes every inherited defect
profile in the Green form

```text
h(z)=delta z+integral_(0,1) G(z,t) d kappa(t),

G(z,t)=min(z,t)(1-max(z,t)),                          (2)
```

where `kappa` is a finite nonnegative Borel measure on `(0,1)`.  The two
sharp support inequalities are exactly

```text
integral (1-t)d kappa(t) <= A_0=L-delta,

integral t d kappa(t)     <= B_0=delta/2.            (3)
```

Call this measure set `K_delta`.  We keep the open interval in the
definition.  Endpoint atoms would consume a budget while contributing
zero to the Green kernel, so including them without normalization would
introduce a fake nonuniqueness unrelated to profiles.

The map `kappa -> h` is affine and injective on `(0,1)`, because

```text
-D^2 h=kappa                                           (4)
```

in distributions.  Thus it is enough to classify the extreme points of
`K_delta`.

Put

```text
M_0=A_0+B_0=L-delta/2,

t_0=B_0/M_0=delta/(2L-delta).                        (5)
```

The point `t_0` is the barycentre forced when both inequalities in (3)
are equalities.

There is also a useful compact normalization.  For an interior measure
`kappa`, write its moments as `A(kappa),B(kappa)` and set

```text
kappa_bar
 =kappa+(A_0-A(kappa))delta_0+(B_0-B(kappa))delta_1. (5a)
```

Then `kappa_bar` lies in the saturated slice

```text
Lambda_delta
 ={lambda>=0 on [0,1]:
   A(lambda)=A_0, B(lambda)=B_0}.                    (5b)
```

Conversely, restriction from `[0,1]` to `(0,1)` reverses (5a).
Because `G(z,0)=G(z,1)=0`, this is an affine identification of the
profile class with `Lambda_delta`.  Every member of `Lambda_delta` has
total mass `M_0` and barycentre `t_0`; the slice is weak-star compact.
This removes endpoint slack without pretending that invisible endpoint
mass is curvature.

## 2. Complete extreme-point atlas

The extreme points of `K_delta` are precisely the following.

### Zero

```text
kappa=0.                                             (6)
```

### One atom

For any `0<t<1`,

```text
kappa=m(t) delta_t,

m(t)=min((L-delta)/(1-t), delta/(2t)).                (7)
```

At least one budget is tight.  More explicitly,

```text
t<t_0:  m(t)=(L-delta)/(1-t);

t=t_0:  m(t)=M_0;

t>t_0:  m(t)=delta/(2t).                              (8)
```

The middle atom is the unique one-atom profile saturating both budgets.

### Two atoms

Choose

```text
0<s<t_0<t<1.                                         (9)
```

Then

```text
kappa=m_s delta_s+m_t delta_t,

m_s=M_0 (t-t_0)/(t-s),

m_t=M_0 (t_0-s)/(t-s).                              (10)
```

Both masses are positive and both budgets in (3) are tight.  Every
two-atom extreme point has this form.

There are no other extreme points.

## 3. Proof

First suppose a feasible measure has support on at least three points.
Choose three pairwise disjoint Borel sets `E_1,E_2,E_3`, each of positive
`kappa`-mass.  Their two moment vectors

```text
v_i=(integral_(E_i)(1-t)d kappa,
     integral_(E_i)t d kappa) in R^2                 (11)
```

are linearly dependent.  Choose nonzero real `c_i` with
`sum_i c_i v_i=0`, and set

```text
nu=sum_i c_i kappa|_(E_i).                           (12)
```

For sufficiently small `epsilon>0`, both
`kappa+epsilon nu` and `kappa-epsilon nu` are nonnegative.  They have
exactly the same two moments as `kappa`, are feasible, and are distinct.
Thus an extreme measure has support size at most two.  The same argument
applies to non-atomic support by first choosing the three positive-mass
sets.

For one atom, if both inequalities are strict, its mass can be perturbed
up and down.  Hence an extreme one-atom measure must have the maximal
mass in (7).  Conversely, if one budget is tight and

```text
m delta_t=(nu_1+nu_2)/2,                             (13)
```

positivity forces both `nu_i` to be supported at `t`.  Tightness forces
both masses to equal `m`, so (7) is extreme.  Comparing the two entries
in (7) gives (8).

For two positive atoms, if fewer than two budgets are tight, the two
masses admit a nonzero small perturbation preserving every tight
constraint and keeping each slack constraint feasible.  Therefore both
budgets must be equalities.  Adding the two equalities gives total mass
`M_0`; the second gives mean `t_0`.  Positive masses exist exactly when
the two support points strictly straddle `t_0`, and solving the two by
two system gives (10).

Conversely, if a measure in (10) is the midpoint of two feasible
measures, positivity restricts both summands to `{s,t}`.  Since both
moments of the midpoint saturate (3), both summands saturate them as
well.  The two moment vectors at distinct `s,t` are independent, so the
two masses are uniquely (10).  The summands equal the midpoint, proving
extremality.

Finally, zero is extreme by positivity.  This exhausts all cases.

## 4. Stronger ambiguity and the rigid endpoints

The earlier note used `kappa=0` and the special atom

```text
M_0 delta_(t_0)                                     (14)
```

to prove nonuniqueness.  The atlas is stronger: for every
`t in (0,1)`, (7) gives a distinct extreme Green profile, and every pair
straddling `t_0` gives another one through (10).  Hence at every interior
`delta` the relaxation has continuum many extreme points, not merely
continuum many convex interpolants between two endpoints.

Explicitly, the extreme profiles are

```text
h_0(z)=delta z;                                      (14a)

h_t(z)=delta z+(A_0/(1-t))G(z,t),   0<t<t_0;        (14b)

h_(t_0)(z)=delta z+M_0G(z,t_0)
          =min(Lz,delta(1+z)/2);                     (14c)

h_t(z)=delta z+(B_0/t)G(z,t),       t_0<t<1;        (14d)
```

together with

```text
h_(s,t)(z)
 =delta z
  +M_0 (t-t_0)/(t-s) G(z,s)
  +M_0 (t_0-s)/(t-s) G(z,t),

0<s<t_0<t<1.                                        (14e)
```

In (14e), the left obstacle is attained for every `z<=s` and the right
obstacle for every `z>=t`:

```text
h_(s,t)(z)=Lz                    for z<=s,

h_(s,t)(z)=delta(1+z)/2          for z>=t.           (14f)
```

The one-atom profiles similarly touch exactly the outer obstacle whose
budget they saturate.

The translation to the inherited mirror-symmetric norm class is affine
and injective.  For `X=|s|,Y=|t|`,

```text
p_h(se+to)
 =2gY,                                      Y>=X;

 =2gX-(X+Y)h((X-Y)/(X+Y)),                 X>=Y.    (14g)
```

Consequently (14a)--(14e), and only those profiles, give the extreme
norms **within this inherited curvature relaxation**.  This is not an
assertion that any nontrivial member is realized by a knot family.

The boundary cases recover the two rigid profiles without a limiting
argument:

```text
delta=0:
  integral t d kappa=0, hence kappa=0 on (0,1);

delta=L:
  integral (1-t)d kappa=0, hence kappa=0 on (0,1).  (15)
```

For `g=1`, `L=0` and there is no interior defect interval.

## 5. Exact scope

This theorem classifies the convex relaxation forced by:

```text
the two calibrated axes;

mirror symmetry and homogeneity;

the complete two-frequency inertia selector polygon;

the two Green moment budgets.                       (16)
```

It does not say that every extreme profile is realized by the stable
Gordian norm of a knot family.  The missing sidecar is geometric
realizability: one needs upper constructions or new lower selectors that
select a particular curvature measure.

THM-2346 gives a complementary warning for catalytic allocation
problems: symmetric cohabitation tensors do not become tournament
orientations merely by taking signs.  Here the faithful object is again
convex rather than tournament-like—the two-moment curvature cone and
its extreme atomic atlas.
