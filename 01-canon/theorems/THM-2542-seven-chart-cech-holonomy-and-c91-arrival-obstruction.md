---
id: THM-2542
title: "Seven-chart Cech holonomy, the C91 mapping torus, and the semantic-arrival obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The theorem computes the finite root-deck local system already present in
  THM-2535.  It does not identify its abstract C91 mapping torus with a
  physical C91 endpoint/current, does not construct the semantic 2-cell, and
  does not remove an LRC(14) row.
source: codex-2026-07-27-holotopy-groupoid
depends_on:
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse
related:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2418-alternating-base-thirteen-septimal-carry-matrix-and-rank-one-boundary
  - THM-2424-coprime-common-root-crt-and-unit-residue-spectrum
  - THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2539-diagonal-cubic-owner-clock-boundary-current
  - THM-2540-weighted-live-event-kakeya-flux-and-transverse-gain-boundary-refinement
  - THM-2541-canonical-typed-row-full-target-plane-support
script: 04-computation/lrc14_seven_chart_cech_holonomy_thm2542.py
output: 05-knowledge/results/lrc14_seven_chart_cech_holonomy_thm2542.out
script_sha256: 97aac5f1f9a600a93ce5d14f68c89402ad6f6633950211803c8850c521edd347
output_sha256: fdc49019fa0a3fde0fe8f8406e5429411ec3da6c26abd20f71f8d7c99d5dc666
hash_basis: working-tree bytes (LF)
---

# THM-2542 -- seven-chart Cech holonomy and the C91 arrival obstruction

THM-2535's seven transported selector/scheduler charts do not merely fail to
share a preferred root origin.  Their overlap maps carry a nonzero Cech class
in

```text
H^1(C_7^graph;F_13)=F_13.                                   (1)
```

Here `C_7^graph` is the seven-cycle nerve and (1) is simplicial/Cech
cohomology of its geometric realization.  It is **not** cyclic-group
cohomology of the group `C_7` with trivial coefficients.

For marker root `a!=0`, the class is `7a`.  It cannot be removed by changing
the seven chart origins.  Its minimal trivializing clock cover has degree
thirteen, and the corresponding root-clock mapping torus is one connected
cycle of length ninety-one.

This is a genuine spatial holonomy theorem.  It is not semantic arrival:
there is an exact equal-mass hostile with all primitive mixed modes, the same
nonzero `7a` chart class, and no target-active vertical edge.  The faithful
carrier is therefore a path/double category: invertible horizontal arrows are
root/clock chart transports, directed vertical arrows are semantic
source-to-arrival transitions, and lawful physical couplings are 2-cells.
One possible root-trivializing 2-cell would have to contribute total
correction `-7a`; after a lawful intertwiner which kills target-neutral roles,
at least one target-active role would be forced.  A physical construction may
instead retain the twisted C91 bundle.  Neither kind of 2-cell is built here.

## 1. The seven-chart nerve

Fix one positive marker class `a in F_13^*` in THM-2535 Section 6.  For every
clock `k in F_7`, let

```text
X_k=E_a intersection P_k,
x_k=measure(X_k)>0.                                         (2)
```

At every sufficiently large admissible delay the seven `X_k` are disjoint
positive rational Boolean events.  THM-2535 evaluates the adjacent clock
boundary in the transported marker-to-target chart as

```text
R_k(v)=x_k 1_(v=a)-x_(k+1)1_(v=0),          v in F_13.       (3)
```

The coefficient/mass of `X_k` is represented at cut-output root `0` in the
incoming chart `k-1` and at cut-output root `a` in the outgoing chart `k`.
These are chart coordinates, not two asserted physical inverse-root locations
of the event.  The overlap transition from the incoming root coordinate to the
outgoing one is therefore

```text
g_k=a in F_13.                                               (4)
```

The nerve is the oriented seven-cycle.  Thus `g=(g_k)` is a Cech one-cochain
with values in the additive root-deck group `F_13`.

Changing the root origin in chart `k` by `h_k in F_13` replaces

```text
g_k -> g_k+h_k-h_(k-1).                                     (5)
```

Hence the cyclic sum

```text
Hol(g)=sum_(k in F_7)g_k=7a in F_13                         (6)
```

is gauge invariant.  Since `a!=0` and `13` does not divide `7`,

```text
Hol(g)!=0.                                                  (7)
```

On the cycle graph, the sum map identifies `H^1(C_7^graph;F_13)` with
`F_13`; therefore
(7) is exactly the obstruction to a common root trivialization.  Equivalently,
the equations

```text
h_k-h_(k-1)=-a,                  k in F_7,                   (8)
```

sum to the contradiction `0=-7a`.

For every nontrivial root character `alpha`, the associated Wilson phase is

```text
zeta_13^(-7 alpha a)!=1.                                   (9)
```

Thus all twelve root characters see the same nontrivial class.

## 2. The minimal trivializing cover is thirteen-fold

Pull (4) back to the `n`-fold cover of the clock cycle.  Its holonomy is

```text
n Hol(g)=7na.                                               (10)
```

It vanishes exactly when `13` divides `n`.  Hence the smallest positive
trivializing degree is thirteen.

There is an equivalent mapping-torus description.  On

```text
F_7 x F_13
```

use the skew successor

```text
(k,r)->(k+1,r+a).                                           (11)
```

The element `(1,a)` has order `lcm(7,13)=91`.  Since the state space has
exactly ninety-one elements, (11) is one connected ninety-one-cycle.  After
seven clock edges its root coordinate has moved by `7a`; only after thirteen
clock turns does it close.

This explains, rather than assumes, why the natural chart carrier is a C91
object.  It does **not** identify (11) with THM-2424's physical CRT packet,
THM-2419's relation shell, a `91`-unit ordinary frequency, or an actual
terminal-current endpoint.  Those objects require their own typing maps.

### 2.1 General two-root/two-clock classification

The calculation is not tied to the adjacent clock or the chord `a->0`.  Fix

```text
eta in F_7^*,                  s,t in F_13, s!=t.             (11a)
```

Use the clock boundary and its affine root chart

```text
D_k(h,r)
 =epsilon 1_(h=s)(1_(r=k)-1_(r=k+eta)),

phi_k(h,r)=h+(t-s)eta^(-1)(r-k).                            (11b)
```

The chart sends `(s,k)` to `s` and `(s,k+eta)` to `t`, so its pushforward is

```text
R_k=epsilon(delta_s-delta_t).                               (11c)
```

The event at clock `k` is at root `t` in the incoming chart `k-eta` and at
root `s` in the outgoing chart `k`.  Hence the mapping-torus step is

```text
(eta,s-t) in F_7^* x F_13^*.                                (11d)
```

Its Cech class is `7(s-t)!=0`, and (11d) again generates one cycle of order
ninety-one.  The primitive transforms are

```text
Dtilde_k(alpha,beta)
 =epsilon zeta^(-alpha s)xi^(-beta k)
          (1-xi^(-beta eta))!=0,

sum_k Rhat_k(alpha)
 =7epsilon(zeta^(-alpha s)-zeta^(-alpha t))!=0.              (11e)
```

Thus the C91 conclusion is intrinsic to every nonzero clock step and every
oriented pair of distinct roots.  The THM-2535 packet is
`(eta,s,t)=(1,a,0)`.

### 2.2 Coprime-cycle holotopy lemma

The mechanism is reusable.  Let `p!=q` be primes.  On a `q`-cycle of charts
with nonzero clock step `eta in F_q` and nonzero overlap translation
`delta in F_p`, the constant transition cochain has

```text
[g]=q delta in H^1(C_q^graph;F_p).                           (11f)
```

It is nonzero, its minimal trivializing base cover has degree `p`, and the
mapping-torus successor

```text
(k,r)->(k+eta,r+delta)                                       (11g)
```

is one cycle of length `pq`.  The proof is exactly (5)--(11): a gauge adds a
telescoping coboundary, an `n`-fold pullback multiplies (11f) by `n`, and
`(eta,delta)` has order `pq` in `F_q x F_p`.  For an oriented pair of distinct
roots, the equal-mass boundary has every primitive mixed character because
neither `1-xi_q^(-beta eta)` nor the difference of two distinct `p`th-root
characters can vanish.  Sections 1--3 are the instance `(p,q)=(13,7)`.

## 3. Full spectral holonomy still need not supply a semantic edge

The distinction has a sharp rational finite hostile.  Take seven disjoint
Boolean intervals of one common rational mass `0<epsilon<=1/7` and put
`x_k=epsilon`.
Before applying the transported root chart, the local signed clock table is

```text
D_k(h,r)
 =epsilon 1_(h=a)(1_(r=k)-1_(r=k+1)).                       (12)
```

With `zeta=zeta_13`, `xi=zeta_7`, its unnormalized mixed transform is

```text
Dtilde_k(alpha,beta)
 =epsilon zeta^(-alpha a)xi^(-beta k)(1-xi^(-beta))!=0      (13)
```

for every `alpha in F_13^*`, `beta in F_7^*`, and every chart `k`.
Nevertheless

```text
sum_k D_k=0.                                                (14)
```

Applying each transported chart first gives

```text
R_k=epsilon(delta_a-delta_0),

sum_k R_k=7epsilon(delta_a-delta_0),                         (15)
```

whose twelve nontrivial root coefficients are all nonzero.  Equations
(12)--(15) simultaneously realize:

```text
all 504 chart-labelled primitive mixed modes for fixed a;
nonzero root-deck holonomy 7a;
all twelve charted root colours;
zero uncharted clock boundary.                              (16)
```

Across all twelve nonzero marker roots this is `6,048` local modes and `144`
charted root modes.

Now attach no semantic vertical edge to any empty head (equivalently, set an
auxiliary target-active incidence indicator to zero).  None of (2)--(16), nor
any axiom used to derive them, changes.  This is the canonical equal-mass
seven-chart hostile to

```text
nonzero root-clock chart holonomy
  => positive target-active semantic arrival.               (17)
```

It is a lawful finite local-system/Boolean-atom hostile, not a claimed full
fourteen-speed scalar cover.  Its purpose is to locate the first missing map:
the current chart class has horizontal spatial monodromy while the semantic
vertical-arrow set may still be empty.  Chronological arrival is directed and
is not being mis-typed as another invertible `H^1` coordinate.

This can be recorded as a bidegree audit.  Give root/clock chart transports
horizontal degree `(1,0)` and a directed semantic source-to-arrival arrow
vertical degree `(0,1)`.  Gauge changes, finite Fourier transforms, Cayley
inversion, scheduler summation, and the chart pushforwards used above all act
inside vertical degree zero.  Their sums and compositions therefore remain in
the linearized horizontal path category (equivalently, the additive
enrichment of the horizontal arrows).  A lawful physical coupling between a
selected wall and an arrival is new data of bidegree `(1,1)`; it cannot be
manufactured by iterating the existing degree-`(*,0)` operators.  The hostile
is the exact nonzero-horizontal/empty-vertical control for this typing
statement.

## 4. What one root-trivializing semantic 2-cell would have to do

The obstruction gives a precise conditional positive statement for one
possible closure architecture.  Suppose a future common-ancestry construction
assigns an additional root correction

```text
c_k in F_13                                                   (18)
```

to the semantic edge over clock overlap `k`, and suppose its filled
root/semantic 2-cell makes the combined transitions a coboundary.  Then

```text
0=sum_k(g_k+c_k)=7a+sum_k c_k,

sum_k c_k=-7a!=0.                                           (19)
```

Consequently at least one semantic vertical edge has `c_k!=0`.

THM-2461 proves that four of the five first-failure roles are target-neutral
and only one is target-active.  To apply (19) physically one still needs a
lawful affine intertwiner

```text
pi: target-role covectors -> root-deck corrections            (20)
```

which sends target-neutral roles to zero and the relevant target-active
direction to a nonzero correction.  If a common-ancestry 2-cell supplies
(20) **and trivializes the horizontal root chart**, then (19) forces at least
one target-active edge.  If the edge weights
are the positive Boolean incidences of THM-2537, that edge has positive mass.

This is sufficient, not necessary: a physical ancestry lift may live directly
on the nontrivial C91 mapping torus instead of killing its class.  Neither
(20), a root-trivializing filled 2-cell, nor a positive vertical lift on the
twisted bundle is constructed here.  In particular, the root deck,
THM-2461's two-dimensional target covector, the visible scheduler
clock, and the inherited THM-2449 deep/source sheet remain distinct typed
coordinates.  Equation (19) is the exact invoice a proposed intertwiner must
pay; it is not permission to identify them.

## 5. Positive common-clock itineraries are not the missing 2-cell

There is one useful positive survivor.  Let `X_0,...,X_6` be arbitrary
positive rational BV subsets of the circle, and let `T(x)=13x`.  There are
arbitrarily separated integers

```text
0=n_0<n_1<...<n_6                                            (21)
```

such that

```text
measure(intersection_k T^(-n_k)X_k)>0.                       (22)
```

Indeed, start with `A_0=X_0`.  If `A_j` has positive measure, mixing of the
expanding map gives

```text
measure(A_j intersection T^(-N)X_(j+1))
  ->measure(A_j)measure(X_(j+1))>0.                           (23)
```

Choose an arbitrarily large `N=n_(j+1)` for which the intersection is
positive and continue.  Every intermediate set remains a rational BV set.

Thus the seven selector/scheduler cells can be visited by one genuine orbit
itinerary.  If one repeats `X_0` at an eighth stop, the same induction gives a
cyclic clock itinerary and positive graph couplings on all seven clock edges.
The exact companion includes the positive cylinder whose first seven
base-thirteen digits are `0,1,...,6`; its measure is `13^(-7)`.

Equation (22) solves only common positive mass on one orbit.  It does not
force the root transition (4) at those times, preserve one old deep sheet,
match ordinary frequencies, identify a target-active role, or fill the
semantic 2-cell of Section 4.  Mixing supplies vertices of the desired
holotopy; it does not supply its labelled face.

## 6. Exact stopping boundary

The proved path is now

```text
seven transported selector/scheduler charts
  -> nonzero Cech class 7a
  -> minimal degree-13 trivializing clock cover
  -> one abstract C91 mapping-torus cycle
  -> exact positive common-orbit itinerary.                  (24)
```

The missing arrow is

```text
abstract C91 chart cycle
  -/-> target-active same-ancestry semantic 2-cell
       retaining the old deep/source carrier.                (25)
```

The equal-mass hostile proves that another Fourier-mode count cannot create
that arrow.  One closing architecture is a flat fill: build the affine map
(20) on the physical packet and a positive target-active overlap whose
boundary cancels `7a`.  The other is a twisted lift: construct a positive
semantic vertical path directly on the C91 mapping torus, preserving the old
deep/source carrier and transporting rather than cancelling its holonomy.
THM-2538's anchored cross-Kakeya program and
THM-2540's transverse product packing remain relevant only if they preserve
that same ancestry/deep face.

No terminal arrival, scalar row exclusion, or proof of LRC(14) follows.

## 7. Exact referee

Run

```bash
python3 04-computation/lrc14_seven_chart_cech_holonomy_thm2542.py
python3 -O 04-computation/lrc14_seven_chart_cech_holonomy_thm2542.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_seven_chart_cech_holonomy_thm2542.out
```

The dependency-free referee works in `F_547`, which contains primitive
seventh and thirteenth roots.  It checks the full general bank of `936`
triples `(eta,s,t)`, all `6,552` basis-gauge changes, minimal trivializing
degree thirteen, all `936` connected ninety-one-cycle mapping tori, all
`471,744` primitive local modes in (11e), all `11,232` charted root modes,
the twelve THM-2535 special classes, all `85,176` uncharted cancellation
cells, the empty semantic vertical-edge control, and the seven-stop rational
digit cylinder.

**QED.**
