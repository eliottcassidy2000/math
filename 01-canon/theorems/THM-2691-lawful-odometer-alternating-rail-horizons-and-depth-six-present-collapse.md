---
id: THM-2691
title: "Lawful odometer alternating rail horizons and depth-six present collapse"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The lawful
  alternating THM-2657 lifts +/-(13^5+1) support, for every finite H, a
  selected H-event THM-2584 rail/clock cylinder of exact length
  1/(7*13^(H+1)).  After dynamically typing and intersecting the THM-2640
  current-present factor at every event, the selected cylinders remain
  positive through H=5 but are empty at H=6.  The H=5 all-event survivor has
  an explicit open atom of length 13^(-11) and future-digit word 7,0,0,0,0.
  The H=6 cylinder lies strictly in the plus-state present-free gap consisting
  of the missing future-digit cells 5 and 6.  This closes the selected central cycle,
  not other lifts, rails, components, mixed D/slope faces, or LRC(14).
source: root/lrc-physical-gluing-probe-2026-07-28
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2689-affine-clock-support-typing-tradeoff-and-odometer-phase-locality
related:
  - THM-2684-three-tooth-rail-envelope-diagonal-arrival-law-and-full-dilation-nilpotence
script: 04-computation/lrc14_lawful_odometer_alternating_rail_horizon_probe.py
output: 05-knowledge/results/lrc14_lawful_odometer_alternating_rail_horizon_probe.out
script_sha256: 39ea72b7387b4ea06a89f958db155dc98f4674df7824cadcbdab116a0f6168fa
output_sha256: 9a4854701a4b348adeeb2b75210ac525a7188a3da15de0bb89cccc537d17ffc3
secondary_script: 04-computation/lrc14_lawful_odometer_alternating_unit_line_boundary.py
secondary_output: 05-knowledge/results/lrc14_lawful_odometer_alternating_unit_line_boundary.out
secondary_script_sha256: 66722656a5e7fb42b87a1c97fce8d483a370e7ae06c72027ec0eb5cbf0fcc9ae
secondary_output_sha256: 5445d7840e66b12d256de5b5a1e0c92d52102fc1b47472ef86f59d252a7b43e4
hash_basis: LF-normalized bytes
---

# THM-2691 -- the lawful central odometer cycle dies exactly at the sixth present layer

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2689 proves that every nonzero physical odometer lift escapes the
unshifted three-event support zero, while no such lift transports the seven
clock labels globally.  The most symmetric local escape is the alternating
central cycle.  This theorem follows that cycle through its actual labelled
THM-2584 rails and then inserts the first missing THM-2640 sidecar.

The outcome is a sharp depth invoice:

```text
rail + intrinsic clock:        positive at every finite horizon;
all dynamically typed
current-present factors:       positive through H=5, empty at H=6.       (1)
```

Thus an arbitrarily long rail germ is not an arbitrarily long physical
packet.  The failure occurs exactly when the repelling dynamics resolves the
depth-six present coordinate.

## 1. Lawful alternating cycle

Put

```text
p=13,              R=p^6,              S=p^5,
K=S+1,             T_-(x)={px-K/R},     T_+(x)={px+K/R}. (2)
```

The two points

```text
x_+=1/2+K/(14R)=4879851/9653618,
x_-=1/2-K/(14R)=4773767/9653618                         (3)
```

satisfy

```text
T_-(x_+)=x_-,                 T_+(x_-)=x_+.             (4)
```

Their intrinsic shallow-to-owner clock edges are

```text
x_+: 4->3,                    x_-: 3->4,                (5)
```

so each owner is the next event's shallow label.  Their predecessor carries
are `7,5`, their future digits at the cycle points are both `6`, and their
THM-2657 quotient root steps are `11,2`.  The lift magnitude `K` is a
thirteenth-unit, so these are lawful nonzero quotient lifts.  By THM-2689 the
cycle is phase-local rather than the restriction of a global clock action.

The two selected THM-2584 rail keys are

```text
R_+=(6,1,3,0),                 R_-=(6,1,4,12).           (6)
```

At the cycle points their unique local support pieces are

```text
I_+=(1195/2366,171/338),
I_-=(167/338,1171/2366),                                (7)
```

with the same positive weight `27,580,222,516`.  They lie inside the shallow
clock cells `4` and `3`, respectively.  The intervals in `(7)` are reflections
of one another.

## 2. Every finite rail horizon is positive

Write the initial point as `x=x_++e`.  As long as the orbit stays in the
local branches `(7)`, its state at time `n` is

```text
x_n=x_+ +13^n e   for n even,
x_n=x_- +13^n e   for n odd.                            (8)
```

The one-sided margins of `x_+` in `I_+` are

```text
M_out=14281/(7R),                M_in=2040/R.            (9)
```

The margins at `x_-` are interchanged.  Since

```text
13 M_in>M_out,                 13 M_out>M_in,            (10)
```

the last pulled-back rail is binding at every horizon.  The exact `H`-event
initial cylinder is

```text
C_H=(x_+-M_out/13^(H-1), x_++M_in/13^(H-1))  if H odd,
C_H=(x_+-M_in /13^(H-1), x_++M_out/13^(H-1)) if H even. (11)
```

Therefore

```text
length(C_H)=1/(7*13^(H+1))>0                   for all H>=1. (12)
```

The cylinders are nested and shrink to `x_+`; the infinite-horizon
intersection is the repelling cycle point, not an open interval.  The exact
companion checks `(11)` through `H=128`.  Independently, a literal weighted
common-grid pullback at `H=3` gives one component with the same endpoints and
weight product.  The packet transported so far consists of the THM-2584
b-word depth-five product profile, source, arrival, owner, absolute deep
label, and explicit shallow clock.  It does not yet contain the THM-2640
present or private-unit data.

## 3. Why the reflected coefficient units do not transport

Both cycle points lie on the future half wall

```text
{R x_+}={R x_-}=1/2.                                    (13)
```

Rebuild the two surrounding THM-2640 rail rows with metadata `(1,3,0)` and
`(1,4,12)`, global content `26`, and carries `7,5`.  The four coefficient
profiles recorded in THM-2689 obey the exact reflected scalar laws

```text
Y_- =10 rho(Y_+)  on the positive side,
Y_- = 4 rho(Y_+)  on the negative side.                 (14)
```

But the private half-edge used in `(14)` misses a sufficiently near
one-sided cycle point.  The opposite edge contains the point and is again a
unit row.  Its normalized profiles are

```text
positive side: (0,0,0,0,9,5,4), (0,3,9,9,0,0,0),
negative side: (0,0,0,0,4,4,10),(0,9,8,4,0,0,0).       (15)
```

For neither pair in `(15)` is the minus profile a nonzero scalar multiple of
the clock reversal of the plus profile.  Thus the scalar covariance and the
local private geometry choose opposite edge charts.  This is already a
pointwise failure of the tempting coefficient transport, although the local
replacement rows remain algebraic units.

## 4. Present-free gap at the plus state

At fixed future digit `h=6`, the plus cycle point is outside all seven
THM-2640 present-clock factors.  Its exact open gap is

```text
G_6=(202105/399854,202131/399854).                       (16)
```

The minus point is not symmetric: it lies in the present factor with clock
`3`.  This is a load-bearing hostile control against reflecting `(16)`.

If `h=6` is frozen, `C_2` still has positive present intersection, while
`C_3` lies strictly inside `(16)` and has none.  Freezing the centre's digit
across a neighbourhood is not physically valid, however.  Restore the
pointwise digit

```text
h(x)=floor(13 {R x}).                                    (17)
```

The dynamically typed present condition at one event is

```text
x belongs to union_(ell in F_7) F_(ell,-h(x)).           (18)
```

Around `x_+`, the exact dynamically typed gap is

```text
G=(x_+-3/(26R),x_++1/(26R))
 =(31719030/62748517,31719032/62748517).                 (19)
```

It is formed precisely by the adjacent missing future digits `5,6`; its
right endpoint is the first `h=7` support boundary.

## 5. All-event present census and sharp sixth-horizon zero

For every `x in C_H`, generate its alternating states `x_0,...,x_(H-1)` and
impose `(18)` at every state.  Exact common-grid intersection gives

```text
H       components       exact positive mass
1          6380          7021099/12298709332
2          4298          687145/22840460188
3          3207          780376/519620469277
4          1645          67780/965009442943
5            68          1259/7168641576148
6             0          0.                              (20)
```

The positive `H=5` row is not an endpoint artifact.  It contains the strict
open atom

```text
(905927272952/1792160394037,
 905927272953/1792160394037),                            (21)
```

of length `13^(-11)`, with future-digit word

```text
(7,0,0,0,0).                                             (22)
```

The entire `H=6` rail cylinder is

```text
(1811854513263/3584320788074,
 12682981649963/25090245516518).                         (23)
```

It lies strictly inside the state-zero gap `(19)`.  Quantitatively, the left
endpoint of the last `H=5` state-zero survivor exceeds the right endpoint of
`C_6` by

```text
171365/25090245516518>0.                                 (24)
```

Thus the sixth-horizon product is empty already at its current present
factor; later event factors cannot restore it.  Nesting of the cylinders
then makes every `H>=6` product empty.

The calculation of `(20)` is exhaustive on an integer grid, not a midpoint
test.  At horizon `H`, use the common grid `13^H T`.  For every state, the
future digit `h` occupies cells

```text
[(13j+h)u,(13j+h+1)u),                 u=grid/(13R).      (25)
```

The companion intersects every cell meeting the pulled-back rail cylinder
with the exact present union for that cell's own digit.  A separate referee
reconstructed the state translations and repeated this all-event
intersection independently, obtaining `(20)`--`(24)`.

## 6. Consequence and exact boundary

The selected lawful central cycle is now closed as a full present-bearing
chronology:

```text
arbitrarily long rail/clock germ
does not imply
arbitrarily long dynamically typed present packet.       (26)
```

The horizon `6` is structural: it matches the odometer depth `R=13^6` at
which the repelling rail cylinder resolves the missing future-digit sidecar.
This suggests testing alternative lawful lifts and heterogeneous edge words
by their depth-matched present gaps, not merely by finite rail support.

The result is sharply local.  It closes only the selected alternating rails
`(6,1,3,0)` and `(6,1,4,12)` and their plus-current present chronology.  It
does not close another lawful lift, another rail/component sequence, the
known mixed dilation/slope face, or arbitrary configuration switching.
Delayed words, private edge/root geometry on the surviving `H<=5` components,
semantic endpoint transport, all `165` rows, and LRC(14) remain open.

There is also an arithmetic reason not to extrapolate from this cycle.  Any
symmetric two-cycle `1/2+u,1/2-u` under opposite lifts of magnitude `K`
satisfies

```text
(13+1)u=K/R,            {R(1/2+u)}={1/2+K/14}.           (27)
```

Thus every such symmetric design is forced onto the speed-`14` resonance;
choosing a larger symmetric lift cannot remove it.  Asymmetric periodic
tails are the first honest alternative.

Run

```bash
python3 04-computation/lrc14_lawful_odometer_alternating_rail_horizon_probe.py
python3 -O 04-computation/lrc14_lawful_odometer_alternating_rail_horizon_probe.py
python3 04-computation/lrc14_lawful_odometer_alternating_unit_line_boundary.py
python3 -O 04-computation/lrc14_lawful_odometer_alternating_unit_line_boundary.py
```

Normal and optimized executions match the stored outputs byte-for-byte.
