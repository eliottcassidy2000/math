---
id: THM-2693
title: "Odometer skew-product three-event escape and uniform delayed depth-four nilpotence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY REPLAYED.  On the canonical
  THM-2584/THM-2640 carrier, the lawful affine lifts +/-(13^5+2) have a
  positive fully labelled three-event interval with literal clock gluing,
  dynamically typed present factors, delayed sector words, predecessor
  carries, private roots, and primitive units.  Nevertheless the terminal
  coordinate y={13^6 x} evolves as y'={13y} for every integer lift.  The
  union of the two raw delayed sectors is positive through three event states
  and empty at four, uniformly over all lift sequences and all clock words.
  The zero already follows from three consecutive D_(13^3) target teeth and
  the fourth-state D_14 safety factor.  This closes the inherited delayed
  carrier, not other delayed grammars, semantic endpoints, rows, or LRC(14).
source: root-long-frontiers/lrc-delayed-skew-product-2026-07-28
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2689-affine-clock-support-typing-tradeoff-and-odometer-phase-locality
related:
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
  - THM-2691-lawful-odometer-alternating-rail-horizons-and-depth-six-present-collapse
script: 04-computation/lrc14_odometer_alternating_lift_labelled_tail_scout_20260728.py
output: 05-knowledge/results/lrc14_odometer_alternating_lift_labelled_tail_scout_20260728.out
script_sha256: 1efcac15f6fbf4cdd14fd909a0ac9ddf07eb34f70a6c16684049b209ed8fc482
output_sha256: 6a930b1093b235a6f5725c8ae6f157a6c5706bce449a5dd079da809ef6d6a365
hash_basis: working-tree bytes (LF)
---

# THM-2693 -- the odometer lift escapes for three events, but cannot steer the delayed tail

**PROVED + VERIFIED-EXACT + INDEPENDENTLY REPLAYED.**

THM-2689 proves that every nontrivial odometer lift restores intrinsic
three-event rail/clock support, while no such lift gives a global seven-clock
action.  THM-2691 follows one central lift through its present factors and
finds a sharp sixth-layer zero.  The missing coordinate is the delayed tail.
It is autonomous: integer carry labels change its odometer fibre but do not
change its base dynamics.

This theorem gives both necessary controls.  A nearby lawful lift supports a
strict positive three-event packet after all current local factors are
inserted.  Yet the inherited delayed language has nilpotence index exactly
four, uniformly in every integer lift, sector choice, and clock word.

## 1. A fully labelled positive three-event escape

Put

```text
p=13,              R=p^6,              S=p^5,
k=S+2,
T_-(x)={p x-k/R},  T_+(x)={p x+k/R}.                     (1)
```

The exact symmetric points

```text
x_+=17079479/33787663,       x_-=16708184/33787663       (2)
```

form a two-cycle under `T_-`, then `T_+`.  Their inherited intrinsic clock
edges are respectively `4->3` and `3->4`, so owner-to-next-shallow gluing is
literal.  The points themselves lie on a delayed wall, but a neighbouring
open interval supports two handoffs and three fully typed event states.

For common source `1`, delayed sector `0`, and the right half-tooth, take

```text
J=[4286847956371849/8480502984583084,
   4286847956371861/8480502984583084).                    (3)
```

Its length is

```text
3/(7*13^13)>0.                                            (4)
```

Every `x in J` has the following state data:

| event | clock edge | `(carry,h,kappa)` | private root | intermediate `D`-root | unit determinant |
|---:|---:|---:|---:|---:|---:|
| 0 | `4->3` | `(7,7,1)` | `2` | `2` | `2` |
| 1 | `3->4` | `(5,9,1)` | `11` | `6` | `1` |
| 2 | `4->3` | `(11,8,1)` | `10` | -- | `1` |

All three states retain eleven active source rails, ten primitive-unit rails,
the dynamically typed present packet, delayed word, predecessor carry,
future half-digit, private root, and unit test.  The carry handoffs are

```text
7 --(h=7,-2)--> 5,          5 --(h=9,+2)--> 11.          (5)
```

The root law must retain the intermediate dilation root:

```text
2  --D-root 2,+9--> 11,
11 --D-root 6,+4--> 10.                                  (6)
```

Adding `9` or `4` directly to the current private root would be ill typed.
There are `16` same-source labelled paths on this interval and `2,000` if the
three sources are allowed to vary independently.  Thus the eventual zero is
not inherited from the unshifted support failure of THM-2682/2684.

## 2. The affine map is a skew product

For an arbitrary integer lift `m`, write

```text
Rx=N+y,        N in Z/RZ,        0<=y<1,
T_m(x)={13x+m/R}.                                      (7)
```

If `a=floor(13y)`, then exactly

```text
y(T_m x)={13y},
N(T_m x)=13N+m+a                 (mod R).                (8)
```

The integer lift acts only in the odometer fibre.  In particular, no choice
of a time-dependent sequence `(m_0,m_1,...)` can steer the delayed base away
from the map

```text
B(y)={13y}.                                               (9)
```

This is the type separation that a carry-only or clock-only automaton loses.

## 3. The raw delayed language dies at four states

Use the canonical half-open danger tooth

```text
I=[-1/14,1/14)
```

in centered circle coordinates, and write

```text
D_q={y : centered(qy) in I}.                              (10)
```

Let `W_0,W_1` be THM-2623's two raw delayed guard sectors before any clock
cut.  Their guard factors are complementary, so their union is exactly

```text
W=W_0 union W_1
 =D_(13^3) intersect D_14^c intersect D_27^c
              intersect D_40^c intersect D_53^c
              intersect D_66^c intersect D_(2*13^5)^c.  (11)
```

Exact half-open interval arithmetic gives

| event states `d` | components of `intersection_(j<d) B^(-j)W` | exact mass |
|---:|---:|---:|
| 1 | `47,484` | `604725613249/11455265301480` |
| 2 | `16,244` | `513351/371664293` |
| 3 | `6,776` | `2662/62748517` |
| 4 | `0` | `0` |

Every one of the `2^d` sector words is positive for `d=1,2,3`; none of the
sixteen sector words at `d=4` is positive.  After the alternating clock cuts
`4,3,4,3`, the component counts are `41,217`, `12,610`, `5,214`, `0`, with
depth-three mass `28677/878479238`.  Hence the refined component graph has
nilpotence index four and no recurrent strongly connected component.

Every clock-cut word is a subset of one raw sector.  Combining this inclusion
with `(8)` proves the uniform statement:

```text
every four-state inherited delayed packet is empty
for every integer-lift sequence and every intrinsic-clock word. (12)
```

## 4. The two-factor analytic certificate

The large census in fact uses far too much data.  Retain only

```text
M=D_(13^3) intersect D_14^c.                              (13)
```

At event states zero, one, and two, put `z=centered(13^3y)`.  The target
conditions say

```text
z in I,          centered(13z) in I,
                 centered(13^2z) in I.                   (14)
```

There is no nonzero wrap branch.  Indeed

```text
13/14=1-1/14,                                        (15)
```

and either apparent wrap would require the excluded right endpoint of `I`.
Applying this observation twice yields

```text
z in I/13^2.                                             (16)
```

At event state three the speed-14 safety factor is evaluated at `13^3y`, so
its centered phase is `14z`.  Since `14<13^2`, equation `(16)` gives
`14z in I`, contradicting safety.  Thus the zero does not use the fourth
target tooth, the other four ordinary unit speeds, the high-speed tooth, a
guard sector, a clock label, or an odometer lift.

As an independent positive/hostile control, the repeated two-factor carrier
`M` has component/mass rows

```text
d=1: (1886, 6/49),
d=2: (1606, 22187/2798978),
d=3: (1452, 1452/2599051),
d=4: (0,0).                                               (17)
```

## 5. General sparse nested-danger lemma

The scale mechanism is not special to `13`.  Let `p>=2`, `r>=2`, and put

```text
I_p=[-1/(p+1),1/(p+1)),       B_p(y)={py}.                (18)
```

For positive integers `a,q` with

```text
a divides q,                 q/a <= p^2,                 (19)
```

one has the factor-sparse identity

```text
D_(a p^r)
 intersect B_p^(-1)D_(a p^r)
 intersect B_p^(-2)D_(a p^r)
 intersect B_p^(-r)D_q^c = empty.                        (20)
```

To prove it, write `z=centered(a p^r y)`.  The first three factors give
`z,pz,p^2z in I_p`.  Since

```text
p/(p+1)=1-1/(p+1),                                      (21)
```

the same endpoint argument as `(15)` removes both wrap branches and gives
`z in I_p/p^2`.  At state `r`, the congruence

```text
a B_p^r(y) = z                  (mod 1)                  (22)
```

and divisibility `a|q` give

```text
q B_p^r(y)=(q/a)z               (mod 1).                 (23)
```

Condition `(19)` puts `(23)` back in `I_p`, contradicting the last factor.

The coefficient threshold in `(19)` is sharp for precisely these four
retained factors.  If `c=q/a>p^2`, choose

```text
1/((p+1)c) < z
 < min(1/((p+1)p^2),1/(2c)),                            (24)
```

and set `y=z/(a p^r)`.  Then the three target factors in `(20)` hold while
`c z` is safely outside `I_p`.  This does not assert survival of any omitted
factor; it identifies the exact boundary of the analytic implication.

The theorem's depth-four stop is `(p,r,a,q)=(13,3,1,14)`.  The same lemma at
`(13,5,2,14)` gives a factor-sparse depth-six obstruction to any fixed word
retaining the target tooth `D_(2*13^5)` and the speed-14 safety factor.  This
explains why the two canonical blocker scales exhibit four- and six-state
stops without identifying their remaining sidecars.

## 6. Exact boundary and reproduction

The proved obstruction closes the inherited delayed carrier.  It does not
close an affine handoff with a different delayed grammar, a phase-shifted or
reversed speed-14 factor, a heterogeneous `D`/slope word, or a carrier that
changes the base dynamics.  The positive interval has no semantic endpoint
current or scalar return.  No row is excluded, and the LRC(14) ledger remains
`165`.

Run

```bash
python3 04-computation/lrc14_odometer_alternating_lift_labelled_tail_scout_20260728.py
python3 -O 04-computation/lrc14_odometer_alternating_lift_labelled_tail_scout_20260728.py
```

The Fraction cylinder recursion and the independent refined-integer-grid
intersection agree internally.  Normal and optimized runs byte-match the
stored output.  An independent replay reproduced the output hash
`6a930b1093b235a6f5725c8ae6f157a6c5706bce449a5dd079da809ef6d6a365`.

QED.
