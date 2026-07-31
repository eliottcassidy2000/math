# The odometer escape has three events, then the delayed tail stops uniformly

**VERIFIED-EXACT SCOUT SYNTHESIS; NOT CANON AND NOT AN LRC(14) CONCLUSION.**

The first odometer reflection isolated the right missing coordinate but left
the decisive physical test open.  The envelope-level cycle with
`k=13^5+1` has correct intrinsic clock edges, yet one of its two exact points
misses the matching present factor.  More generally every symmetric exact
integer-lift cycle lands its delayed coordinate on the speed-14 danger grid,
so the exceptional points themselves are not full THM-2640 atoms.

The nearest present-side repair is

```text
R=13^6,              S=13^5,              k=S+2,
x_+=1/2+k/(14R),     x_-=1/2-k/(14R),
T_-(x)={13x-k/R},    T_+(x)={13x+k/R}.
```

Then `T_-(x_+)=x_-` and `T_+(x_-)=x_+`.  The inherited stored clock edges
are `4->3` and `3->4`, so owner-to-next-shallow gluing is literal.  Both
points meet eleven central rails and the correct present packet.  The exact
points still lie on a delayed boundary, but a nearby positive interval
supports three fully labelled event states.

For source label `s=1`, delayed sector zero, and the right half-tooth, one
such path is

| event | shallow -> owner | `(c,h,kappa)` | private root | intermediate `D`-root | unit determinant |
|---|---|---|---|---|---|
| 0 | `4 -> 3` | `(7,7,1)` | `2` | `2` | `2` |
| 1 | `3 -> 4` | `(5,9,1)` | `11` | `6` | `1` |
| 2 | `4 -> 3` | `(11,8,1)` | `10` | -- | `1` |

The source interval is exactly

```text
[4286847956371849/8480502984583084,
 4286847956371861/8480502984583084),
```

of width `3/2120125746145771`.  Each event has ten unit source rails for
each half-edge.  Allowing independent sources gives `2*10^3=2000` labelled
three-event paths; imposing one common source leaves `16`.

The root typing is a useful correction to the tempting additive shortcut.
The THM-2657 increment acts on the root after dilation, not directly on the
current private root:

```text
r=2  --D-root 2 --(+9)--> r=11,
r=11 --D-root 6 --(+4)--> r=10.
```

The second line is the visible H-drift.  Forgetting the intermediate
`D`-root would incorrectly predict `11+4=2`.

## The uniform stopping coordinate

Put `y={Rx}`.  For every integer `m`, not merely `m=+/-k`,

```text
R*{13x+m/R} = 13Rx+m        modulo integers,
```

and therefore

```text
y({13x+m/R})={13y}.                                  (1)
```

The odometer lift changes the predecessor carry by `m mod 13` and the root
by `2m mod 13`, but it is invisible in (1).  Thus the correct object is a
skew product: the delayed base evolves autonomously by multiplication by
thirteen, while the odometer integer is a cocycle sidecar.

Let `W_0,W_1` be the two raw delayed guard-sector words before any clock cut,
and put `W=W_0 union W_1`.  Exact half-open interval arithmetic gives

| event states `d` | components of `intersection_(j<d) D^(-j)W` | exact mass |
|---:|---:|---:|
| 1 | 47,484 | `604725613249/11455265301480` |
| 2 | 16,244 | `513351/371664293` |
| 3 | 6,776 | `2662/62748517` |
| 4 | 0 | `0` |

The large census conceals a two-factor proof.  The safe and danger guard
sectors partition the guard circle exactly, so their union deletes the guard
condition:

```text
W = D_(13^3)
    intersection D_14^c intersection D_27^c intersection D_40^c
    intersection D_53^c intersection D_66^c
    intersection D_(2*13^5)^c.
```

In fact only `D_(13^3)` and `D_14^c` are needed.  Use the centered half-open
danger interval `I=[-1/14,1/14)` and put `z={13^3 y}` in centered
coordinates.  Membership in the target tooth at event states zero, one, and
two says

```text
z in I,       13z mod 1 in I,       13^2 z mod 1 in I.
```

There is no wrap branch.  Indeed `13/14=1-1/14`; a nonzero wrap would require
`|z|>=1/14`, and the only boundary candidate maps to the excluded right
endpoint of `I`.  Applying this twice gives

```text
z in (1/13^2) I.
```

At event state three the speed-14 unit factor is evaluated at `13^3 y`, so
its argument is `14z`.  Since `14<13^2`, the preceding inclusion forces
`14z in I`, contradicting the required speed-14 safe factor.  Thus the
depth-four stop does **not** use the fourth target tooth, the other four unit
speeds, the `2*13^5` tooth, a guard sector, a clock label, or an odometer
lift.  The exact two-factor superset has component/mass rows

| event states `d` | components | exact mass |
|---:|---:|---:|
| 1 | 1,886 | `6/49` |
| 2 | 1,606 | `22187/2798978` |
| 3 | 1,452 | `1452/2599051` |
| 4 | 0 | `0` |

This is the underlying toothpick mechanism: three consecutive target teeth
contract the centered `13^3` phase by `13^2`, and the fourth-time image of
the smallest unit speed, `13+1`, fits strictly back inside the forbidden
tooth.  The huge interval automaton is a refinement of this scale mismatch,
not its cause.

The same argument is an elementary all-parameter lemma.  For `p>=2` and
`r>=2`, let `D(y)={p y}`, let
`I=[-1/(p+1),1/(p+1))`, and consider any carrier contained in

```text
M_(p,r) = {y : p^r y mod 1 is in I,
                (p+1)y mod 1 is not in I}.
```

Then

```text
intersection_(j=0)^r D^(-j) M_(p,r) = empty.             (2)
```

Indeed the target conditions at `j=0,1,2` contract
`z={p^r y}` to `z in I/p^2`, while the safe condition at `j=r` asks
`(p+1)z` to avoid `I`; this is impossible because `p+1<p^2`.  Equation (2)
is a scale-gap stopping lemma.  The present exact result is its sharp
`(p,r)=(13,3)` instance: the depth-three intersection is positive, so the
fourth state is the first forced zero for this minimal carrier.

Every binary sector word of lengths one, two, and three is positive.  Every
one of the sixteen length-four sector words is empty.  Consequently the
refined component transition graph has nilpotence index four and no recurrent
strongly connected component.  A sector-only two-vertex graph would conceal
this: the load-bearing state is the interval component (or an equivalent
length-three stalk), not merely the sector bit.

This is stronger than an alternating-clock accident.  Every clock-cut word
`Q_(sector,ell)` is a subset of `W_sector`; hence the raw depth-four zero
holds uniformly for every intrinsic-clock word and every sequence of integer
odometer lifts.  On the witness itinerary `4,3,4,3`, the component counts are
`41,217`, `12,610`, `5,214`, `0`, with depth-three mass
`28677/878479238`.

The positive three-event interval therefore genuinely escapes the fixed-`D`
clock diagonal, but the inherited delayed packet stops it before a fourth
event.  This is not a row exclusion: it says that one proposed chronology
cannot iterate inside these inherited packets.  Endpoint transport, scalar
closure, and an LRC(14) exit remain absent.

Exact replay:

```bash
python3 04-computation/lrc14_odometer_alternating_lift_labelled_tail_scout_20260728.py
python3 -O 04-computation/lrc14_odometer_alternating_lift_labelled_tail_scout_20260728.py
```
