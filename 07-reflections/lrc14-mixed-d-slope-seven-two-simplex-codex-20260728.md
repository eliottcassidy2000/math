# D-nilpotence is not mixed nilpotence

**Status: FINITE-EXACT POSITIVE SCOUT + INDEPENDENT DIRECT REPLAY.**  This is
an unnumbered research bridge, not proved canon and not an LRC(14) closure.

THM-2680 proves that a physical two-event dilation fibre product exists.
THM-2682/2684 prove that its inherited `D(x)={13x}` grammar is nilpotent at
three events.  THM-2672, meanwhile, proves honest positive overlaps between
the slope-seven THM-2640 carry charts.  The cheapest faithful cross-test is
therefore not a third `D` event.  It is one `D` edge followed by one physical
carry/root configuration switch.

## 1. The heterogeneous correspondence

Put

```text
R=13^6,                    D(x)={13x},
tau_delta=7 delta/R.
```

Let `E_0,E_1` be THM-2680 atoms with the typed clock/digit covariance, and
let `P_(e,c)` be a full unintegrated THM-2640 packet.  At `z=D(x)` define

```text
M_delta = E_0(x)
          intersect E_1(z)
          intersect P_(e,c)(z)
          intersect P_(e,c+7delta)(z+tau_delta).          (1)
```

Every factor in (1) is a physical nonnegative packet before integration.
The last factor is exactly the slope-seven carry/root lift:

```text
c -> c+7delta,                 r -> r+delta.              (2)
```

Thus (1) is a mixed component-refined two-simplex: its first edge is the
proved dilation handoff and its second edge is a translated configuration
correspondence.  It is not a three-event `D` chain.

## 2. Exact positive witness

Take THM-2680's first displayed positive pair:

```text
clock path                 3 -> 1 -> 0,
current atom               (source,rail,j,h,eps,kappa)
                           =(1,2,5,2,0,0),
following atom             =(1,0,2,6,1,1).               (3)
```

The first exact open component found by the primary route contains

```text
x =649039434905733/1304692766858936,
z ={13x}=46873542509301/100360982066072.                 (4)
```

At `z`, retain the THM-2640 configuration

```text
(rail,sector,edge,local-kappa,h,shallow)=(0,0,0,1,6,1),
(carry,root)=(2,6).                                      (5)
```

Both (3) and (5) hold strictly at the rational point (4), including the
rail, present factor, safe-sector delayed word, predecessor carry, future
half-digit, private half-tooth, and primitive unit test.  For every

```text
delta in {1,2,3,4,5,6,9,10,11,12},                      (6)
```

the shifted point `z+7delta/R` lies strictly in the corresponding unit
packet with labels `(c+7delta,r+delta)`.  In particular,

```text
delta=1:             (2,6) -> (9,7),                     (7)
```

so `M_1` has a positive open component.  The two excluded labels have exact
and different mechanisms:

```text
delta=7:             shifted root is zero;
delta=8:             shifted row is not a THM-2640 unit. (8)
```

The primary script reconstructs the D-edge component and discovers (4)--(8).
The second script freezes (4) and independently checks every factor directly,
then exhausts the twelve nonzero slope labels.  Both normal and optimized
executions match the stored transcripts.

## 3. Why this evades the clock trap

THM-2682's obstruction compares the shallow clocks of `z=Dx` and `Dz` on a
central-arrival return cylinder.  The fixed-side identity forces those two
clocks equal, while a stored D-edge requires them to differ.

Equation (1) never asks for `Dz`.  Its second edge keeps shallow clock `1`
and changes the predecessor carry/private root by (2).  The central-return
identity is therefore irrelevant.  This is a genuine support bypass, but it
works by changing the edge grammar, exactly one of THM-2682's allowed escape
routes.  It does not refute or weaken D-chain nilpotence.

The stronger convenient guess also failed its cheap probe: none of the
thirteen selected *first-component midpoints* of THM-2672's canonical
twelve-chart facet bank refined the displayed D edge.  This is only a
midpoint spot-check, not an exclusion of all components or configurations.

## 4. Connection contract

```text
source:
  one positive THM-2680 D-edge component plus THM-2640 slope charts

target:
  a component-refined heterogeneous nerve with one D edge and one
  carry/root-switch edge

map:
  x -> (x, z={13x}, z+7delta/13^6)

preserved:
  one physical common point; D owner/shallow and digit covariance;
  sources, rails, carry, future half-digit, private root, and unit flags

not supplied:
  THM-2365 target action; canonical owner/endpoint current; a second D-like
  chronology edge; an order-thirteen action; terminal phase transport

needed sidecar:
  identify the translated private packet as a lawful next-event endpoint
  carrier, or replace scalar D by a labelled affine handoff and retain its
  odometer digit under composition
```

The incoming affine-cocycle scout
`lrc14-the-odometer-digit-is-a-handoff-coordinate-20260728.md` gives the
sharp next test.  Its phase-level maps with lifts `k=-14,+14` alternate the
central clock forever at a distinguished two-cycle.  The next decisive
physical computation is the corresponding labelled rail product with every
source, carry, half-edge, unit, and endpoint factor retained.  Our witness
shows that the first mixed physical face exists; it does not prove that the
affine faces concatenate.

## 5. Reproduction

```bash
python 04-computation/lrc14_mixed_d_slope7_two_simplex_scout.py
python -O 04-computation/lrc14_mixed_d_slope7_two_simplex_scout.py
python 04-computation/lrc14_mixed_d_slope7_two_simplex_independent_referee.py
python -O 04-computation/lrc14_mixed_d_slope7_two_simplex_independent_referee.py
```

The matching outputs are

```text
05-knowledge/results/lrc14_mixed_d_slope7_two_simplex_scout.out
05-knowledge/results/lrc14_mixed_d_slope7_two_simplex_independent_referee.out.
```

Their LF-normalized SHA-256 hashes are

```text
primary script   48113885218444c724d5d89588c467c113fbed40c362456085da235d13069f77
primary output   939e6e4fa307bb2f52242211bfe0ac9747492178066cdb57fbca888ddc189910
referee script   1215b285ebf5eb026facf2e63c595b90a75d9f0a733d2eae99d8066193f74d3a
referee output   e87a211d1eb1c981ecd57fc2492250fd5dc57bec39ff9934941bacf2c918887d
```

No scalar row exclusion, endpoint transition, holonomy trivialization, or
LRC(14) conclusion follows.
