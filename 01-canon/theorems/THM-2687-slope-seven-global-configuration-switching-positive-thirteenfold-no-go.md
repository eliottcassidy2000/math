---
id: THM-2687
title: "Slope-seven global configuration-switching positive thirteenfold no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the full
  THM-2640 slope-seven private-root carrier, allow every target label to choose
  its rail, base cell, present factor, delayed clock/sector, and private edge
  independently.  No positive open common component supports all thirteen
  labels.  A necessary fine-coordinate envelope forgets compatibility among
  those choices and therefore only enlarges every physical chart, yet covers
  at most twelve of thirteen carry values on every elementary future-coordinate
  cell.  The exact open-cell census is 47,512/39,948/7,536 at coverage
  0/11/12; the maximum cells miss only carry 0, 6, or 12.  Since
  delta -> c0+7delta permutes F_13, arbitrary labelwise configuration switching
  cannot evade the cap.  Isolated half-open boundary contacts are not
  classified set-theoretically and no LRC(14) row is excluded.
source: root/lrc-physical-gluing-probe-2026-07-28
depends_on:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
related:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2684-three-tooth-rail-envelope-diagonal-arrival-law-and-full-dilation-nilpotence
script: 04-computation/lrc14_slope7_global_configuration_switching_probe.py
output: 05-knowledge/results/lrc14_slope7_global_configuration_switching_probe.out
script_sha256: 851a5d0ac866b38a775209c8b8d97c3e1fed0dade3092a2e5021fa2e4edf8ca1
output_sha256: b5778f7e7dfba026f3e7a0fcabaca0f8a54b0c219b1ff4c30516cf621d405ca3
hash_basis: LF-normalized bytes
---

# THM-2687 -- configuration switching still misses one carry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2672 proves a root-zero cap only after fixing
`(rail,sector,edge,kappa,h)`.  It therefore leaves a precise loophole: choose
those data separately for each target label and ask whether the thirteen
resulting unions have a positive common component.  This theorem closes that
loophole inside the complete THM-2640 slope-seven carrier.

The proof is intentionally by an enlarged necessary envelope.  It forgets
several physical coordinates and all compatibility among their witnesses.  A
negative result for that superset is consequently stronger than a direct scan
of one gauge, while the converse is neither used nor asserted.

## 1. Common coordinate and carry orbit

Use the THM-2640 notation

```text
p=13,                   R=13^6,                   y={Rx}.       (1)
```

On the interior of a future half-digit, `y` fixes the pair

```text
(h,kappa) in F_13 x F_2.                                  (2)
```

For a source carry `c0` and target label `delta`, the slope-seven lift has

```text
c_delta=c0+7delta mod 13.                                 (3)
```

Multiplication by seven permutes `F_13`.  Thus thirteen target labels at one
common `y` require all thirteen carry values, in some order.

## 2. The enlarged fine-coordinate envelope

For each fixed `(y,h,kappa,c)`, retain only the following necessary tests:

1. some delayed clock in one of the two sectors has positive prefix support;
2. one of the two literal private half-edges contains `y`, with its exact root;
3. some one of the `162` THM-2640 rails has the exact unit flag for that
   `(sector,edge,c,kappa,h)`.

Existentially union over the delayed clock, sector, edge, rail, and base cell.
Forget the present factor, the component identity, and every compatibility
condition saying that the separate existential witnesses come from one
physical packet.  Denote the resulting option set by `E_c(y)`.

Every physical label chart satisfies

```text
physical option for carry c at y  =>  E_c(y) is nonempty. (4)
```

Every deletion used to define `E_c` enlarges support.  In particular, (4) has
the direction needed for a no-go; no converse is claimed.

The private-edge test has a two-threshold closed form.  Write

```text
d=2h+kappa=13b+a,             y=(d+v)/26,   0<v<1.       (5)
```

Edge one is available throughout the half-digit except when `a=12`, where it
requires `v<1/14`.  Edge zero is available throughout except when `a=0`,
where it requires `v>13/14`.  A literal circular-tooth reconstruction agrees
with this formula in all `780` carry/edge test cases.

## 3. Exact open-cell cap

Partition the future circle at every delayed-prefix endpoint, future
half-digit endpoint, and private-edge threshold from (5).  All option sets
`E_c(y)` are constant on each resulting open elementary interval.  Exact
enumeration gives

```text
number of covered carries       elementary intervals
             0                         47,512
            11                         39,948
            12                          7,536
            13                              0.            (6)
```

The corresponding exact `y`-measures are

```text
coverage 0:   10850539688231 / 11455265301480,
coverage 11:     28263039347 /   636403627860,
coverage 12:       258735593 /    30876725880.             (7)
```

They sum to one.  Every maximum cell has a unique missing carry, with census

```text
missing 0: 2067,              missing 6: 3402,
missing 12: 2067.                                        (8)
```

After rebasing by `c0`, every target label is the unique missing label in
exactly `7,536` maximum-cell instances.  Hence the asymmetry in (8) cannot be
removed by changing the source carry.

Suppose a physical thirteen-label intersection contained a positive open
component.  By (4), every carry option `E_c` would be nonempty throughout
that component.  Removing the finite partition walls leaves a nonempty open
subinterval contained in one elementary cell, contradicting the cap twelve
in (6).  This proves the theorem.

## 4. Positive and hostile controls

The canonical physical twelve-chart midpoint from THM-2672 lies strictly in
one maximum cell with `(h,kappa)=(6,0)` and missing carry six.  Thus the cap is
sharp and is not an artifact of an empty relaxed envelope.

The theorem concerns positive open components, exactly the support notion in
THM-2658's physical nerve.  It does not classify isolated half-open boundary
contacts.  Such points cannot create a positive component: any positive
intersection crossing a finite wall contains an adjacent open subinterval.

Run

```bash
python3 04-computation/lrc14_slope7_global_configuration_switching_probe.py
python3 -O 04-computation/lrc14_slope7_global_configuration_switching_probe.py
```

Both executions byte-match the stored output and the declared hashes.  An
independent referee rebuilt the superset implication, the carry permutation,
the open-wall argument, and the normal/optimized replay.  The audit confirmed
that forgotten rail, base, present, clock, sector, edge, and component
compatibilities only enlarge support.

This closes THM-2672's arbitrary configuration-switching thirteenfold question
inside the THM-2640 slope-seven private-root family.  It does not exclude any
of the `165` LRC(14) rows, provide a target transition, or close LRC(14).

QED.
