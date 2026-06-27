---
id: HYP-2967
title: LRC14 apex-aperture comb trichotomy
status: PROOF-INTERFACE / exact local certificate plus residual trichotomy; not a proof
source: codex-2026-06-24-S154
script: 04-computation/lrc14_apex_aperture_comb_certifier_codex_s154.py
result: 05-knowledge/results/lrc14_apex_aperture_comb_certifier_codex_s154.out
related:
  - HYP-2962
  - HYP-2961
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2960
  - HYP-2955
  - HYP-2953
  - HYP-2910
  - HYP-2909
  - HYP-2908
  - HYP-2905
  - THM-523
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-2967: LRC14 Apex-Aperture Comb Trichotomy

This pass attacks the actual residual after THM-571.  The many-multiples branch
is already closed:

```text
|M14(S)| >= 7  =>  M(S) > 1/14.
```

So the live covering-strictness branch has at most six multiples of `14`.  Write

```text
S = C union 14Q,
C contains no multiple of 14,
1 <= |Q| <= 6.
```

Rebased over the incoming HYP-2964 Moon-core skeleton, HYP-2965 boundary-gap
packet bridge, and HYP-2966 NORK pinch-template atlas, HYP-2967 should be read
as a local reducer inside the Moon core: it certifies rows with a surviving
denominator-14 aperture and names the two local obstruction types that still
need off-apex packet machinery.

## Apex-Aperture Lemma

Fix a unit `a in (Z/14Z)^*`.  At `t=a/14`, every `c in C` is closed-safe:

```text
||c a/14|| >= 1/14.
```

Let `[ac]_14` be the representative in `{1,...,13}`.  For the right side
`t=a/14+h`, a strict `C`-safe aperture exists until

```text
U_a^+(C) = min_c (13 - [ac]_14)/(14c).
```

If some `[ac]_14=13`, then `U_a^+=0`.  Similarly the left aperture is

```text
U_a^-(C) = min_c ([ac]_14 - 1)/(14c),
```

with `[ac]_14=1` as the left-side blocker.

For the `14Q` speeds, the danger comb in aperture coordinate `h` is

```text
D_Q = union_{q in Q, k in Z}
      [ (k-1/14)/(14q), (k+1/14)/(14q) ].
```

If either one-sided interval `(0,U_a^+)` or `(0,U_a^-)` is not covered by
`D_Q`, then `S` has a strict lonely time at threshold `1/14`: choose any
rational `h` in the uncovered gap and set `t=a/14 +/- h`.

## Trichotomy

After q-witness discharge and THM-571, every remaining strict counterexample
must be in one of these local-apex obstruction classes.

`Aperture-certified`.
: Some denominator-14 unit has a one-sided core aperture not covered by the
  `14Q` combs.  Then `M(S)>1/14`, constructively.

`Full unit-support / AP-GW-balanced core`.
: Every one-sided aperture is first-order blocked.  Equivalently, the core `C`
  hits every unit residue class modulo `14`.  This is the AP/GW-like local
  skeleton and should route through HYP-2960, source-spectrum labels, or
  HYP-2908/THM-572 state lift.

`Comb-saturated tiny aperture`.
: At least one one-sided aperture exists, but every such aperture is covered by
  the danger combs from at most six `14`-multiples.  This gives explicit
  rational inequalities between the smallest `q in Q` and the boundary speed of
  `C`, so the next theorem should force scale separation or a bounded atlas.

## Computation

Script:

```text
04-computation/lrc14_apex_aperture_comb_certifier_codex_s154.py
```

Stored output:

```text
05-knowledge/results/lrc14_apex_aperture_comb_certifier_codex_s154.out
```

Global S154 readout over named rows plus the S151 AP-neighborhood banks:

```text
total_rows: 68375
live_low_m14_rows: 18909

aperture-comb-certified:           12548
all-apertures-first-order-blocked:  4661
all-apertures-comb-saturated:       1700
```

The local certificate proves about two thirds of the low-multiple live rows by
explicit rational witnesses.  The remainder are exactly the two obstruction
types above.

## Negative Signal

The named covering rows

```text
12->84,
12->168,
drop(4,6)->add(19,42),
drop(2,6)->add(17,42)
```

are first-order blocked at every denominator-14 apex.  S151 knows they are
positive-open, so their witness is off-apex.  A proof of LRC14 cannot rely only
on one-sided denominator-14 local apertures.

## Tournament Analysis

For each row, S154 forms a six-vertex tournament on the denominator-14 unit
apexes.  The pairwise observable is the best uncovered aperture measure after
removing the `14Q` danger combs.  Continuous data is made discrete by exact
rational comparison; ties follow

```text
1 > 3 > 5 > 9 > 11 > 13.
```

Every row lands in the same six-vertex transitive isomorphism class:

```text
transitive6:000000000000000
```

That says this quotient is a proof gate, not a final classifier.  The proof-gate
tournament is also transitive:

```text
qdiv-direct
> THM-571-apex-majority
> aperture-comb-certified
> all-apertures-first-order-blocked
> all-apertures-comb-saturated
> state-lift-or-bounded-core-residual.
```

## Proof Target

The next theorem is:

> A comb-saturated tiny aperture with `1<=|Q|<=6` is scale-separated or belongs
> to a bounded finite-core atlas; a full unit-support core is AP/GW-skeleton and
> routes through HYP-2960/HYP-2908 rather than through a new analytic family.

Combined with THM-571, this would shrink covering-strictness to a finite
AP/GW-skeleton/state-lift problem.
