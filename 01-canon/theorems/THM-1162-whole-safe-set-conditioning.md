---
id: THM-1162
title: Whole-safe-set conditioning scout, corrected by the two-run continuum model
status: CORRECTED FINITE-EXACT EVIDENCE plus a proved consecutive-limit profile.  The original single-run description is false.  Exact finite scans over all 495 cores on the named killer rows find a good gap, and the consecutive continuum bad locus has two mirror runs of total measure 2/21.  The all-spacing 2/21 ceiling and the finite-to-continuum/discrete-gap transfer remain open; uniform r=5 and LRC(14) remain open
source: kind-pasteur-S128c73, corrected by S128c74--c76 and codex-S75
depends_on:
  - THM-1167
  - THM-1146
  - THM-1147
  - THM-1148
related: [THM-1160, THM-1161, THM-1094, THM-1097]
script: 04-computation/bad_zone_width_kps_S128c73.py, 04-computation/whole_safe_set_kps_S128c73.py
---

# THM-1162 — corrected whole-safe-set conditioning scout

The useful global question is not whether every `k1`-gap is good.  THM-1167
refutes that certificate.  It is whether the core-safe set `S(P)` contains at
least one *complete* `k1`-gap whose three foreign teeth leave a component
longer than `1/(7k4)`.

## 1. What the finite exact bank establishes

For each of the following killer quadruples, the computation enumerates all
495 eight-speed cores and maximizes over the complete `k1`-gaps contained in
`S(P)`:

| killers | minimum over all cores of `7*k4*best` |
|---|---:|
| `(157,158,159,160)` | `1.97771` |
| `(197,198,199,200)` | `4.50000` |
| `(300,301,302,303)` | `4.45500` |
| `(317,318,319,320)` | `4.47161` |
| `(394,395,396,397)` | `4.45051` |
| `(371,374,377,379)` | `2.23010` |

A strided consecutive sweep through `k1=420` has minimum `1.74074` at
`(191,192,193,194)`.  These are exact finite certificates on the stated
rows, not a quantifier over all quadruples.

## 2. Correction: the bad set has two runs

The original scout misread “twelve bad indices, longest run six” as one run.
For consecutive offsets `(1,2,3)`, the normalized tooth positions are driven
by `{0,u,2u,3u}`.  Reflection `u -> 1-u` preserves that set, so the bad locus
has two mirror components, near `u=1/4` and `u=3/4`.

In the continuum model the first run is exactly

```text
[5/21,2/7]
```

and has width `1/21`; reflection gives total bad measure `2/21`.  This
corrects the former one-run/`0.0457` claim without invalidating the finite
all-core certificates above.

For general offset triples, the proportional family `m*(1,2,3)` has the same
total `2/21` by period rescaling.  A complete finite sweep over
`1<=d2<d3<d4<=16` found no larger value, but maximality of the proportional
family is not proved.

## 3. What a uniform proof still needs

The attractive inequality

```text
bad phase measure <=2/21 < 0.164 <= measure(S(P))
```

is not yet a complete discrete selection proof.  It needs three uniform
suppliers:

1. the all-offset continuum ceiling `bad<=2/21`;
2. a finite-scale error bound transferring that ceiling to actual teeth; and
3. an eroded core-atlas bound ensuring that a safe phase outside the bad set
   contains a *complete* `k1`-gap, rather than merely a safe point.

This is the corrected component-aware route.  It is substantially narrower
than the original six-linear-functions problem, but none of these three
suppliers is claimed here.  Uniform `r=5` remains open.

## 4. Carrier/tournament audit

The runner-order tournament is transitive and forgets the predicate.  The
proof-bearing carrier has vertices at gap indices (or wall-crossing events),
with the normalized phase `u`, three labelled tooth centres, core-safe
component ownership, and erosion width as sidecars.  It preserves “this
complete gap lies in `S(P)` and survives all three teeth”; quotienting to
tooth order alone destroys metric width and core availability.  This is the
assumption challenged by the two-run correction: bad phases, not runners,
are the natural vertices.
