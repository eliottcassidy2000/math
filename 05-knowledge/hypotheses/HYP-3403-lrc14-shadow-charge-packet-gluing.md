---
id: HYP-3403
title: LRC14 shadow-charge packet gluing
status: EVIDENCE / actual-packet controlled-forgetting test; not an LRC14 proof
source: codex-2026-06-28
tangent: T1364
technique: LTI-364
tournament_technique: LTT-264
script: 04-computation/lrc14_shadow_charge_packet_gluing_codex_20260628.py
result: 05-knowledge/results/lrc14_shadow_charge_packet_gluing_codex_20260628.out
reflection: 07-reflections/lrc14-shadow-charge-packet-gluing-codex-20260628.md
related:
  - HYP-3400
  - HYP-3401
  - HYP-3402
  - HYP-3311
  - HYP-3310
  - HYP-3301
  - HYP-3300
  - HYP-3266
  - HYP-3265
  - HYP-3260
  - HYP-3259
  - HYP-3258
  - HYP-3256
  - HYP-3254
  - HYP-3252
  - HYP-3249
  - HYP-3248
  - HYP-3247
  - HYP-3246
  - THM-523
  - OPEN-Q-108
---

# HYP-3403: LRC14 Shadow-Charge Packet Gluing

## Claim

The current LRC14 proof packet should glue three shadows, but not treat all
three as equally terminal:

```text
unit C3 / index shadow
quadratic Q(sqrt(-7)) sign shadow
nonunit covering residue + magnitude shadow
```

On the actual HYP-3311 theorem-facing bank, the ambient index/C3 and
`Q(sqrt(-7))` shadows are descriptive but not separating.  The first low-cost
sidecar that separates theorem exits is the nonunit covering residue word.
The `v2` magnitude word alone is too weak on this bank.  The full three-shadow
packet is a good gluing interface, while exact height is a debt detector for
same-residue moves rather than a proof quotient.

This is not a proof of LRC14.  It is a sharper rule for which creative
reframes are allowed to become proof coordinates.

## Exact Readout

The script reuses the `31` actual-packet rows from HYP-3311 and computes
controlled-forgetting fibers for theorem `kernel_flag`, route, and transfer
labels.

```text
kernel_flag_hist =
  q-witness-discharge      12
  AP/GW-zero-open-equality  2
  unit-petal-named          4
  named-K33-state-lift      3
  positive-Haar-open       10

q_bucket_hist =
  lt14 12
  eq14 12
  gt14  7
```

Sidecar table:

```text
sidecar                              fibers kernel route transfer maxK maxR maxFiber
raw_coarse_sheaf_base                     6      1     1        6    2    2       10
ambient_index_unit_c3                    12      1     1        9    2    2        6
imaginary_quadratic_balance               8      1     1        7    2    2        7
c3_plus_qsqrt_binding_packet             16      1     1        8    2    2        6
covering_residue_sheaf                   25      0     0        3    1    1        5
two_adic_magnitude_gate                  24      1     1        4    2    2        4
residue_magnitude_cover_signature        29      0     0        1    1    1        3
three_shadow_packet                      29      0     0        1    1    1        3
exact_height_hinge_oracle                31      0     0        0    1    1        1
analytic_lifting_ledger                  31      0     0        0    1    1        1
```

Thus:

```text
C3 + Q(sqrt(-7)) binding packet still leaves 1 mixed kernel fiber.
nonunit residue word leaves 0 mixed kernel fibers and 0 mixed route fibers.
nonunit v2 word alone still leaves 1 mixed kernel fiber.
```

The `qdiv>14` slice remains positive-open on this bank:

```text
qdiv_gt14_rows = 7
qdiv_gt14_kernel_hist = {positive-Haar-open: 7}
```

The same-residue height warning is still real.  There are `3`
same-residue height-debt fibers, led by the divisor-loaded tail family
`84, 84*2, 84*3, 84*4, 84*5`.  All are positive-Haar-open here, so this bank
does not show a theorem-exit failure, but it records exactly where an enlarged
bank can make residue-only data unsafe.

## Interpretation

The concurrent HYP-3401 three-coordinate exactness, HYP-3402 owner-current
tropical wall, S294 `Tropical_Current_Ledger` coordination, index-theorem,
imaginary-quadratic, C6, sheaf, and analytic-lifting threads now fit into a
stricter packet:

```text
index / C3:
  describes the unit binding saddle and its three antipodal slots.

Q(sqrt(-7)):
  organizes the quadratic sign split inside the unit and floor residue data.

covering residue:
  first actual sidecar that separates theorem exits on the current packet bank.

v2 / height:
  not first repair here, but necessary debt for same-residue enlarged-family
  moves and for the 12->24 hinge story.

transfer / analytic lifting:
  useful as a final ledger after residue and height sidecars are named, but
  too expensive to be the first proof quotient.

Tropical_Current_Ledger:
  downstream owner/current and wall-crossing audit for failures of
  residue-word exactness; not a replacement for naming the first low-cost
  separating sidecar.
```

So the creative reframe is:

```text
do not ask whether the index theorem or Q(sqrt(-7)) proves LRC14;
ask which packet coordinate each shadow preserves, which coordinate it
destroys, and whether the covering-residue repair plus height debt can be
globalized.
```

## Tournament Analysis

Vertices are proof sidecars and shadow-charge carriers, not runners or arcs.

Pairwise observable:

```text
(mixed kernel fibers, mixed route fibers, payload cost,
 mixed transfer fibers, preserved charge)
```

Switch/gauge:

```text
A -> B iff A has fewer mixed theorem-exit fibers, then lower payload debt.
```

Fingerprint:

```text
vertex_count = 10
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
priority_path =
  covering_residue_sheaf
  -> residue_magnitude_cover_signature
  -> three_shadow_packet
  -> analytic_lifting_ledger
  -> exact_height_hinge_oracle
  -> raw_coarse_sheaf_base
  -> imaginary_quadratic_balance
  -> ambient_index_unit_c3
  -> two_adic_magnitude_gate
  -> c3_plus_qsqrt_binding_packet
```

The path is not a beauty ranking.  It says the first proof-facing repair on
the actual bank is covering residue; the more elaborate shadows become
sidecars or debt once the first repair is retained.

## Assumption Challenge

Alternate vertex sets considered:

```text
runners
gaps
fixed circle sections
section boundaries
wall-crossing events
residues
cover arcs
Fourier/cyclotomic modes
matroid contact circuits
proof sidecars / shadow-charge carriers
```

Chosen vertices: proof sidecars.  The quotient preserves the LRC predicate
only as theorem-exit separability for LRC14 packet rows.  It destroys exact
height, endpoint owners, and off-grid floor unless those fields are retained,
discharged, or named as debt.

## Next Proof Obligations

A. Enlarge the HYP-3311 actual-packet bank toward HYP-2963 and record the
first row where nonunit residue no longer separates kernel or route.

B. If residue fails, classify the first missing repair as `v2`, exact height,
endpoint owner, transfer/state lift, or off-grid floor.

C. Keep C3/index and `Q(sqrt(-7))` as binding/residue shadows, but do not let
either replace the covering residue/height debt.

D. Prove a theorem interface:

```text
coarse sheaf base
+ covering residue
  separates theorem exits on the finite packet bank;
same-residue enlarged rows
  either remain positive-open, are repaired by v2/height/owner sidecars,
  or become named residual debt.
```

E. Re-run the tournament on every enlarged bank and require the first
theorem-facing repair to stay below the exact-height oracle.
