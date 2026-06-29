---
id: HYP-3475
title: LRC14 colored gate mirror-orbit audit
status: EVIDENCE / exact mirror-orbit quotient audit; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3474, HYP-3473, HYP-3472, HYP-3471, HYP-3461, and HYP-3438
tangent: T1435
technique: LTI-435
tournament_technique: LTT-335
script: 04-computation/lrc14_colored_gate_mirror_orbit_codex_20260629.py
result: 05-knowledge/results/lrc14_colored_gate_mirror_orbit_codex_20260629.out
reflection: 07-reflections/lrc14-colored-gate-mirror-orbit-codex-20260629.md
related:
  - HYP-3474
  - HYP-3473
  - HYP-3472
  - HYP-3471
  - HYP-3470
  - HYP-3462
  - HYP-3461
  - HYP-3460
  - HYP-3459
  - HYP-3458
  - HYP-3457
  - HYP-3456
  - HYP-3455
  - HYP-3454
  - HYP-3453
  - HYP-3451
  - HYP-3450
  - HYP-3438
  - HYP-3436
  - HYP-3425
  - HYP-2595
  - HYP-2594
  - HYP-2593
  - THM-523
  - OPEN-Q-108
---

# HYP-3475: LRC14 Colored Gate Mirror-Orbit Audit

## Claim

HYP-3471 turned low-rank survivor gates into typed endpoint-color words.  The
incoming HYP-3472/HYP-3473/HYP-3474 stack then added the boundary-current
sidecar, Lean producer interface, and quotient-legality partition lattice.
The missing extension coordinate from HYP-3461 is the mirror orbit:

```text
gate interval I  <->  1-I
branch0          <->  branch1
left/right debt  <->  right/left debt
```

On the current `135`-row HYP-3438/HYP-3453 bank, every survivor gate belongs to
an exact two-gate mirror orbit.  There are no fixed gates and no unpaired gates.
This is the right finite carrier for colored gate-extension rows: it retains
typed endpoint residues, branch masks, adjacency, cover-delta vectors, and the
mirror partner before any scalar gate count is used.

## Exact Readout

Executable scout:

```text
04-computation/lrc14_colored_gate_mirror_orbit_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_colored_gate_mirror_orbit_codex_20260629.out
```

Aggregate mirror ledger:

```text
rows_audited=135
survivor_gates=8702
mirror_orbits=4351
orbit_size_hist={2:4351}
fixed_orbits=0
unpaired_or_duplicate_gates=0
low_rank_mirror_orbits=4333
high_rank_mirror_orbits=18
dead_rows_with_e_branch_low_rank_orbit=130/130
```

The quotient compresses HYP-3471's gate colors:

```text
typed_mirror_orbit_palette_size=186
structural_mirror_orbit_palette_size=82
delta_orbit_hist={2:3372,3:500,4:354,5:88,6:29,7:7,8:1}
```

The AP84 four-color packet becomes a two-orbit packet:

```text
('B0:5|E:0','E:0|B1:5')
('B1:7|E:0','E:0|B0:7')
```

It appears in `68/135` rows and `67/130` dead rows, exactly matching
HYP-3471's AP packet count after mirror quotienting.

## Hard Orbit Debt

The large cover-delta debt is finite and small:

```text
hard_orbits_delta_ge_7=8
hard_orbit_rows=[
  random_covering_022,
  random_covering_031,
  random_covering_049,
  random_covering_078,
  random_covering_080,
  random_covering_085,
  random_covering_113
]
```

`random_covering_022` has two hard orbits, one with the unique delta `8`.
The rest of the hard set has delta `7`.  HYP-3455's `random_covering_031`
mirror pair is therefore not isolated as the only hard orbit; it is one member
of a seven-row hard-orbit family.  That strengthens the proof target by naming
the complete current-bank family rather than a single example.

## Proof Pull

HYP-3475 upgrades HYP-3471:

```text
dead_components(row)>0
  => row has a low-rank E/branch mirror orbit.
```

The next finite lemma should discharge the eight hard same/cross-branch
two-sided mirror orbits.  Acceptable exits are:

```text
HYP-3455/HYP-3451 gluing or conductance route
endpoint-spine/wall lift
owner-current imbalance
two-adic descent
HYP-3460 phase-branch bypass
signed-SPEC/Rprime debt
```

The AP84 side is now a closed mirror packet through
HYP-3462/HYP-3470/HYP-3461/HYP-3460/HYP-3459/HYP-3458/HYP-3454/HYP-3456/HYP-3457.
HYP-3472 supplies the dead-cover current sidecar, HYP-3473 names the formal
producer obligations, and HYP-3474 guards which quotients may forget this
payload.  The non-AP transfer should focus on the seven hard random rows
above.

## Tournament Analysis

Vertices are mirror-orbit proof carriers, not runners, individual gates, or raw
counts.

```text
pairwise_observable =
  predicate retention + mirror payload + typed color payload
  + hard-delta localization + scalar firewall
score_hist={6:1,41:1,54:1,56:1,60:1,62:1,64:2}
directed_3cycles=0
hamiltonian_path =
  M00_exact_mirror_orbit_quotient
  -> M01_hard_delta_orbit_ledger
  -> M02_dead_row_e_branch_orbit
  -> M03_typed_mirror_color_pair
  -> M04_ap84_two_orbit_packet
  -> M05_phase_branch_hit_sidecar
  -> M06_single_gate_color_word
  -> M07_raw_gate_count
```

Assumption challenge: considered runners, residues, typed endpoint colors,
individual survivor gates, mirror gate pairs, fixed circle sections, component
boundaries, cover arcs, hard delta gates, and proof obligations.  The chosen
quotient preserves gate-extension gluing debt plus low-rank/E-branch escape
availability under `t -> 1-t`.  It destroys which member of a mirror pair
carries branch orientation, so the ordered typed/structural pair must stay as a
sidecar.
