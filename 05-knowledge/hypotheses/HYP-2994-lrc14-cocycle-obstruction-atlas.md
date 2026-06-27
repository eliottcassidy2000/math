---
id: HYP-2994
title: LRC14 cocycle obstruction atlas
status: PROOF-INTERFACE / cochain-cocycle atlas and quotient guardrail; not a proof
source: codex-2026-06-24-S166
artifacts:
  - 04-computation/lrc14_cocycle_obstruction_atlas_codex_s166.py
  - 05-knowledge/results/lrc14_cocycle_obstruction_atlas_codex_s166.out
  - 07-reflections/lrc14-cocycle-obstruction-atlas-codex-s166.md
related:
  - HYP-2993
  - HYP-2992
  - HYP-2991
  - HYP-2990
  - HYP-2989
  - HYP-2988
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2984
  - HYP-2981
  - HYP-2979
  - HYP-2978
  - HYP-2969
  - HYP-2963
  - HYP-2956
  - HYP-2953
  - HYP-2887
  - HYP-2595
  - HYP-2594
  - THM-572
  - OPEN-Q-108
---

# HYP-2994: LRC14 Cocycle Obstruction Atlas

HYP-2994 is the cochain/cocycle lift of the recent Haar and zipper work.
HYP-2991 identifies the local fixed-margin cocycle

```text
zeta(T)=T00-T01-T10+T11.
```

HYP-2992 supplies the dyadic Haar tile interaction atlas, and HYP-2993 records
the zipper theorem pattern.  This pass asks what kind of obstruction each
carrier is: exact coboundary, closed cocycle, torsion/period class, Cech gluing
class, gauge with boundary stop, or named harmonic residual.

## Computation

Script:

```text
04-computation/lrc14_cocycle_obstruction_atlas_codex_s166.py
```

Stored output:

```text
05-knowledge/results/lrc14_cocycle_obstruction_atlas_codex_s166.out
```

The script records a proof-interface cochain ledger:

```text
C0: packet labels, owner potentials, safe sections, Fejer centers,
    exact-period residues, state-lift route names.
C1: certificate handoff arrows, endpoint transfers, smoothing gauges,
    Ramanujan phase shifts, source-spectrum pullbacks, section restrictions.
C2: fixed-margin Haar zeta squares, Haar tile curls, tope/cocircuit curls,
    color-resonance squares, boundary-moment curls, state-lift faces.
H2: unpaired mixed Haar mode, no-hidden-kernel survivor, Johnson-harmonic F7,
    THM-572 state-lift obstruction.
```

It then scores `15` carriers:

```text
fixed_margin_haar_zeta
haar_tile_interaction_curl
certificate_handoff_delta
tope_cocircuit_boundary_class
exposure_kernel_cech_class
fejer_interval_dual_coboundary
ramanujan_period_torsor
smoothing_homotopy_gauge
boundary_moment_curl
color_resonance_square
source_spectrum_pullback
apex_sheaf_gluing_class
octahedral_hodge_current
fixed_margin_johnson_F7
raw_scalar_shadow
```

The classification counts are:

```text
kinds={'1_coboundary_atlas': 1,
       'boundary_coboundary_or_stop': 1,
       'cech_1_cocycle': 2,
       'closed_2_cocycle': 1,
       'dual_coboundary': 1,
       'gauge_1_cochain_with_boundary': 1,
       'global_2_cocycle': 1,
       'harmonic_2_residual': 1,
       'moment_2_cocycle': 1,
       'named_harmonic_residual': 1,
       'negative_control': 1,
       'observer_1_cocycle': 1,
       'torsion_1_cocycle': 1,
       'typed_2_cocycle': 1}
cochain_levels={'C0-shadow': 1, 'C1': 7, 'C1/C2': 1, 'C2': 4, 'H2': 2}
```

The sparsely preserved labels are informative:

```text
exact_scale: 4 carriers
mixed_haar_sign: 3 carriers
phase_period: 4 carriers
```

Those are exactly the labels most likely to be accidentally destroyed by a
scalar quotient.

## Tournament Analysis

Vertices are cocycle carriers/proof obligations, not runners.  The pairwise
observable is retained LRC predicate payload plus obstruction-control vector:
predicate, fiber, exactness, closedness, torsion/period data, residual name,
formal checkability, anti-scalar guardrail, and computable next test.

Stored fingerprint:

```text
vertices=15
score_hist={0:1,1:1,2:1,3:1,4:1,6:3,8:1,9:1,10:1,11:1,12:1,13:1,14:1}
directed_3cycles=1
SCC_sizes=[3,1,1,1,1,1,1,1,1,1,1,1,1]
Hamiltonian_path_count=3
```

The score leader path is:

```text
certificate_handoff_delta
> fixed_margin_haar_zeta
> exposure_kernel_cech_class
> tope_cocircuit_boundary_class
> haar_tile_interaction_curl
> fejer_interval_dual_coboundary
> ramanujan_period_torsor
> smoothing_homotopy_gauge
> color_resonance_square
> fixed_margin_johnson_F7
> boundary_moment_curl
> source_spectrum_pullback
> octahedral_hodge_current
> apex_sheaf_gluing_class
> raw_scalar_shadow
```

The nontrivial `3`-carrier SCC is the important signal.  Certificate handoff,
local `zeta`, and exposure gluing are mutually entangled: a local mixed cocycle
is theorem-safe only when the handoff atlas records where it goes and the
exposure cover confirms no hidden kernel survives.

## Obstruction Face Dictionary

HYP-2994 names six reusable faces:

```text
fixed_margin_square      closed local 2-cocycle
haar_fejer_dual_square   primal/dual coboundary check
endpoint_tope_boundary   boundary coboundary or declared stop
period_smoothing_square  torsion gauge
source_exposure_cech     Cech 1-cocycle
octahedral_harmonic      harmonic residual
```

This is the dictionary future agents should use when deciding whether a new
quotient is exact, closed, torsion, or residual.

## Readout

1. Exact coboundaries are the only safe way to forget packet potentials.
2. Closed but non-exact local squares are useful only when a dual certificate
   pairs them or a named residual receives them.
3. Torsion/period classes are real proof data; squarefree and gcd quotients
   erase live prime-power packets.
4. `F7` should be defined as a preserved harmonic residual sector, not as the
   complement of known tests.
5. Raw scalar shadows preserve the slogan but destroy nearly every obstruction
   coordinate.

## Assumption Challenge

Candidate tournament vertices considered: runners, endpoint owners, fixed
circle sections, section boundaries, wall-crossing events, residues, cover
arcs, Fourier/Haar modes, Ramanujan packets, Fejer atom banks, matroid/tope
circuits, source-spectrum pullbacks, Cech overlaps, harmonic residuals, and
proof obligations.

Chosen vertices: proof carriers / cocycle obligations.

Preserved LRC predicate:

```text
strict witness vs AP/GW equality vs named state-lift residual
```

Destroyed information if scalarized:

```text
exact scale, mixed Haar sign, phase period, endpoint-owner topology,
packet family, dual-certificate identity, and residual route.
```

The challenged assumption is that one scalar invariant or one runner-based
tournament can close the proof.  The cocycle atlas says the proof object is a
labelled cochain complex, and every allowed quotient must declare whether it
is taking a coboundary, a closed cocycle, a torsion class, or a residual.

## Next Work

1. Attach `zeta` signatures to HYP-2963 packet rows and count independent
   color-compatible mixed modes.
2. Tensor Ramanujan `c_q` labels with Haar tile classes to catch `q=25,27,36,
   63,84,98,168,280,4312` side channels.
3. Group Fejer interval certificates by identical cocycle signature before
   interval generation.
4. Build a no-hidden-kernel Cech check: every exposure cover overlap glues or
   emits a named residual.
5. Define an executable `F7` record with fields:

```text
packet_fiber
mixed_cocycle
harmonic_sector
state_lift_route
failed_exits
```
