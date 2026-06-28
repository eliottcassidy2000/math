---
id: HYP-3421
title: LRC14 off-grid resonance transparency and Rprime closure
status: SYNTHESIS / exact-rational transparency scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1382
technique: LTI-382
tournament_technique: LTT-282
script: 04-computation/lrc14_offgrid_resonance_transparency_rprime_closure_codex_20260628.py
result: 05-knowledge/results/lrc14_offgrid_resonance_transparency_rprime_closure_codex_20260628.out
reflection: 07-reflections/lrc14-offgrid-resonance-transparency-rprime-closure-codex-20260628.md
related:
  - HYP-3418
  - HYP-3417
  - HYP-3416
  - HYP-3415
  - HYP-3414
  - HYP-3412
  - HYP-3410
  - HYP-3310
  - HYP-3266
  - HYP-3265
  - HYP-3255
  - HYP-3140
  - HYP-3136
  - HYP-3129
  - HYP-3125
  - HYP-3124
  - HYP-2896
  - THM-523
  - OPEN-Q-108
---

# HYP-3421: LRC14 Off-Grid Resonance Transparency And Rprime Closure

## Claim

The "resonant survivor" obstruction should be recast as a finite
full-optimum off-grid transparency ledger feeding the covering floor.  This
is the exact-rational companion to HYP-3415, now corrected by HYP-3418:
HYP-3415 reduces the LRC14 proof to the covering-case decorrelation floor
`R' > 0`; HYP-3421 attacks the resonant part of that floor by replacing it
with transparent off-grid packet checks plus a symbolic `Rprime` constant
chase; HYP-3418 warns that this does **not** justify a coprime-to-14 or
non-resonant-only reduction, because the odd subpacket prefers `t=1/2` where
the even speeds die.  The remaining floor is 2-adic/even-speed descent plus
signed-SPEC/fiber-PGF decorrelation.

The on-grid core is the six-unit equioscillation set

```text
(1,13), +(3,11), +(5,9) at a/14, a in (Z/14)^*.
```

Multiples of `14` attack exactly this core: at every unit grid point `a/14`
each `14Q` runner sits at the observer.  The useful witnesses then move
off-grid.  On the checked packets, every resonant speed, meaning every speed
divisible by `2` or `7` and in particular every `14Q` tip, is already safe at
the selected off-grid optimum.

So the proof route becomes:

```text
finite off-grid transparency / positive packet classifier
  + 2-adic even-speed descent
  + Rprime signed-SPEC constant chase
  + AP/GW finite equioscillation rigidity
  + owner/state-lift terminal router
  => LRC14
```

This is not a completed proof.  It sharpens the open O12 obligation from
HYP-3266, the HYP-3415 critical-path map, and the HYP-3418 2-adic correction:
prove the transparency classifier for every residual packet at the full
optimum, then close the even-speed descent and HYP-3129's all-packet
`Rprime >= c` constant chase.

## Exact Scout

The script
`04-computation/lrc14_offgrid_resonance_transparency_rprime_closure_codex_20260628.py`
performs an exact rational pair-crossing audit.  Stored output is in
`05-knowledge/results/lrc14_offgrid_resonance_transparency_rprime_closure_codex_20260628.out`.

Named rows:

| case | `M` | selected `t*` | off `14`-grid | resonant distance range | `14Q` min |
|---|---:|---:|---:|---:|---:|
| q12 seed core `{1..11,13}` | `1/12` | `1/12` | yes | `1/6..1/2` | n/a |
| even covering row `{2..14}` | `1/8` | `1/16` | yes | `1/8..1/2` | `1/8` |
| many-`14Q` row | `1/9` | `1/9` | yes | `1/9..4/9` | `1/9` |
| many-`14Q` row | `1/8` | `15/112` | yes | `1/8..1/2` | `1/8` |
| canonical `{1..11,13,84}` | `7/89` | `37/89` | yes | `7/89..44/89` | `7/89` |

The canonical one-tail tower verifies exactly through `m=8`:

```text
S_m = {1,2,3,4,5,6,7,8,9,10,11,13,84m}
t_m = (35m+2)/(84m+5)
M(S_m) = 7m/(84m+5)
active pair = (5,84m)
```

Thus the most resonant one-tail family is not an obstruction.  The `14Q`
runner is one of the binding flanks, and the margin is

```text
7m/(84m+5) - 1/14 > 0.
```

## Proof Obligations

The scout names six proof-facing obligations.

1. `O12a finite_offgrid_transparency`: every O12 residual must reach one of
   the denominator floors, the `84m` formula, a positive owner packet, or a
   named residual row.
2. `O12b divisor_grid_localization`: `14Q` tips kill the `14`-grid but are
   safe on certified off-grid cells.
3. `O12c signed_SPEC_constant_chase`: turn HYP-3129's per-row exact low part
   plus Parseval tail into a closed-form all-packet `Rprime >= c` theorem.
4. `O12d fiber_PGF_translation`: keep HYP-3140's identity
   `Rprime = E[N_R | Q]/E[N_R]` before scalarizing.
5. `O12e two_adic_even_descent`: incorporate HYP-3418's correction that
   coprime-to-14 reduction fails and the covering floor is driven by even
   speeds under the halving map `u=2t`.
6. `O15 finite_equioscillation_core`: finish AP/GW rigidity from the three
   antipodal unit binding pairs plus blind sidecars.
7. `GLUE terminal_factorization`: combine transparency or
   `L=Rprime*R_safe*Q_lonely` with owner/state-lift terminal routers.

## Geometry Reframe

The right picture is not "resonance survives."  It is:

```text
14-grid core      = finite equioscillation / AP-GW boundary
off-grid bulk     = full-optimum transparency / Rprime decorrelation
even speeds       = 2-adic binding floor under u=2t
14Q resonance     = grid-local danger, off-grid guard flank
```

This matches the `14=2*7` split from HYP-3310 but also incorporates
HYP-3418's warning about which prime is proof-critical.  The unit `C3`
skeleton and apex-7 field organize the off-path equality/census regime, while
the covering floor is driven by the `2`-adic even-speed layer.  The
`12 -> 24` hinge remains a 2-adic magnitude sidecar; it must not be forgotten
as a residue-only move.

## Tournament Analysis

Tournament vertices are proof carriers, not runners, arcs, raw residues, or
floor denominators.  Pairwise observable: majority over retained LRC
predicate, finite exactness, resonance transparency, `Rprime` interface,
core-rigidity interface, formalization readiness, quotient legality, and
failure guardrails.

The scout reports:

```text
vertices = 8
edges = 28
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles = 0
hamiltonian_path_count = 1
selected path =
  finite_offgrid_transparency_lemma
  -> canonical_84m_binding_formula
  -> two_adic_even_speed_descent
  -> signed_SPEC_Rprime_constant_chase
  -> edge_witness_tail_tip_packet
  -> fiber_PGF_conditional_moment
  -> owner_cut_terminal_router
  -> unit_core_equioscillation_rigidity
  -> raw_resonant_survivor_worry
```

The leading vertex says what to prove next: an all-packet transparency
classifier.  The third vertex says what remains after that: the symbolic
closed-form `Rprime` constant chase.

## Status

This is evidence and proof-route synthesis.  It does not prove the all-packet
transparency classifier, and it does not replace HYP-3129's remaining
closed-form SPEC constant chase.  It does identify a stronger completion route:
retire raw resonant-survivor positivity as a standalone analytic worry and
replace it with exact off-grid packet transparency plus the signed-SPEC/fiber
PGF/edge-witness `Rprime` floor and HYP-3418's 2-adic even-speed descent.
