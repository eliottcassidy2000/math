---
id: HYP-3459
title: LRC14 coloring/discrepancy bridge
status: SYNTHESIS / exact AP84 color-quotient audit; not an LRC14 proof
source: codex-2026-06-29 continuation of prior coloring/discrepancy work plus HYP-3438, HYP-3441, HYP-3456, and HYP-3457
tangent: T1419
technique: LTI-419
tournament_technique: LTT-319
script: 04-computation/lrc14_coloring_discrepancy_bridge_codex_20260629.py
result: 05-knowledge/results/lrc14_coloring_discrepancy_bridge_codex_20260629.out
reflection: 07-reflections/lrc14-coloring-discrepancy-bridge-codex-20260629.md
related:
  - HYP-3458
  - HYP-3457
  - HYP-3456
  - HYP-3455
  - HYP-3454
  - HYP-3453
  - HYP-3452
  - HYP-3441
  - HYP-3438
  - HYP-3437
  - HYP-3436
  - HYP-3431
  - HYP-2991
  - HYP-2990
  - HYP-2263
  - HYP-2247
  - HYP-1802
  - THM-523
  - OPEN-Q-108
---

# HYP-3459: LRC14 Coloring/Discrepancy Bridge

## Claim

The useful common lesson from the repo's prior coloring work is:

```text
a coloring is a quotient, and a quotient is legal only when it preserves the
predicate being proved or emits the lost coordinate as named debt.
```

HYP-3459 applies that rule to the current AP84 tail package for

```text
S_m={1,2,...,11,13,84m}.
```

It connects:

- centered tournament edge variables `s_e=A_e-1/2` and the W-polynomial
  product rule;
- distance-graph LRC as a regular circular coloring problem;
- the old colored discrepancy/CRT finite-placement layer;
- HYP-2991's Haar zipper cocycle, where margins need the local `zeta`
  side-channel;
- HYP-2247's coloring-extension rank warning;
- HYP-3438's canonical mod-`35` survivor-gate law;
- HYP-3441's incident `C3/Qsqrt(-7)` residue router;
- HYP-3456/HYP-3457's AP84 floor-count and finite-transient packets.

The exact result is a small but important guardrail: the AP84 gate palette,
floor-correction word, and endpoint phase are distinct colors.  A proof may
use one only if it also carries or legally routes the others.

## Exact Readout

Script:

```text
04-computation/lrc14_coloring_discrepancy_bridge_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_coloring_discrepancy_bridge_codex_20260629.out
```

For residues `r=1..35`, HYP-3438's canonical gate law gives:

```text
outer gate present iff 7 does not divide r
inner gate present iff 5 does not divide r
```

so the gate-bucket histogram is

```text
{
  'both_outer_inner': 24,
  'outer_only_7_gate': 6,
  'inner_only_5_gate': 4,
  'clean_only_lcm35': 1
}
```

HYP-3456 gives the floor correction

```text
N(m)=floor(12m/35)+d[(m-1) mod 35],
```

with correction histogram

```text
{0:1, 1:22, 2:12}.
```

The collision audit is:

```text
gate_bucket -> correction values:
  both_outer_inner: [1,2]
  clean_only_lcm35: [0]
  inner_only_5_gate: [1]
  outer_only_7_gate: [1]

correction -> gate buckets:
  d=0: ['clean_only_lcm35']
  d=1: ['both_outer_inner','inner_only_5_gate','outer_only_7_gate']
  d=2: ['both_outer_inner']

mixed_gate_fibers_count=1
mixed_correction_fibers_count=1
mixed_gate_plus_correction_fibers_count=0
```

Thus neither the gate color nor the correction color alone is a legal
replacement for the other, while the pair is clean for the AP84 local target.

There is a further nonperiodic endpoint-phase sidecar:

```text
m=1 and m=36 have the same residue r=1,
but m=1 is a finite mixed endpoint transient,
while m=36 is in HYP-3454's rank-one E/E tail phase.
```

So the legal AP84 packet must retain the residue gate, floor word, linear
height level, endpoint phase, branch mirror, incident-core color, and Haar
zipper side-channel.

## Discrepancy Readout

The period-`35` word is not flat.  Exact cyclic-interval and dyadic-Haar
contrasts include:

```text
outer_gate       max cyclic discrepancy 6/7
inner_gate       max cyclic discrepancy 4/5
both_gate        max cyclic discrepancy 64/35
clean_gate       max cyclic discrepancy 34/35
correction       max cyclic discrepancy 93/35
correction_eq_2  max cyclic discrepancy 88/35
```

These are not bad news; they identify where a scalar equidistribution proof
would need the HYP-2991 zipper coordinate or a named sidecar.  The bridge is
therefore a colored discrepancy theorem, not a raw discrepancy theorem.

## Candidate Labelled Packet Theorem

For the AP84 tail, a legal proof packet should carry at least:

```text
1. residue gate color: outer iff 7 not divide m, inner iff 5 not divide m;
2. floor word: N(m)=floor(12m/35)+d_r;
3. linear height color: floor(12m/35), not only r mod 35;
4. endpoint phase: m=1..4 mixed transient versus rank-one E/E tail;
5. branch mirror color: C1/C0 and branch0/branch1 ownership;
6. incident core color: C3 slot plus Qsqrt(-7) sign from HYP-3441;
7. Haar zipper color: every collapsed margin retains or kills its cocycle.
```

Candidate theorem:

```text
Any AP84 splice into HYP-3439 is legal if the map from local survivor gates to
global branch-union escapes is a homomorphism for this color-packet product,
or else the first failed color is routed to HYP-3438/HYP-3453/HYP-3455,
owner-current, two-adic descent, exact-period/state-lift, or SPEC debt.
```

This theorem would not finish LRC14 by itself.  It would make the AP84 splice
and survivor-gate transfer proof-safe, removing one of the current quotient
ambiguities.

## Tournament Analysis

Vertices are quotient/color proof carriers, not runners, arcs, or raw
residues.

```text
pairwise_observable =
  preserved AP-tail predicate + lost-coordinate sidecar + finite-checkability
switch_gauge =
  higher retained color-packet score; ties by labelled-packet priority
score_hist={20:1,36:1,52:1,58:1,62:1,65:1,67:1,68:1,80:1}
directed_3cycles=0
hamiltonian_path_count=1
hamiltonian_path =
  labelled_color_packet_theorem
  -> residue_gate_plus_floor_word
  -> haar_zipper_cocycle_repair
  -> endpoint_phase_sidecar
  -> incident_C3_Qsqrt_router
  -> branch_mask_discrepancy_word
  -> distance_graph_regular_coloring
  -> raw_mod35_gate_color
  -> raw_scalar_escape_count
```

Assumption challenge: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, Haar
rectangles, centered edge signs, endpoint owners, incident cores, survivor
gates, and proof obligations were considered.  The chosen quotient is the
color-packet proof obligation.  It preserves AP-tail gate availability, floor
count, endpoint phase, and sidecar debt declaration, but destroys raw runner
identity, non-AP geometry, and arbitrary primitive-row adjacency.  The
challenged assumption is that a coloring quotient is harmless when its classes
look balanced; HYP-3459 shows balance is insufficient unless floor, endpoint,
branch, and zipper colors remain reconstructible.
