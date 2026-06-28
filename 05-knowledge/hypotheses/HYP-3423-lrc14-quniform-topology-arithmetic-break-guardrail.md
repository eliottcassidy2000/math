---
id: HYP-3423
title: LRC14 q-uniform topology / q-specific arithmetic break guardrail
status: SYNTHESIS / executable proof-route guardrail; not an LRC14 proof
source: codex-2026-06-28 synthesis of HYP-3312, HYP-3411, HYP-3413, HYP-3415, HYP-3416, HYP-3417, S259/HYP-3418, and HYP-3419
tangent: T1384
technique: LTI-384
tournament_technique: LTT-284
script: 04-computation/lrc14_quniform_topology_arithmetic_break_codex_20260628.py
result: 05-knowledge/results/lrc14_quniform_topology_arithmetic_break_codex_20260628.out
reflection: 07-reflections/lrc14-quniform-topology-arithmetic-break-codex-20260628.md
related:
  - HYP-3422
  - HYP-3421
  - HYP-3420
  - HYP-3419
  - HYP-3418
  - HYP-3417
  - HYP-3416
  - HYP-3415
  - HYP-3414
  - HYP-3413
  - HYP-3411
  - HYP-3312
  - HYP-3311
  - HYP-3310
  - HYP-3406
  - HYP-3405
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3423: LRC14 Q-Uniform Topology / Q-Specific Arithmetic Break Guardrail

## Claim

The C2 = Borsuk-Ulam naming in HYP-3312 is not merely a nice dictionary
entry.  It is a proof guardrail:

```text
C2/Borsuk-Ulam topology is q-uniform.
Therefore it can close only q-uniform residue/equioscillation obligations.
Any q-specific magnitude conclusion needs an arithmetic or finite packet.
```

In the LRC14 language:

```text
13 runners = two regular C6 orbits + fixed apex 7
           = units U + even shadow 2U + {7}.
```

The two subfield packets carry the two uniform proof halves:

```text
Q(cos 2pi/7): real cubic, C3 trace, cap/equioscillation residue half
Qsqrt(-7): imaginary quadratic, C2/BU sign, Gauss/floor orientation
```

But the Goddyn-Wong magnitude break is not uniform in q.  HYP-3413 records the
current switch:

```text
canonical GW doubling (n-2)->2(n-2) is ON iff q == 1 mod 3.
```

Thus the requested contrast is structural:

```text
q=4,7: magnitude break ON
q=5,6: magnitude break OFF
```

The C2/BU charge is present in all four rows.  It cannot be the switch.

## Executable Readout

Script:

```text
04-computation/lrc14_quniform_topology_arithmetic_break_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_quniform_topology_arithmetic_break_codex_20260628.out
```

The script prints the mod-14 orbit decomposition:

```text
[1,3,5,9,11,13]    unit/binding orbit
[2,4,6,8,10,12]    even/covering orbit
[7]                fixed apex
```

It then scans `q=3..22`.  The C2/BU residue charge is present on all
`20/20` rows, while the GW magnitude switch is ON on only `7/20` rows:

```text
q == 1 mod 3: q=4,7,10,13,16,19,22
```

For prime `q`, this is the Eisenstein/C3 gate from HYP-3413.  For composite
rows it is retained as the same verified q-mod-3 switch, but not promoted to a
prime-field unit-group claim.

## The Labelled Packet Theorem Target

The theorem target is not "topology proves the magnitude hinge."  The theorem
target is a legality rule for proof quotients:

```text
If a quotient uses C2/BU or C6 topology, it may certify a residue or
equioscillation obligation.

If the same quotient claims a q-specific magnitude conclusion, it must retain
or restore at least one of:
  1. HYP-3413 q mod 3 / Eisenstein arithmetic switch;
  2. HYP-3417 labelled owner-current packet for local mixed fibers;
  3. S259/HYP-3418 two-adic covering-floor descent;
  4. HYP-3415 decorrelation floor / Rprime inequality.
```

This is the guardrail that HYP-3416's recursive quotient ladder needed.  A
forgetful map is legal only when the forgotten coordinate is not needed by
the next theorem obligation, or when a named sidecar resurrects it.

Rebase correction: S259 already claimed HYP-3418 for the sharper statement
that the covering floor is two-adic, with even speeds as the binding
obstruction; incoming S297/S298 claimed HYP-3420/T1381/LTI-381/LTT-281 for
owner-cut chiral synthesis; and later mainline claimed HYP-3421/T1382 and
HYP-3422/T1383 for off-grid resonance transparency and two-adic relocation.
This guardrail is therefore HYP-3423/T1384/LTI-384/LTT-284.  It treats
S259/HYP-3418, HYP-3421, and HYP-3422 as arithmetic/floor sidecars that a
magnitude proof must retain, and HYP-3420 as the adjacent owner-cut/chiral
sidecar lane.

## Integration With HYP-3417

Incoming HYP-3417 found a local frontier current:

```text
{2:g2, 11:g1, 13:g1}
```

This is exactly one even-cover label plus two binding labels.  HYP-3423 says
how to use that signal: it is a local owner-current certificate compatible
with the residue/magnitude split, not a global substitute for the q-mod-3
arithmetic switch and not a substitute for the HYP-3415 covering floor.

After S259/HYP-3418 and HYP-3419, that same frontier current has a more
specific reading: `2:g2` is the even-cover / two-adic label.  Any future
owner-cut theorem should track whether hard frontiers systematically require
an even-cover label before importing apex-7 or topological language.

In other words:

```text
C2/BU topology: residue charge
HYP-3413 q mod 3: GW magnitude census switch
HYP-3417 owner current: local labelled finite discharge
S259/HYP-3418: two-adic covering-floor descent signal
HYP-3415 Rprime floor: critical-path inequality
```

The proof may move between these only when the packet labels are retained.

## Tournament Analysis

Vertices are proof-route obligations, not runners:

```text
HYP3415_decorrelation_floor_Rprime
HYP3417_labelled_owner_current_packet
HYP3413_q_mod_3_GW_arithmetic_switch
HYP3416_recursive_quotient_guardrail
HYP3411_C6_two_orbits_fixed_apex_packet
HYP3312_C2_BU_topological_charge
real_cubic_C3_trace_equioscillation
raw_topology_closes_magnitude_false_route
```

Pairwise observable:

```text
preserved proof coordinate plus forbidden forgetting debt
```

Switch/gauge:

```text
higher route score; ties by declared priority
```

Fingerprint:

```text
vertices = 8
score_hist = {-27:1, 13:3, 14:1, 19:1, 30:1, 50:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
selected_path =
  HYP3415_decorrelation_floor_Rprime
  -> HYP3417_labelled_owner_current_packet
  -> HYP3413_q_mod_3_GW_arithmetic_switch
  -> HYP3416_recursive_quotient_guardrail
  -> HYP3411_C6_two_orbits_fixed_apex_packet
  -> HYP3312_C2_BU_topological_charge
  -> real_cubic_C3_trace_equioscillation
  -> raw_topology_closes_magnitude_false_route
```

The ranking deliberately puts the critical covering floor first, then the
finite owner-current sidecar, then the q-specific arithmetic switch.  The
uniform topology route is useful but lower-ranked because it cannot select
the ON/OFF magnitude rows.

## Assumption Challenge

Alternative vertices considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing
events, residues, C6 orbits, subfields, q-rows, owner labels, floor packets,
matroid circuits, Fourier modes, proof obligations, and quotient policies.
```

Chosen vertices:

```text
proof-route obligations.
```

This quotient preserves the LRC predicate each route can certify.  It destroys
raw runner identity and raw q-row detail, but it records the destroyed
coordinate as forbidden forgetting debt.  The challenged assumption is that a
q-uniform topological invariant can prove a q-specific magnitude hinge.

## Status

This is not a proof of LRC14.  It is a rigor guardrail and a work allocator.
It says:

```text
Residue/equioscillation half: topology and Galois charges are legal.
Magnitude/census/floor half: arithmetic, owner-current packets, or Rprime floor
must enter.
```

The highest-leverage terminal path remains HYP-3415's decorrelation floor
unless the HYP-3417 owner-current packet can be shown to feed that inequality.
