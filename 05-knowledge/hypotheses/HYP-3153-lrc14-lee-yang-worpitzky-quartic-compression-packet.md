---
id: HYP-3153
title: LRC14 Lee-Yang/Worpitzky/quartic compression packets may make the k=8 dip a bounded-degree certificate
status: SYNTHESIS / exact finite packet scout; not a proof
source: codex-2026-06-28
tangent: T1218
technique: LTI-279
tournament_technique: LTT-177
script: 04-computation/lrc14_lee_yang_worpitzky_quartic_packet_codex_20260628.py
result: 05-knowledge/results/lrc14_lee_yang_worpitzky_quartic_packet_codex_20260628.out
reflection: 07-reflections/lrc14-lee-yang-worpitzky-quartic-packet-codex-20260628.md
related:
  - HYP-3152
  - HYP-3151
  - HYP-3150
  - HYP-3149
  - HYP-3147
  - HYP-3142
  - HYP-3139
  - HYP-3136
  - HYP-3109
  - HYP-3108
  - HYP-3099
  - THM-577
  - OPEN-Q-108
---

# HYP-3153: Lee-Yang/Worpitzky/Quartic Compression Packet

## Claim

This lane tests whether HYP-3151's executable function-compression calculus,
HYP-3152's Lee-Yang circle/radius web, and the mac-mini HYP-3150 parity split
can be fused into a proof-facing packet:

```text
Lee-Yang circle zeros / Pascal mass / cap
+ off-circle phi4 correction / dip
+ Worpitzky odd edge kernel
+ quartic Newton-Maclaurin extremality at the AP
+ odd-ear witness sidecar
= bounded-degree k=8 certificate candidate.
```

The executable target is not a proof of LRC14.  It should produce a finite
certificate saying exactly which parts are algebraic identities, which parts
are measured signals, and which parts remain proof obligations.

## Evidence

The scout
`04-computation/lrc14_lee_yang_worpitzky_quartic_packet_codex_20260628.py`
verifies these packet components exactly, except where explicitly marked as a
numeric root-radius sidecar:

```text
q0 = q6 * R^6 exact by Vieta for consec k=8..13
root-radius spread numeric ratios = 1.1427, 1.2226, 1.2602, 1.3267, 1.3322, 1.3629
pair_mass - cap dips = 1081/76440, 1/4004, 0, 0, 0, 0
```

For the k=8/9/10 bimodality functional,

```text
L_y = q0 + q6 + q3/10 <= cap
k=8 margin = 683/29400
k=9 margin = 106901/2102100
k=10 margin = 69/910
```

For k=8 the scaled identity is exact:

```text
10q0 + q3 + 10q6
= 10S0 - 10S1 + 10S2 - 9S3 + 6S4
= 2633/735.
```

The correction split is:

```text
base = 1744/147
odd  = -9S3 = -2973/245
even =  6S4 = 944/245
|odd/even| = 3.149364
```

This keeps the HYP-3150/S71 message precise: the even side is the solvable
biquadratic face, while the odd side is the larger Worpitzky/ear sidecar.

The Worpitzky packet reuses HYP-3151's exact data:

```text
K3 C/T edge-flip kernel = [[0,1],[1/3,2/3]]
stationary(C,T) = (1/4,3/4)
nontrivial eigenvalue = -1/3
Eulerian rows = (1,4,1) and (1,11,11,1)
```

The even algebraic packet is:

```text
g(t) = (t-1)(t-2)(t-4)(t-5)
u=t-3 gives u^4 - 5u^2 + 4
v=u^2 gives v^2 - 5v + 4, discriminant 9.
```

The Newton/Maclaurin packet records all k=8 normalized moment defects as
negative:

```text
-103487/2160900,
-3670987/972405000,
-1099703/1944810000,
-226031/1944810000,
-109/2881200.
```

The ear sidecar is not promoted to a theorem proof, but the scout verifies the
finite witness grammar:

```text
C3 is factor-critical
C3 plus a length-3 odd ear remains factor-critical
directed C3 is strongly connected
directed C3 plus a directed ear remains strongly connected
```

## Tournament Fingerprint

Tournament Analysis uses proof packets as vertices.  The selected path is:

```text
HYP3152_leeyang_radius_root_curve
-> HYP3151_function_legality_packet
-> k8_Ly_bimodality_certificate
-> odd_worpitzky_ear_packet
-> even_biquadratic_galois_fold
-> newton_ap_violation_packet
-> pascal_cap_dip_table
-> generic_quintic_wall_alarm
-> raw_scalar_p0_only
```

Fingerprints:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
SCCs=singletons
hamiltonian_path_count=1
```

## Assumption Challenge

The tournament vertices for this lane should be proof obligations:

```text
Lee-Yang circle packet, Pascal cap packet, phi4 dip packet,
Worpitzky odd kernel, biquadratic even fold, Newton quartic AP inequality,
odd-ear witness packet, radius/product sidecar, raw scalar warning.
```

This deliberately challenges the assumption that a tournament must be built
from runners, raw arcs, or score classes.  The quotient is only legal when it
preserves the LRC predicate or emits the destroyed coordinate as a sidecar:
root radius, root angle spread, odd orientation, ear endpoint, fiber PGF, or
finite-address debt.

## Guardrail

The Abel-Ruffini language remains only a route alarm.  The exact statement to
test is smaller: the visible k=8 hard row decomposes into audited packets whose
effective algebraic degrees are at most four, with the even part already
folding to degree two.

## Remaining Proof Obligations

The scout leaves three load-bearing tasks:

1. Bound the off-circle dip/lambda, not merely describe it.
2. Prove the AP/consec `L_y` extremality in the required bank or route each
   exception to an existing finite-address/observer-gluing packet.
3. Carry the ordered/Worpitzky/ear sidecars through the LRC predicate before
   scalarizing to `p0`, `cap`, or a root radius.
