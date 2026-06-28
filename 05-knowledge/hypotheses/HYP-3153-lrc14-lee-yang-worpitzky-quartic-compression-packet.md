---
id: HYP-3153
title: LRC14 Lee-Yang/Worpitzky/quartic compression packets may make the k=8 dip a bounded-degree certificate
status: RESERVED / executable packet scout pending; not a proof
source: codex-2026-06-28
tangent: T1218
technique: LTI-279
tournament_technique: LTT-177
script: 04-computation/lrc14_lee_yang_worpitzky_quartic_packet_codex_20260628.py
result: 05-knowledge/results/lrc14_lee_yang_worpitzky_quartic_packet_codex_20260628.out
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

## Reserved Claim

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

## Planned Measurements

1. Reuse the known k=8..13 cap/dip rows and record the transition
   `|P|=4 -> 5` as the first large off-circle correction.
2. Treat the cap as pair-normalized Pascal / de Moivre-Laplace mass and the
   dip as the fourth-order off-circle correction.
3. Audit the reciprocal-root relation `q0 = q6 * prod |root_i|` and the
   circle approximation `q0 ~= q6 * R^6` on named rows when available.
4. Package the k=8 dual line
   `L_y = p0 + p6 + (1/10)p3`
   as a bimodality functional, and compare it against the cap row.
5. Recheck the even/odd split: even side through the biquadratic
   `u^4 - 5u^2 + 4`, odd side through the K3/Worpitzky `-1/3` kernel.
6. Add the graph-theoretic witness language:
   strongly connected ear decompositions, factor-critical odd ears, and
   odd-cycle/Omega sidecars.
7. Run Tournament Analysis on proof carriers rather than runners, arcs, or
   raw class names.

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
