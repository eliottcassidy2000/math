---
id: HYP-3135
title: LRC14 resolvent-packet synthesis: the 120/320 De Moivre clue is a middle elementary-symmetric payload, and the remaining proof should be packaged as a bounded-core plus signed-SPEC packet theorem
status: SYNTHESIS / proof-target refinement; not a proof
source: codex-2026-06-27-S272
script: 04-computation/lrc14_resolvent_packet_synthesis_codex_s272.py
result: 05-knowledge/results/lrc14_resolvent_packet_synthesis_codex_s272.out
extends:
  - HYP-3136
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3131
  - HYP-3130
  - HYP-3129
  - HYP-3128
  - HYP-3127
related:
  - HYP-3125
  - HYP-3124
  - HYP-3122
  - HYP-3116
  - HYP-3111
  - HYP-3110
  - OPEN-Q-108
external: De Moivre quintic, resolvent quartic, Lee-Yang, Asano contraction, phi4, elementary symmetric functions
---

# HYP-3135: LRC14 Resolvent-Packet Middle-Layer Synthesis

## Claim

The user's quintic example should be read as a packet lesson, not as a
standalone numerology.  For the resolvent roots

```text
2, -4, 8, -16
```

the quartic coefficients

```text
x^4 + 10*x^3 - 120*x^2 - 320*x + 1024
```

come from elementary symmetric layers:

```text
e1 = -10, e2 = -120, e3 = 320, e4 = 1024.
```

Thus `120` and `320` are the pair and triple branch-interaction payloads of
the resolvent orbit.  The real root is the real fifth-root branch sum

```text
fifthroot(2) - fifthroot(4) + fifthroot(8) - fifthroot(16),
```

which satisfies `x^5 + 20*x^3 + 20*x^2 + 30*x + 10 = 0`.

The LRC14 analogue is direct: the remaining proof should not chase a single
raw scalar (`p0`, raw Lee-Yang radius, raw A000568 count, raw circuit size, or
raw De Moivre residual).  It should preserve the middle packet:

```text
(Q-floor constants,
 signed SPEC low part,
 SPEC Parseval tail,
 bounded-core Lee-Yang/phi4 status,
 far-push-out monotonicity,
 edge tail/tip deletion sectors,
 finite-address / observer-gluing exit).
```

## Incoming work integrated

HYP-3136 assembles the multi-far floor closure: the covering factorization
`L(S)=Rprime*meas(R-safe)*meas(Q-lonely)` now has the `Q` block, the `R-safe`
wide-V reduction, and the signed-SPEC `Rprime` certificate in place, with the
remaining work reduced to a finite constant chase.  HYP-3132 sharpens the
bounded-core residue: the whole remaining core is the k=8 dip, whose
gK8/Delsarte dual quartic `(t-1)(t-2)(t-4)(t-5)` folds under `u=t-3` to the
biquadratic `u^4 - 5u^2 + 4`, with discriminant `9`.  That is the local
phi4/De Moivre payload this HYP-3135 packet is trying to preserve.

HYP-3134 (incoming A000568 edge-envelope quotient) sharpens HYP-3133: raw
four-sector edge decks are a lower envelope, paired tail/tip child packets are
a safe upper envelope, and A000568 sits between them as a global-consistency
quotient.  That is exactly the quotient discipline this resolvent packet needs:
do not forget a local payload unless a named gluing rule preserves the LRC
predicate.

HYP-3133 (incoming A000568 edge-sandwich) adds the adjacent finite quotient
shadow: sector words, A000568 one-extension shadows, and paired tail/tip child
decks are a local stratifier for the HYP-3129/HYP-3132 constant chase.  This
resolvent packet is the algebraic version of the same warning: keep the middle
payload before quotienting.

HYP-3131 changes the shape of the multi-far problem.  Far elements are not the
obstruction; they push the miss-PGF zeros outward from the good bounded base.
The multi-far Lee-Yang region therefore reduces to bounded-core Lee-Yang plus
far-pushes-out monotonicity.

HYP-3130 closes the `Q`/apex block by a true smooth minorant and exact
cap-matching constants `c_2..c_6`.  It also proves the negative result that
absolute envelopes cannot close the coupling: `B1/h0 > 1`, so sign
cancellation is essential.

HYP-3129 supplies the signed mechanism.  The `Rprime` floor is an elementary
resonance-lattice/SPEC calculation, not an Elliott-Halberstam or
Bombieri-Vinogradov problem.  The existing certificate gives
`Rprime >= 0.64178` on the tested multi-far family, with the honest remaining
work a closed-form constant chase over all `(R,Q)` with `|Q| <= 6`.

HYP-3128 explains why naive Asano cannot certify the joint floor alone: the
`Q` block is Lee-Yang, but the `R` block with at least seven speeds is not
unit-disk zero-free.  Asano is still useful as a packet/far-push language, but
the coupling must pass through signed SPEC or an equivalent middle-layer
certificate.

## Where 120 and 320 have appeared before

The repo has seen `120` and `320` mostly as finite middle-layer signals:

- `05-knowledge/variables/cycle-counts.md` has the n=9 coefficient formula
  containing `-120*bc33`, a cycle-pair correction term.
- HYP-2852 in `00-navigation/SESSION-LOG.md` recorded a provable `0/320`
  bounded resonant-center floor check.
- HYP-2608 recorded `320` random six-support vectors for the residue-addressed
  reciprocal-sum correction, where signed cancellation was the real content.
- The cap/extremality ledgers record a `+0.320` slack signal in high rows:
  useful margin, not a theorem by itself.
- `120` recurs as `5!` labelled tournament mass, fixed-path profile scale,
  n=15 torsion moat/missed-cell count, and q^120 modular support horizon.

The shared reading is not that `120` or `320` are magic.  They usually mark a
finite middle layer: pair/triple interactions, finite stress horizons, or a
bounded check where cancellation has not been safely quotiented away.

## Improved proof target

The current endgame should be stated as a packet theorem.

1. Non-covering or direct-witness branches stay routed through the existing
   packet sheaf: q-witness, AP/GW boundary atom, C27/K33/state-lift, Fejer or
   Ramanujan certificates, or named debt.
2. In the covering branch write `S = R union 14Q`, `2 <= |Q| <= 6`.
3. Use HYP-3130/HYP-3128 to close `meas(Q-lonely) >= c_r > 0`.
4. Use HYP-3131 to remove far placements as independent obstructions: far
   additions increase coverage and push zeros outward, so the binding case is
   the bounded core.
5. Use HYP-3129 to make the coupling explicit:

```text
Rprime = 1 + SPEC / (meas(R-safe)*meas(Q-lonely)),
SPEC = exact low resonance on 14Z + Parseval-controlled signed tail.
```

6. Finish by proving the closed-form version of the HYP-3129/HYP-3132 constant
   chase and the S70 k=8 biquadratic coefficient bound behind the bounded-core
   Lee-Yang/phi4-to-`Rprime` bridge.

This is the resolvent lesson in LRC language: the proof is allowed to forget
raw branches only after the pair/triple payload has been recorded.  For LRC14
that payload is the signed SPEC/Lee-Yang/edge-witness packet.

## Challenged assumption

Do not assume the proof vertices are runners, arcs, primes, roots, or
tournament classes.  For this route the proof vertices are packet coordinates:
`Q-floor`, `signed_SPEC_low`, `SPEC_tail`, `bounded_core_rho`,
`far_push_status`, `edge_tail_child`, `edge_tip_child`, and
`terminal_exit_or_named_debt`.  Raw scalar quotients destroy the same middle
information that the resolvent quartic keeps in `-120` and `320`.

## Next concrete work

1. Turn the HYP-3129 per-row certificate into a closed-form uniform lower
   bound by bounding `SPEC_low` and the Parseval tail at a fixed `M`.
2. Prove far-pushes-out for all legal far placements, not only the S69 scan.
3. Formalize the bounded-core implication `rho > 1 => Rprime >= c`, using the
   phi4/kappa4 row as the finite hard case.
4. Make the S70 reflection fold rigorous: k=8 dip bound equals a degree-2
   bound on the biquadratic `u^2` coefficient after the `s -> 6-s` symmetry.
5. Attach the HYP-3124/HYP-3125 edge-witness fields to the covering packet so
   any collapsed sector remembers tail/tip deletion children and the
   cross-sector orientation word.
