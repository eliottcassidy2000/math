---
id: HYP-3134
title: LRC14 integrated multi-far floor closure
status: SYNTHESIS / reduction to finite constant chase; renumbered from an index-only HYP-3132 row after S70 claimed file-backed HYP-3132
source: codex-2026-06-27 rebase integration
related:
  - HYP-3133
  - HYP-3132
  - HYP-3131
  - HYP-3130
  - HYP-3129
  - HYP-3128
  - HYP-3127
  - HYP-3126
  - HYP-3125
  - HYP-3124
  - THM-546
  - THM-530
  - OPEN-Q-108
---

# HYP-3134: Integrated Multi-Far Floor Closure

This card preserves the earlier index-only "integrated multi-far floor"
claim after incoming S70 claimed the file-backed HYP-3132 namespace for the
k=8 De Moivre biquadratic resolvent.

For a covering packet `S = R union 14Q`, with `R` 14-free and
`2 <= |Q| <= 6`, the proof target factors as

```text
L(S) = Rprime(S) * meas(R-safe) * meas(Q-lonely).
```

The current synthesis assigns the three factors to three proof tools:

- `meas(Q-lonely) >= c_r > 0`: HYP-3130 supplies the Gaussian/Beurling-
  Selberg minorant and cap-matching Q/apex floor, modulo the THM-546 peel.
- `meas(R-safe) > 0`: HYP-3125/HYP-3126 supply wide-V Gaussian/equidistribution
  decoupling to a bounded core.
- `Rprime(S) >= c`: HYP-3129 supplies the retained-edge signed SPEC
  certificate, with representative-family floor `>= 0.64178` and remaining
  closed-form finite constant chase.

HYP-3131 says far tips push Lee-Yang zeros outward, so the far side is
stabilizing rather than binding.  HYP-3132 sharpens the bounded-core tail:
the remaining hard node is the k=8 phi4/De Moivre biquadratic resolvent.
HYP-3133 adds an A000568 middle shadow for stratifying the finite SPEC rows
before escalating to paired tail/tip children or named resonance debt.

Edge-witness reading: the tip child carries Q/apex minorant, zero-free, and
far-push certificates; the tail child carries the bounded core, SPEC lattice,
and k=8 resolvent obligation; the cross-sector orientation word carries the
signed coupling that becomes `Rprime`.
