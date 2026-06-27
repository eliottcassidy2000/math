---
id: HYP-3140
title: LRC14 fiber-PGF certificate: the multi-far Rprime floor is a conditional first-moment inequality over the 14 lifts of u=14t
status: EVIDENCE / exact generating-function scout; not a proof
source: codex-2026-06-27-S273
tangent: T1205
technique: LTI-266
tournament_technique: LTT-164
script: 04-computation/lrc14_fiber_pgf_certificate_codex_s273.py
result: 05-knowledge/results/lrc14_fiber_pgf_certificate_codex_s273.out
extends:
  - HYP-3136
  - HYP-3137
  - HYP-3135
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3129
  - HYP-3125
  - HYP-3124
  - HYP-3122
  - HYP-3112
related:
  - HYP-3078
  - HYP-3063
  - HYP-3009
  - HYP-2758
  - HYP-2523
  - OPEN-Q-108
external: generating functions, sheet-count PGF, Delsarte/MacWilliams transforms, q-Pochhammer tails, Moser-de Bruijn, fibbinary partial cubes, Lee-Yang
---

# HYP-3140: LRC14 Fiber-PGF Rprime Certificate

## Claim

For every residual LRC14 row of the form

```text
S = R union 14Q
```

the open multi-far coupling

```text
Rprime = meas(R-safe and Q-lonely)/(meas(R-safe)*meas(Q-lonely))
```

can be rewritten exactly as a finite generating-function moment.  Put
`u = 14t` and define

```text
N_R(u) = #{a in {0,...,13}: (u+a)/14 is R-safe}.
```

Because `Q-lonely` is a function of `u` alone,

```text
Rprime = E[N_R(u) | Q-lonely] / E[N_R(u)].
```

So the HYP-3129 signed-SPEC floor is the Fourier transform of a finite
14-sheet coefficient problem.  The next proof target is:

```text
E[N_R | Q-lonely] >= c * E[N_R]
```

for every legal residual row, with the HYP-3134 global-consistency sidecar
recording when individual sheets may be forgotten.

This is a child of HYP-3136's integrated multi-far factorization

```text
L(S) = Rprime(S) * meas(R-safe) * meas(Q-lonely).
```

HYP-3136 leaves the `Rprime` factor as the remaining signed-SPEC finite
constant chase; HYP-3140 replaces that one scalar factor by a finite
14-sheet coefficient/moment packet.

It also instantiates the completed HYP-3137 generating-function payload atlas:
the retained carrier is the coefficient and first-moment layer of
`F_R(y), F_R,Q(y)`, while raw scalar evaluation at `y=1` is deliberately not
the proof payload.

## Evidence

The scout `04-computation/lrc14_fiber_pgf_certificate_codex_s273.py` computes
the sheet-count PGF

```text
F_R(y) = sum_n meas{u: N_R(u)=n} y^n
```

and the Q-masked PGF

```text
F_R,Q(y) = sum_n meas{u: Q-lonely and N_R(u)=n} y^n.
```

It verifies, by exact rational interval arithmetic, that

```text
meas(R-safe) = E[N_R]/14
meas(Q-lonely) = meas(Q)
meas(R-safe and Q-lonely) = E[N_R 1_Q]/14
Rprime = E[N_R | Q]/E[N_R].
```

Representative exact readout:

| row | `Rprime` | `E[N_R]` | `E[N_R | Q]` | Q zero-sheet mass | Q sheet range |
|---|---:|---:|---:|---:|---:|
| drop 12, `Q={1}` | `1.16667` | `0.17023` | `0.19860` | `0.74286` | `0..2` |
| drop 13, `Q={1,2}` | `0.70147` | `0.47742` | `0.33489` | `0.52258` | `0..1` |
| `R={1..11}`, `Q={1,2,3}` | `0.96733` | `0.78867` | `0.76290` | `0.28925` | `0..3` |
| `R={1..10}`, `Q={1,2,3,4}` | `1.05790` | `1.93175` | `2.04359` | `0` | `1..4` |
| `R={1..9}`, `Q={1..5}` | `1.03957` | `2.53492` | `2.63522` | `0` | `1..5` |
| `R={1..8}`, `Q={1..6}` | `0.95869` | `3.72143` | `3.56771` | `0` | `2..7` |
| even `R`, `Q={1,3,5}` | `0.90074` | `6.40000` | `5.76471` | `0` | `4..12` |

The known worst targeted row becomes visibly finite and low-degree:

```text
R={1,...,12}, Q={1,2}
F_R(y)   = 7243/13860*y^0 + 6617/13860*y^1
F_R,Q(y) = 7243/13860*y^0 + 521/1980*y^1
Rprime   = 51058/72787 = 0.701471...
```

Thus minimum sheet count is too crude: `Q` still sees zero-sheet mass.  The
right theorem is a conditional first-moment lower bound on the two coefficients
`y^0,y^1`, not a pointwise positive-sheet theorem.

## Cross-Disciplinary Connections

- **Lee-Yang / PGF roots:** HYP-3112 keeps the miss-count PGF because roots
  move only after the one-runner ear payload is retained.  HYP-3140 supplies a
  simpler finite PGF for the same coupling: a sheet-count PGF over the 14
  lifts of the quotient coordinate.
- **Delsarte / MacWilliams / weight enumerators:** this is a transform problem
  between coefficient enumerators and Fourier/SPEC data.  The HYP-3129
  resonance lattice is the dual transform of `F_R,Q`.
- **q-Pochhammer / modular tails:** HYP-3078's lesson is that tail estimates
  become legal only after the principal coefficient packet is named.  Here the
  principal packet is `F_R(y), F_R,Q(y)`, and HYP-3129's L2 tail becomes a
  transform-side certificate.
- **Moser-de Bruijn / fibbinary / partial cubes:** HYP-3063 warns that
  automatic sequence membership is useful only as a typed carrier.  The 14
  sheet indicators are exactly such typed coordinates: they form a finite
  cube/partial-cube sidecar before quotienting to a scalar `Rprime`.
- **A000568 edge envelope:** HYP-3134 says local edge data may be forgotten
  only through a global-consistency quotient.  HYP-3140 gives the analogous
  analytic quotient: individual sheets may be forgotten only after the
  coefficient PGF and global-consistency class are recorded.

## Proof Target

Replace the open scalar floor

```text
Rprime >= c
```

with a coefficient theorem:

```text
For every legal residual packet (R,Q),
F_R,Q'(1) / F_R,Q(1) >= c * F_R'(1) / F_R(1).
```

Equivalently, `Q-lonely` cannot bias the 14-sheet count below a fixed
fraction of its global mean.  HYP-3129 already proves this numerically through
signed SPEC on the tested family; HYP-3140 identifies the finite generating
function whose coefficient inequalities would make the constant chase
symbolic.

## Assumption Challenge

Tournament/LRC vertices need not be runners.  In this lane the vertices are
`u`-fibers, 14 sheet coordinates, and PGF coefficients.  The preserved
predicate is the LRC14 residual floor `Rprime > 0`.  The destroyed information
is the identity of individual safe sheets.  The repair sidecars are
`fiber_pgf_coefficients`, `Q_masked_fiber_pgf`, `global_consistency_class`,
`SPEC_resonance_lattice_status`, and `edge_child_gluing_status`.

## Next Steps

1. Enumerate the legal finite residual packet family after HYP-3131 far-push
   reduction and compute its `F_R,Q` coefficient table.
2. Try to prove the worst `r=2` row by a two-coefficient inequality; it is the
   first row where `Q` zero-sheet mass is large enough to make pointwise
   positivity false.
3. Map HYP-3129's low-frequency SPEC terms to derivatives or coefficient
   moments of `F_R,Q`, so the existing `Rprime >= 0.64178` certificate becomes
   a symbolic finite coefficient bound.
4. Add `fiber_pgf_word` and `Q_masked_fiber_pgf_word` to HYP-3125 edge-floor
   packets before quotienting through HYP-3134.
