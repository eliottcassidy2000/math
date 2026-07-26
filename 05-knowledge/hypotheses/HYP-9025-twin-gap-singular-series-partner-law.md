---
id: HYP-9025
title: "Twin-gap singular-series partner law"
status: OPEN (quantitative empirical law, verified shape at centers <= 10^8)
source: kind-pasteur-2026-07-26-S131
related:
  - THM-2443-twin-rank-parent-parity-margins-and-boundary-crossing-bridge
  - THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry
  - HYP-1994-twin-goldbach-necklace
script: 04-computation/twin_rank_parity_margins_bridge_thm2443.py
output: 05-knowledge/results/twin_rank_parity_margins_bridge_thm2443.out
---

# HYP-9025 -- which previous terms combine: the singular-series partner law

The owner's A014574 puzzle asks for "the pattern of which previous
terms combine to make each new term". THM-2422 (30) shows the
immediate-predecessor rule is exactly the self-gap condition
`g_i = k_i - k_{i-1} in K`, holding on `185,644 / 440,310` census
transitions. This hypothesis is the quantitative law behind the
remaining structure: **the frequency of a consecutive-rank gap `g`
is a smooth decaying availability profile times an explicitly
computable Hardy-Littlewood singular-series weight**

```text
w(g) = prod_{p >= 5, p | g} (p-2)/(p-4)
     * prod_{p >= 5, p | 9g^2-1} (p-3)/(p-4).                  (1)
```

Derivation sketch: consecutive ranks `m` and `m+g` require the four
numbers `6m +- 1, 6(m+g) +- 1` prime. For `p >= 5` the generic number
of forbidden residues of `m mod p` is `4`; it drops to `2` when
`p | g` (the two pairs of forbidden classes coincide) and to `3` when
`p | 3g-1` or `p | 3g+1` (one cross-coincidence; jointly
`p | 9g^2-1`). Ratioing against the generic case gives (1).

Evidence (census at centers `<= 10^8`, output lines 33-73): the raw
gap counts oscillate wildly (e.g. `4768, 12568, 9153, 5730, 17563`
for `g = 1..5`), while the normalized column `count/w(g)` is smooth
and monotone-decaying to the resolution of the table
(`4768.0, 4713.0, 4576.5, 4512.4, 4390.8, ...`). In particular the
partner ranking is **smoothness-ranked, not size-ranked**: the
dominant "other summands" in raw-center units are `30, 42, 12, 60,
18, 72, 180, ...` -- precisely the primorial-adjacent centers whose
gap weights (1) are largest. This is the same arithmetic-extremizer
phenomenon as the repo's standing mistake-card ("random sampling
misses arithmetic extremizers"), realized inside a natural sequence.

Predictions to test:

1. `count(g) / w(g) = A(x) * exp(-c(x) * W(g)) * (1 + o(1))` for a
   Cramer-style cumulative availability `W(g) = sum_{g' <= g} w(g')`,
   uniformly in dyadic center windows (the residual dips at
   `g = 9, 19, 24, 29, 31, 34, 36` in the flat-normalization table
   should be absorbed by the cumulative form).
2. The fraction of transitions satisfying the local rule (30)
   (`~ 42.2%` at `10^8`, window-decaying `0.55 -> 0.43` between
   `10^6` and `10^7`) tends to `0` slowly, since the gap distribution
   shifts right as `log^2 x` while `K` thins as `1/log^2 x`; the
   decay rate should follow prediction 1.
3. The repair-depth distribution recorded in THM-2422's companion
   (`0:185643, 1:67517, 2:41034, ..., max 123`) is asymptotically
   geometric with window-dependent parameter
   `q(x) = P(gap value in K near x)`, with the fat tail at
   `k = 1133072, 6803528, 13153018` explained by locally
   singular-series-poor gap runs, not by any new obstruction.

Cheapest decisive test: recompute the gap table in disjoint dyadic
windows of `k` and fit prediction 1; a single window violating the
cumulative form beyond sampling error refutes the law as stated.
