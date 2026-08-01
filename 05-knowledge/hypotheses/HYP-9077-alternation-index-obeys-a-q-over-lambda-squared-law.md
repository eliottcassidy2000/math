---
id: HYP-9077
title: "The alternation index of a binary linear recurrence is governed by Q/lambda^2 -- and golden is not maximal"
status: >
  VERIFIED-EXACT (numeric, 7 recurrences). For a_k = P a_{k-1} + Q a_{k-2} the
  Simson identity gives a_{k-1}a_{k+1} - a_k^2 = (-Q)^{k-1} c, so the Newton
  circuit alternates in sign with magnitude decaying at rate Q/lambda^2,
  lambda = (P+sqrt(P^2+4Q))/2. The alternation index of THM-3009 sec 12.2 is
  MONOTONE INCREASING in Q/lambda^2. Consequence: within the metallic family
  Q = 1 the golden ratio maximises it (since lambda is least there), but over
  ALL binary recurrences it does NOT -- the DLO-counting sequence
  a_k = 2a_{k-1} + 5a_{k-2} of the Z-set-theory problem beats it, and (1,3)
  beats that. Corrects the 'golden maximal' claim of THM-3009 sec 12.2, whose
  scope was the metallic family only.
source: opus-2026-07-31-amm12592-writeup
related:
  - THM-3010  # klein: metallic recurrences attain maximal circuit alternation
  - THM-3009  # sec 12.2, where the index was defined and first tested
  - THM-3004  # circuit sign-change cluster law
script: 04-computation/alternation_index_q_over_lambda2_law.py
output: 05-knowledge/results/dlo_count_and_alternation_law.out
---

# HYP-9077 -- the Q/lambda^2 law

## 1. Statement

For `a_k = P a_(k-1) + Q a_(k-2)` the Simson/Catalan identity reads
`a_(k-1)a_(k+1) - a_k^2 = (-Q)^(k-1) c`, so with `lambda = (P+sqrt(P^2+4Q))/2`

```text
R_k - 1 ~ +- (Q/lambda^2)^k :
the circuit ALTERNATES in sign and its MAGNITUDE decays at rate Q/lambda^2.
```

Since the alternation index measures alternation *with sustained magnitude*
(THM-3009 sec 12.2), it should be an increasing function of `Q/lambda^2`.
Measured (window `N = 16`):

```text
 (P,Q)    lambda    Q/lambda^2    index
 (3,1)   3.3028     0.091673     0.570654    bronze
 (2,1)   2.4142     0.171573     0.607902    silver
 (1,1)   1.6180     0.381966     0.700362    GOLDEN
 (2,5)   3.4495     0.420204     0.723772    DLO count
 (2,7)   3.8284     0.477592     0.746447
 (1,2)   2.0000     0.500000     0.756588
 (1,3)   2.3028     0.565741     0.785948
```

Monotone in `Q/lambda^2` across all seven.

## 2. What this corrects

THM-3009 sec 12.2 reported that the index is a strict TIE-BREAKER inside
klein's maximal-alternation class with **golden maximal**. That is true only
inside the METALLIC family `Q = 1`, where `Q/lambda^2 = 1/lambda^2` and
`lambda` is least at `phi`. Over all binary recurrences golden is not
maximal, and the counterexample arrived from outside the repo: the
DLO-counting sequence `a_k = 2a_(k-1) + 5a_(k-2)` (growth `1+sqrt6`) of the
Z-set-theory problem has index `0.7238 > 0.7004`.

## 3. Caveat

The index depends on the window length: silver reads `0.590` at `N = 14` and
`0.608` at `N = 16`. It should always be quoted with `N`. The ORDERING was
stable across both windows tested.

## 4. Why it matters here

THM-3010's criterion (maximal sign-change count) is degenerate on this whole
family -- every one of these attains it. The index resolves them, and the
`Q/lambda^2` law says exactly what it resolves them by. Whether that
quantity has a meaning in klein's Newton-circuit setting (where `Q` is not a
free parameter but comes from a ballot column) is untested.
