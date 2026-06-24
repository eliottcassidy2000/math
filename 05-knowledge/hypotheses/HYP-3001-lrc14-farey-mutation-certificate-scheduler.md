---
id: HYP-3001
title: LRC14 Farey mutation certificate scheduler
status: EVIDENCE / proof-route scheduler; not a proof
source: codex-2026-06-24-S164
tangent: T1086
script: 04-computation/lrc14_farey_mutation_certificate_scheduler_codex_s164.py
result: 05-knowledge/results/lrc14_farey_mutation_certificate_scheduler_codex_s164.out
related:
  - HYP-2931
  - HYP-2932
  - HYP-2940
  - HYP-2947
  - HYP-2974
  - HYP-2978
  - HYP-2979
  - HYP-2981
  - HYP-2982
  - HYP-2983
  - HYP-2908
  - THM-572
  - OPEN-Q-108
---

# HYP-3001: LRC14 Farey Mutation Certificate Scheduler

## Claim

The prompt's four Farey mutations are useful only after the exact Farey address is retained.

For `M(S)=p/q`, keep

```text
e = 14p - q,       M(S)-1/14 = e/(14q).
```

Then read the four literal mutations as:

```text
product-value      (p*q)/q = p
sum-value          (p+q)/q = 1 + p/q
denominator-power  p/(q^p)
numerator-power    (p^q)/q
```

The new useful reframe is:

```text
After e=1, product-value p is a certificate scheduler.
```

It is not an ordering scalar and not a theorem denominator.  But on the unit-excess layer it routes:

```text
p=1      q-parent / right-neighbor branch
p=2      C=27 petal or two-block branch
p>=3     K33 / state-lift / Fejer packet branch
```

This integrates HYP-2931/HYP-2940 with the newer HYP-2981/HYP-2982/HYP-2983 stack: Farey mutations are front-end schedulers, while Fejer/Ramanujan/Kaczynski packets are the certificate back end.

## Computation

Script:

```text
04-computation/lrc14_farey_mutation_certificate_scheduler_codex_s164.py
```

Stored output:

```text
05-knowledge/results/lrc14_farey_mutation_certificate_scheduler_codex_s164.out
```

The script reuses the S130 AP/GW/petal/single-replacement bank (`749` rows) and selected HYP-2981 Fejer packet rows.

## Low-Frontier Result

In the S130 bank:

```text
M <= 3/41:
  AP, GW 12->24                    p=1, boundary
  near-miss 12->36                 p=3, K33/state-lift

M <= 2/27:
  AP, GW 12->24                    p=1, boundary
  swap 10->20, swap 13->26         p=2, C27 petal/two-block
  near-miss 12->36                 p=3, K33/state-lift
```

Thus the low known frontier is cleanly scheduled by `(e,p)`:

```text
e=0            AP/GW boundary
e=1,p=2        C27 petal/two-block branch
e=1,p=3        K33 branch
```

## Negative Guardrail

The literal product-value `p` is not order-safe on the full row bank.  Pairwise order agreement against the true key `M-1/14`:

```text
product_value       agree= 88234 flip= 71462 tie=120430 score=21207/40018
sum_value           agree=258329 flip=     0 tie=21797  score=538455/560252
denominator_power   agree=170095 flip= 88234 tie=21797  score=361987/560252
numerator_power     agree=148543 flip=109786 tie=21797  score=318883/560252
```

So:

```text
sum-value = theorem-safe affine check, because it is 1+M;
product-value = branch scheduler only after e=1;
power mutations = magnitude stress tests, not proof denominators.
```

## Selected Packet Rows

The selected HYP-2981 Fejer rows land as expected:

```text
near/K33 12->36                       M=3/41,  e=1,  p=3   -> K33/state-lift
P10+GW and P10+K33                    M=2/27,  e=1,  p=2   -> C27/two-block
single swap 6->63                     M=2/25,  e=3,  p=2   -> nonunit excess
covering 12->168                      M=14/173,e=23, p=14  -> nonunit excess
two drop(12,13)->add(14,29)           M=1/12,  e=2,  p=1   -> nonunit excess
```

This is the useful split: `p` alone does not classify nonunit-excess packets, but it exactly identifies the known unit-excess proof route.

## Tournament Analysis

Vertices are the four literal mutations:

```text
sum_value
product_value
numerator_power
denominator_power
```

Pairwise observable: retained proof-role vector

```text
order safety,
branch scheduling,
packet side-channel,
magnitude stress,
false-quotient detection.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1}
directed_3cycles=0
Hamiltonian path:
  sum_value > product_value > numerator_power > denominator_power
```

Alternate vertices considered: runners, Farey fractions, denominators, endpoint owners, Fejer packet fibers, Ramanujan phases, and proof obligations.  The chosen quotient preserves the binding gap only after `M=p/q` and `e=14p-q` are retained.  It destroys wall geometry, so packet labels must be restored.

## Proof Route

Current best use:

```text
1. Keep exact M=p/q and e=14p-q first.
2. If e=0, route to AP/GW boundary equality.
3. If e=1, use product-value p as scheduler:
   p=1 -> q-parent/right-neighbor;
   p=2 -> C27 petal/two-block;
   p>=3 -> K33/state-lift/Fejer.
4. If e>1, do not trust p; route through q-threshold, Fejer/Toeplitz,
   Ramanujan/exact-period, Kaczynski boundary, or state-lift packets.
5. Use power mutations only to detect magnitude-blind false quotients.
```

LRC14 is not proved.  The contribution is a sharper dispatch rule for the labelled packet proof stack.
