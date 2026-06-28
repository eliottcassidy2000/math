---
id: HYP-3432
title: LRC14 Euler-Mascheroni harmonic wall-budget sidecar
status: SYNTHESIS / exact rational harmonic sidecar audit; not an LRC14 proof
source: codex-2026-06-28 continuation of HYP-3430/HYP-3429/HYP-3428/HYP-3427 under the Euler-Mascheroni prompt
tangent: T1393
technique: LTI-393
tournament_technique: LTT-293
script: 04-computation/lrc14_euler_mascheroni_wall_budget_codex_20260628.py
result: 05-knowledge/results/lrc14_euler_mascheroni_wall_budget_codex_20260628.out
reflection: 07-reflections/lrc14-euler-mascheroni-wall-budget-codex-20260628.md
related:
  - HYP-3430
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3419
  - HYP-3418
  - HYP-3417
  - HYP-3415
  - HYP-3412
  - HYP-3407
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3432: LRC14 Euler-Mascheroni Harmonic Wall-Budget Sidecar

## Claim

Euler's constant should not be imported as a magic constant into the LRC14
floor.  Its useful role is structural: `H_N - log N` says that a reciprocal
tail can be renormalized, but the finite LRC14 proof still needs an exact
survivor interval, a branch, and endpoint wall labels.

HYP-3430 already proves the finite intercept `H_N - log N` is a scalar
firewall: it tracks scale but does not determine endpoint-spine class.
HYP-3432 is the complementary exact budget audit.  Attach to every HYP-3427
wall window and HYP-3429 endpoint spine the exact rational budget

```text
sum_{endpoint owner v} 1/v.
```

This harmonic/Mertens-like budget can rank endpoint debt, but it cannot replace
the HYP-3427/HYP-3429 certificate.  The exact audit shows why: scalar harmonic
budgets collide massively across distinct wall signatures.

## Exact Readout

Script:

```text
04-computation/lrc14_euler_mascheroni_wall_budget_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_euler_mascheroni_wall_budget_codex_20260628.out
```

The script uses only exact rational arithmetic.  No floating logarithm or
floating Euler-Mascheroni value is used.

Finite normalizers:

```text
H_14 = 1171733/360360
reciprocal budget of 2,7,14 = 5/7
reciprocal budget of units mod 14 = 11662/6435
```

On the HYP-3429 best endpoint spines:

```text
rows audited = 150
distinct best-spine harmonic budgets = 144
best-spine kind histogram =
  {('B0',):1, ('B0','B1'):2, ('B0','E'):14,
   ('B1',):3, ('B1','E'):14, ('E',):116}
scalar budget shape-collisions = 2
minimum best-spine budget = 1/4032 at canonical_84m_ext_48, label E:4032
maximum best-spine budget = 89/420 at covering_AP_with_84, labels B1:5,E:84
minimum endpoint reciprocal share = 715/8928419
maximum endpoint reciprocal share = 137847056256/411769125631
```

On the full HYP-3427 wall-window atlas:

```text
survivor windows audited = 5524
distinct wall-owner budgets = 1291
budgets with more than one wall signature = 1197
```

The tight row `{1..11,13,84}` already shows the lesson in miniature:

```text
W1: O1:7 + E:84 budget 13/84, branch b1
W2: E:84 + O1:5 budget 89/420, branch b1
W3: O0:5 + E:84 budget 89/420, branch b0
W4: E:84 + O0:7 budget 13/84, branch b0
```

The harmonic scalar sees `13/84` and `89/420`; the proof still needs branch
orientation and the `O0/O1/E` wall word.

## Finite Lemma Target

Use the harmonic budget as a priority sidecar in the endpoint-spine lemma:

```text
Every primitive covering row has either an E-only free component, a rank-2
mixed endpoint spine, or a named exception; among competing certificates,
the reciprocal endpoint budget ranks which wall-owner debt should be tried
first, but the certificate is accepted only with its exact interval, branch,
and endpoint labels.
```

This is the Euler-Mascheroni-compatible reframe: renormalized harmonic tails
may organize the search for finite wall debt, while exact endpoint geometry
does the proving.

## Tournament Analysis

Tournament vertices are proof carriers and tail sidecars, not runners.

```text
axes =
  predicate_retention, wall_exactness, tail_budget_signal,
  two_adic_compatibility, collision_safety, proof_readiness
score_hist = {8:1, 20:1, 33:1, 39:1, 51:1, 53:1, 56:1}
directed_3cycles = 0
hamiltonian_path =
  endpoint_wall_certificate
  -> component_spine_rank2
  -> two_adic_loss_ledger
  -> harmonic_owner_budget_sidecar
  -> mertens_reciprocal_tail_rank
  -> gamma_shadow_normalizer
  -> named_constant_scalar_shortcut
```

## Assumption Challenge

Alternate vertices considered included runners, gaps, fixed circle sections,
section boundaries, wall-crossing events, endpoint owners, reciprocal budgets,
Mertens tails, gamma shadows, component spines, and proof obligations.  The
chosen quotient preserves the LRC predicate only when paired with an actual
survivor interval and branch.  Used alone, it destroys branch orientation,
endpoint type, interval width, and the exact wall word.
