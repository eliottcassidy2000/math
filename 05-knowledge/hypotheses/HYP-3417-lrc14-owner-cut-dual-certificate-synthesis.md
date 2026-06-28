---
id: HYP-3417
title: LRC14 owner-cut dual certificate synthesis
status: SYNTHESIS / exact owner-current certificate scout; not an LRC14 proof
source: codex-2026-06-28 synthesis of HYP-3408, HYP-3409, and HYP-3410
tangent: T1378
technique: LTI-378
tournament_technique: LTT-278
script: 04-computation/lrc14_owner_cut_dual_certificate_synthesis_codex_20260628.py
result: 05-knowledge/results/lrc14_owner_cut_dual_certificate_synthesis_codex_20260628.out
reflection: 07-reflections/lrc14-owner-cut-dual-certificate-synthesis-codex-20260628.md
related:
  - HYP-3416
  - HYP-3415
  - HYP-3414
  - HYP-3413
  - HYP-3412
  - HYP-3411
  - HYP-3410
  - HYP-3409
  - HYP-3408
  - HYP-3407
  - HYP-3406
  - HYP-3405
  - HYP-3404
  - HYP-3402
  - HYP-3311
  - HYP-3310
  - HYP-3265
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3417: LRC14 Owner-Cut Dual Certificate Synthesis

## Claim

The strongest current proof-facing obligation is no longer just "find the
owner cut."  HYP-3410 already found small owner cuts.  The next obligation is:

```text
Every surviving mixed residue/height fiber should have a labelled
owner-current certificate of bounded size and positive margin, or emit named
owner/height/off-grid/state-lift debt.
```

The certificate is a zero-level signed current on endpoint-owner labels:

```text
unit-island current:
  a label hits the unit-petal row and no positive-open row

positive-debt current:
  a small label set avoids the unit-petal row and hits every positive-open row
```

The named prompt ideas become audit columns around this finite certificate.
Krasner means the common owner core is not a stable theorem-exit root packet.
Sophie-Germain splits the cut-size/core-size obstruction into two exact
quadratic channels.  Meissel-Mertens chooses among equal-size cuts by
reciprocal owner budget.  Ramanujan-Soldner declares the zero level of the
current.  HLW is the guardrail forbidding replacement of labelled finite
packets by scalar shadows.

## Exact Readout

Script:

```text
04-computation/lrc14_owner_cut_dual_certificate_synthesis_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_owner_cut_dual_certificate_synthesis_codex_20260628.out
```

Over the three exact HYP-3410 mixed fibers:

```text
height_leak_12_family:
  common core = {6:g2, 7:g7}
  selected current = positive debt {5:g1}
  margin = 1
  Sophie channels (core,cut)=(2,1): plus=10, minus=2

persistent_owner_leak_26_40_54_family:
  common core = {12:g2}
  selected current = unit island {1:g1}
  margin = 1
  positive-debt current = {11:g1, 13:g1}
  positive-debt Sophie channels (1,2): plus=13, minus=5

height_persistent_owner_leak_10_20_drop_add_family:
  common core = {6:g2}
  selected current = positive debt {2:g2, 11:g1, 13:g1}
  margin = 1
  Sophie channels (1,3): plus=25, minus=13
```

Thus all selected certificates on the current substrate have positive margin
`1`, and the largest selected cut has size `3`.

## The New `13`-Channel Signal

The most useful niche connection is in the positive-debt orientation:

```text
(core, cut) = (1,2) -> Sophie channels 13 and 5
(core, cut) = (1,3) -> Sophie channels 25 and 13
```

In the persistent owner leak, the Mertens-cheapest positive-debt cut is

```text
{11:g1, 13:g1}
```

In the `(72,20)` `10->20` frontier, the Mertens-cheapest positive-debt cut is

```text
{2:g2, 11:g1, 13:g1}
```

This is exactly "one even-cover label plus two binding labels."  The same
frontier has top finite-BDH variance `13:g1=49/50`, and the Sophie channel
also emits `13`.

This is not a proof.  It is an exact finite signal saying that the next theorem
should keep the labelled owner packet while testing whether a repeated
`13`-channel quadratic factor is the algebraic shadow of the first bounded
owner-current cut.

Rebase integration with incoming S257: HYP-3411/HYP-3413 says the global
residue/magnitude story is controlled by the `C3` residue gate and the
Goddyn-Wong `q == 1 mod 3` criterion.  HYP-3417 is narrower.  It shows that the
current frontier owner cut has the same "one even-cover magnitude label plus
binding labels" shape, but it remains a local certificate target until the
bounded-current theorem is proved on the enlarged bank.

S258/HYP-3415 is a stronger warning: the critical path for closing LRC14 is the
decorrelation floor inequality.  HYP-3417 should therefore be used as a sidecar
certificate for mixed owner fibers, not as a replacement for that floor theorem.
Incoming HYP-3416 supplies the recursive quotient-ladder context where this
sidecar can live: the owner current is one certificate layer in the ladder, not
the ladder itself.

## Krasner / Owner-Root Gate

For all three fibers:

```text
krasner_core_projection_exit_pure = False
```

The labels common to every row in a mixed fiber do not determine the theorem
exit.  Therefore a Krasner-style local stability statement cannot stabilize
only the common p-adic or owner core.  The first non-core owner current is the
root-instability payload.

## Proof Target

The next theorem shape is:

```text
For every mixed residue/height fiber in the enlarged HYP-2963 bank, either:
  1. a unit-island owner current separates unit-petal exit from open exit;
  2. a positive-debt owner hitting current separates open exit from unit-petal;
  3. the common owner/contact core becomes theorem-exit pure after a legal
     charal recursion step;
  4. the fiber routes to AP/GW, strict-open mass, q-witness, K33/H7/state-lift,
     off-grid floor, exact-period/BDH exception, or named debt.
```

If the positive-debt current stays bounded by `3` past the current frontier,
the owner-support repair in HYP-3406/HYP-3410 may be promotable from an
observed sidecar into a finite Menger/Farkas owner-current lemma.

## Tournament Analysis

Vertices are proof obligations and owner-current certificates, not runners or
raw constants.

Pairwise observable:

```text
exact mixed-fiber discharge plus retained labelled sidecars
```

Switch/gauge:

```text
higher weighted certificate score; ties by declared priority
```

Fingerprint:

```text
vertices = 7
score_hist = {-23:1, 20:1, 25:1, 33:1, 34:1, 57:1, 59:1}
directed_3cycles = 0
hamiltonian_path_count = 1
priority_path =
  owner_cut_dual_current_certificate
  -> positive_debt_hitting_set_theorem
  -> krasner_owner_core_instability_gate
  -> sophie_germain_channel_audit
  -> finite_mertens_budget_selector
  -> c3_7adic_2adic_group_split
  -> raw_named_constant_scalar
```

## Assumption Challenge

Alternative vertices considered:

```text
runners, residues, owner labels, owner cuts, signed currents, Sophie channels,
Mertens budgets, Krasner cores, and proof obligations.
```

The chosen vertices are proof obligations and labelled owner-current
certificates.  This quotient preserves theorem-exit separation and destroys
raw row order.  The challenged assumption is that a famous scalar, raw
p-adic closeness, or an unlabelled Sophie channel can replace the endpoint
owner packet.
