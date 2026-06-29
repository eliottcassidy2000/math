---
id: HYP-3441
title: LRC14 incident tax-core residue router
status: EVIDENCE / quotient bridge audit; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3437/HYP-3436/HYP-3311
tangent: T1402
technique: LTI-402
tournament_technique: LTT-302
script: 04-computation/lrc14_incident_tax_core_residue_router_codex_20260629.py
result: 05-knowledge/results/lrc14_incident_tax_core_residue_router_codex_20260629.out
reflection: 07-reflections/lrc14-incident-tax-core-residue-router-codex-20260629.md
related:
  - HYP-3455
  - HYP-3454
  - HYP-3453
  - HYP-3452
  - HYP-3451
  - HYP-3450
  - HYP-3440
  - HYP-3439
  - HYP-3438
  - HYP-3437
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3431
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3422
  - HYP-3418
  - HYP-3311
  - THM-523
  - OPEN-Q-108
---

# HYP-3441: LRC14 Incident Tax-Core Residue Router

## Claim

HYP-3437 is the strict one-branch overlap-tax theorem object: a rescue subset
must carry its own internal overlap tax.  In that strict model the audited
negative naive-slack rows have ranks:

```text
{2:7, 4:2, 5:48, 6:2}
```

HYP-3441 tests a deliberately different quotient.  An incident core collects
the full tax of every overlap atom it touches.  This is not a legal proof
replacement because it forgets which other owners make the touched atom
taxable.  It is a router for the missing bridge:

```text
strict HYP-3437 rescue core
  -> incident tax core
  -> endpoint atom bridge or HYP-3436/HYP-3438 survivor-gate word bridge
  -> HYP-3453 rank-2 gate transversal / HYP-3451 conductance route.
```

The useful fact is that all `59` negative rows collapse to incident core size
at most `2`.

## Exact Readout

Script:

```text
04-computation/lrc14_incident_tax_core_residue_router_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_incident_tax_core_residue_router_codex_20260629.out
```

Aggregate bridge audit:

```text
rows_audited=150
negative_naive_slack_rows=59
strict_rescue_rank_hist={2:7, 4:2, 5:48, 6:2}
incident_core_size_hist={1:10, 2:49}
strict_to_incident_rank_gap_hist={1:7, 3:49, 4:3}
strict_rank_larger_than_incident_size=59/59
unit_only_incident_cores=57/59
ramified_apex7_incident_cores=2/59
```

Residue and Galois sidecar signature:

```text
actual_core_hist_top=[((11,13),48), ((9,13),1), ((11,),1), ...]
residue_core_hist_top=[((11,13),48), ((7,),2), ((3,),2), ((1,),2), ...]
role_signature_hist={
  ('ramified_apex7',):2,
  ('unit_binding',):8,
  ('unit_binding','unit_binding'):49
}
C3_slot_signature_hist={
  ('apex7',):2,
  ('slot_pm1',):3,
  ('slot_pm3',):3,
  ('slot_pm3','slot_pm1'):48,
  ('slot_pm5',):2,
  ('slot_pm5','slot_pm1'):1
}
Qsqrt_minus7_signature_hist={
  ('qr_minus',):4,
  ('qr_plus',):4,
  ('qr_plus','qr_minus'):49,
  ('ramified_0',):2
}
```

Endpoint-only transfer is not enough:

```text
core_on_captured_atom_endpoint=58/59
all_core_owners_on_captured_atom_endpoint=51/59
endpoint_bounded_tax_covers_deficit=6/59
any_survivor_spine_core_attachment=10/59
best_survivor_spine_core_attachment=1/59
mixed_survivor_spine_core_attachment=1/59
both_branch_spine_core_attachment=2/59
```

Thus the next bridge cannot be a pure endpoint-tax theorem.  It must use the
HYP-3436/HYP-3438 bad-core and survivor-gate words.  After the concurrent
HYP-3453 mainline update, the strongest formulation is a component/gate
transversal: endpoint-tax failures should be routed through rank-`<=2`
survivor gates whenever a dead-cover obstruction is present.  HYP-3454 now
sharpens HYP-3452 into the AP-tail endpoint-clock clause for the canonical
tower, while HYP-3455 names the current noncanonical rank-`6` gate-gluing
clause as finite debt rather than an uncontrolled family.

## Canonical Stress Row

For `covering_AP_with_84`, HYP-3437 needs strict rank `6`:

```text
strict_subset=(3,5,7,9,11,13)
deficit=18586/315315
overlap_tax=4055/63063
```

The incident quotient finds:

```text
incident_core=(11,13)
incident_tax=3769/63063
incident_margin=37/45045
tax/deficit=1.013935
```

But endpoint-bounded tax gives only `0.908143` of the deficit.  That is the
proof-facing obstruction in one line: the AP rank-6 core is mostly a two-owner
residue/incident packet, but the legal transfer needs the canonical
corridor-fence or survivor-gate bridge.

## Interpretation

The user's factorization prompt is supported in a precise but limited form.

```text
14 = 2 * 7
7-adic residue layer: unit binding pairs plus ramified apex-7 exceptions
2-adic magnitude layer: endpoint/survivor bridge and boundary exits
```

The `C3` quotient organizes unit binding pairs.  The `Qsqrt(-7)` quadratic
character is transverse: the dominant two-owner incident core has one
`qr_plus` and one `qr_minus` owner.  The two ramified-apex rows are explicit
covering-layer exceptions, not the main proof engine.  The `12->24` hinge
remains a two-adic boundary issue, not a 7-adic census proof.

Euler-Mascheroni, Meissel-Mertens, Ramanujan-Soldner, Sophie Germain,
Krasner, and Hermite-Lindemann-Weierstrass are retained only as sidecar
discipline: priority budgets, balance-root metaphors, quartic height splits,
p-adic branch stability, and no-free-slider guards.  None can replace strict
overlap atoms, endpoint labels, or survivor-gate words.

## Candidate Bridge Lemma

For every primitive covering row with negative one-branch naive slack, either:

1. the strict HYP-3437 rescue subset already has bounded rank and is
   discharged directly;
2. the incident core of size at most `2` transfers legally to strict tax
   through endpoint atoms plus HYP-3436/HYP-3438 survivor-gate words, with
   HYP-3453 as the current rank-`<=2` gate-transversal candidate;
3. the row is canonical corridor-fence and routes to HYP-3431;
4. the row is the finite noncanonical gate-gluing clause isolated by HYP-3455;
5. the row emits named owner-current, exact-period, state-lift, signed-SPEC,
   two-adic loss, or new residual debt.

## Tournament Analysis

Vertices are proof obligations and quotient carriers, not runners or arcs.

```text
pairwise_observable =
  strict predicate retained + incident compression + lost-coordinate debt
switch_gauge =
  higher weighted proof-facing score; ties by declared carrier order
score_hist={13:1, 40:1, 52:1, 58:1, 59:1, 60:2, 61:1}
directed_3cycles=0
hamiltonian_path_count=1
hamiltonian_path =
  H3436_survivor_gate_word_bridge
  -> endpoint_atom_bridge_obligation
  -> two_adic_boundary_exit
  -> C3_Qsqrt_minus7_residue_packet
  -> incident_tax_core_router
  -> strict_H3437_overlap_core
  -> Mertens_gamma_priority_budget
  -> raw_constant_or_scalar_route
```

Assumption challenge: runners, odd blockers, even-half gates, fixed circle
sections, section boundaries, wall crossings, residues, cover arcs, Fourier
modes, matroid circuits, strict overlap subsets, incident tax owners, endpoint
labels, survivor-gate words, and proof obligations were considered.  The
chosen quotient preserves the deficit/tax comparison as an inequality target
but destroys internal co-owner multiplicity.
