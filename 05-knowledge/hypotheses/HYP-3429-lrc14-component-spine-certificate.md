---
id: HYP-3429
title: LRC14 component-spine endpoint certificate for the two-adic floor
status: EVIDENCE / exact endpoint-spine audit; not an LRC14 proof
source: codex-2026-06-28 extension of HYP-3428/HYP-3427/HYP-3425 after HYP-3426
tangent: T1390
technique: LTI-390
tournament_technique: LTT-290
script: 04-computation/lrc14_component_spine_certificate_codex_20260628.py
result: 05-knowledge/results/lrc14_component_spine_certificate_codex_20260628.out
reflection: 07-reflections/lrc14-component-spine-certificate-codex-20260628.md
related:
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3420
  - HYP-3419
  - HYP-3418
  - HYP-3417
  - HYP-3415
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3429: LRC14 Component-Spine Endpoint Certificate

## Claim

HYP-3425 reduced the two-adic covering-floor relocation target to

```text
relocation_good = E_safe minus (B0_odd cap B1_odd).
```

HYP-3426 removes branch ambiguity by the mirror involution `u -> 1-u`.
HYP-3427 records the branch mask and exact wall signature of survivor windows.
HYP-3428 records the two-adic descent loss ledger and exception vocabulary.
HYP-3429 pushes the next proof compression: in exact audits, positive
`relocation_good` is not merely a surviving mass or wall word.  It contains a
survivor window whose endpoints are controlled by at most two active walls:

```text
E  = even-safe wall, ||e*u|| = 1/14
B0 = branch-0 odd wall, ||o*u/2|| = 1/14
B1 = branch-1 odd wall, ||o*u/2|| = 3/7.
```

Thus the proof target can be stated as an endpoint-spine lemma rather than an
unstructured interval-union estimate.

## Exact Readout

Script:

```text
04-computation/lrc14_component_spine_certificate_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_component_spine_certificate_codex_20260628.out
```

The script audits `150` rows: the `62` HYP-3425 rows, the canonical
`{1..11,13,84m}` tower extended through `m=48`, two-tail extensions through
`m=18`, and forty deterministic primitive covering rows.

Aggregate readout:

```text
positive survivor window rows:        150/150
rows with all windows endpoint labels: 150/150
endpoint-labelled survivor windows:   15576/15576
best endpoint-spine rank <= 2:         150/150
mixed even/odd endpoint spine exists: 148/150
both-branch survivor exists:          149/150
best-rank histogram:                  {1: 47, 2: 103}
```

The smallest best spine in this bank is still explicit:

```text
row = canonical_84m_ext_48
best window = [4621/28224, 4633/28224]
length = 1/2352
branches = (1,)
labels = E:4032
```

The two rows without a mixed even/odd spine are not failures.  They are easier
E-only windows:

```text
canonical_84m_ext_35: E:2940, length 1/1715
two_tail_ext_15: E:1260 and E:1470, length 71/61740
```

## Candidate Lemma

For every primitive covering row `S = O union 2E`, at least one component of

```text
E_safe minus (B0_odd cap B1_odd)
```

has an endpoint-spine certificate of one of the following forms:

1. an E-only free component whose endpoints are even-safe walls and whose
   interior avoids the two-color odd bad core;
2. a mixed `E/B0` or `E/B1` endpoint spine of rank at most `2`;
3. a two-color odd wall `B0/B1` spine of rank at most `2`, followed by a
   branch-both witness;
4. a named owner-current or exact-period exception.

This lemma would turn the HYP-3425 Helly theorem into a finite endpoint
certificate problem.  It is stronger than a mass lower bound because it
retains the actual branch choice and the active walls that create the witness.

## Proof Use

The component-spine audit suggests the covering floor may admit the following
local proof route:

```text
two-adic split S=O union 2E
-> E_safe components
-> endpoint spine of rank <= 2
-> branch t=u/2 or t=(u+1)/2
-> signed-SPEC / Rprime handoff if the endpoint route becomes analytic.
```

The key move is to avoid scalarizing `good_union` into total measure.  The
proof should show that at least one even-safe component has an endpoint
configuration that cannot be swallowed by the two-color odd core.

## Guardrail

HYP-3429 does not claim LRC14.  It also does not replace HYP-3421 Rprime,
HYP-3422 relocation, HYP-3423 topology guardrail, HYP-3424 transfer, or
HYP-3426 mirror reduction, and it should be read as a rank-compressed
companion to HYP-3428's loss ledger and HYP-3427's wall-signature atlas.  It
only sharpens the finite-ruler part of the two-adic branch.  Owner-current
labels remain names for exceptions, not the floor proof itself.

## Tournament Analysis

Tournament vertices are proof carriers/endpoints, not runners or raw
components.

```text
score_hist={26:1, 51:1, 57:1, 58:2, 60:1, 63:1}
directed_3cycles=0
hamiltonian_path=
  two_endpoint_spine_certificate
  -> mixed_even_odd_wall_certificate
  -> even_only_free_component
  -> helly_pair_piercing_bound
  -> branch_both_relocation_window
  -> owner_current_exception_router
  -> raw_good_measure_scalar
```

## Assumption Challenge

Alternate vertices considered:

```text
rows, runners, E_safe components, survivor windows, interval endpoints,
endpoint labels, odd-pair walls, branch choices, proof obligations.
```

The chosen compression uses endpoint labels and proof carriers.  This preserves
the LRC covering-floor predicate only because each spine still selects an
actual `u`-window and a legal branch `t`.
