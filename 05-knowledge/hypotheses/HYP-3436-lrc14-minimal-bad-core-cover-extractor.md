---
id: HYP-3436
title: LRC14 minimal bad-core cover extractor
status: EVIDENCE / exact interval-cover classification; not a proof
source: codex-2026-06-29 continuation of HYP-3435 and the HYP-3422/HYP-3425 two-adic bad-core chain
tangent: T1397
technique: LTI-397
tournament_technique: LTT-297
script: 04-computation/lrc14_minimal_bad_core_cover_extractor_codex_20260629.py
result: 05-knowledge/results/lrc14_minimal_bad_core_cover_extractor_codex_20260629.out
reflection: 07-reflections/lrc14-minimal-bad-core-cover-extractor-codex-20260629.md
related:
  - HYP-3437
  - HYP-3435
  - HYP-3434
  - HYP-3433
  - HYP-3432
  - HYP-3431
  - HYP-3430
  - HYP-3429
  - HYP-3428
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3436: LRC14 Minimal Bad-Core Cover Extractor

## Claim

HYP-3436 executes the inversion requested by HYP-3435.  Instead of first
selecting a surviving branch witness in

```text
E_safe(1/14) cap (branch0_good union branch1_good),
```

it classifies the complementary two-color bad core

```text
E_safe(1/14) cap B0_odd cap B1_odd.
```

For every bad-core component in the audited primitive covering bank, the
extractor records:

```text
minimal branch-0 odd-owner cover
minimal branch-1 odd-owner cover
left/right endpoint wall labels
parent E_safe component status
survivor-gap data
```

The audit finds survivors in every row, so it does not produce an LRC14 proof.
It does make the next finite lemma sharper: a counterexample would have to glue
local minimal bad-core covers across all even-safe components while retaining
the branch owner words and even endpoint gates.  Raw bad-core measure,
Euler-Mascheroni/Mertens tails, and topology-only summaries are illegal
shortcuts unless they reconstruct or route those sidecars.

## Exact Readout

Script:

```text
04-computation/lrc14_minimal_bad_core_cover_extractor_codex_20260629.py
```

Stored result:

```text
05-knowledge/results/lrc14_minimal_bad_core_cover_extractor_codex_20260629.out
```

Aggregate audit:

```text
rows_audited=135
structured_rows=15
random_rows=120
rows_with_survivor=135/135
total_E_components=17164
E_components_fully_bad=3770
E_components_mixed=6228
E_components_clean=7166
bad_core_components=11670
survivor_components=15868
max_cover_total=6
max_endpoint_support_size=3
```

Minimal bad-core cover histogram:

```text
(1,1):10288, (1,2):531, (2,1):531, (2,2):184,
(1,3):35, (3,1):35, (2,3):22, (3,2):22,
(3,3):6, (1,4):5, (4,1):5, (1,5):2, (5,1):2,
(2,4):1, (4,2):1
```

Here `(a,b)` means `a` branch-0 odd owners and `b` branch-1 odd owners are
needed to cover the component.  Most local bad-core pieces are single-owner on
both branches, but the hardest random components need six total odd owners.

Endpoint support remains small:

```text
endpoint_support_size_hist={1:2204, 2:9432, 3:34}
```

## Tight Canonical Row

For the tight row

```text
S={1,2,3,4,5,6,7,8,9,10,11,13,84}
```

the extractor gives

```text
even_safe=107/245 in 26 components
bad_core=314/735 in 26 components
survivor=1/105 in 4 components
E components: full_bad=22, mixed=4, clean=0
min_survivor_gap=1/588
max_cover_total=3
cover_pair_hist={(1,1):14, (1,2):6, (2,1):6}
```

The largest displayed bad-core components have length `1/49` and are bounded
by the moving even wall `E:84`; examples include

```text
[43/588,55/588]   B0=(1,)     B1=(11,13)
[71/588,83/588]   B0=(1,)     B1=(7,9)
[169/588,181/588] B0=(7,13)   B1=(3,)
[211/588,223/588] B0=(5,11)   B1=(3,)
```

This is compatible with HYP-3431's corridor-fence proof: most even-safe
components are swallowed locally, but four mixed components leave the positive
relocation windows.

## Proof Target

The next lemma should not be stated as a scalar lower bound on branch measure.
It should be stated as a gluing impossibility:

```text
No primitive covering row can glue the local minimal B0/B1 odd-owner covers
across all E_safe components without leaving a survivor gap or emitting a
named sidecar/debt.
```

The legal exits are the existing ones: HYP-3431 corridor-fence reproduction,
HYP-3432 endpoint-budget ordering, HYP-3429 endpoint-spine compression,
HYP-3428 two-adic loss ledger, HYP-3427 wall words, HYP-3426 one-branch
mirror, HYP-3424 transfer, HYP-3423 topology-to-magnitude guardrail,
HYP-3417/HYP-3420 owner-current routing, HYP-3437 overlap-tax Menger cuts, and
HYP-3129/HYP-3421 signed-SPEC/Rprime.

After rebase, HYP-3437 is the graph-cut sibling of this packet.  HYP-3436
supplies the local bad-core atoms and minimal owner-cover signatures; HYP-3437
should use them as the incidence layer for overlap-tax/no-gluing certificates.

## Tournament Analysis

Vertices are proof obligations / exact interval-cover ledgers, not runners.

```text
pairwise_observable =
  predicate retention + cover exactness + endpoint payload
  + two-adic compatibility + scalar-firewall safety

switch =
  higher proof-facing weighted score; ties use declared code order

score_hist={-3:1, 66:1, 77:1, 78:1, 84:1, 91:1, 95:1, 97:1, 99:1}
directed_3cycles=0
hamiltonian_path =
  B00_local_bad_core_component_ledger
  -> B01_minimal_odd_owner_subcovers
  -> B02_endpoint_gate_wall_labels
  -> B03_survivor_gap_tax_certificate
  -> B04_overlap_tax_bridge
  -> B05_owner_current_exception_router
  -> B06_two_adic_descent_loss_ledger
  -> B07_topology_magnitude_legality_guard
  -> B08_scalar_bad_measure_shortcut
```

Assumption challenge: runners, gaps, fixed sections, section boundaries,
wall-crossing events, residues, endpoint owners, branch bad intervals,
even-safe components, bad-core components, survivor gaps, and proof obligations
were considered.  The chosen quotient preserves whether `E_safe` is swallowed
by `B0_odd cap B1_odd`; scalarization destroys odd-owner covers, even endpoint
gates, and survivor-gap placement.
