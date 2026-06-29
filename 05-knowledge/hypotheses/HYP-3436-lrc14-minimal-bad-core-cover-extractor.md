---
id: HYP-3436
title: LRC14 minimal bad-core cover extractor
status: EVIDENCE / exact local bad-core cover classification; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3435 and the HYP-3422/HYP-3425 two-adic bad-core chain, using HYP-3434 overlap tax and HYP-3429 endpoint-spine guardrails
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

HYP-3436 inverts HYP-3435's positive two-adic branch-cover certificate.

Instead of only reporting a surviving branch witness in

```text
E_safe(1/14) cap (branch0_good union branch1_good),
```

this packet will classify the complementary two-color obstruction

```text
E_safe(1/14) cap B0_odd cap B1_odd.
```

For each `E_safe` component, the executable extractor records exact bad-core
components, the minimal branch-0 and branch-1 odd-owner subsets that cover each
bad component, endpoint owners, survivor gaps, and a Tournament Analysis over
proof obligations rather than runners.  The proof-facing hope is now sharper:
a full counterexample would have to glue local minimal covers into a global
bad-core cover, while HYP-3426/HYP-3427/HYP-3429/HYP-3434 already show that the
needed overlap tax and endpoint sidecars are tightly constrained.

## Executable Readout

Script:

```text
04-computation/lrc14_minimal_bad_core_cover_extractor_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_minimal_bad_core_cover_extractor_codex_20260629.out
```

The script reuses the HYP-3435 row bank: `15` structured rows plus `120`
deterministic random primitive covering rows.  The exact identity

```text
E_safe minus (B0_odd cap B1_odd)
= E_safe cap (branch0_good union branch1_good)
```

holds by measure in every row.

```text
rows_audited=135
branch_union_measure_identity=135/135
positive_survivor_rows=135/135
total_even_safe_components=17164
total_fully_bad_safe_components=3770
total_positive_safe_components=13394
total_bad_core_components=11670
total_survivor_components=15868
max_minimal_two_color_owner_count=6
```

The local obstruction is much smaller than the raw component count suggests.
The minimal cover-size histogram is:

```text
(1,1): 10288
(1,2) and (2,1): 531 each
(2,2): 184
all other signatures: 136 total
```

Thus `10288/11670` bad-core intervals are already covered by one odd owner in
branch `0` and one odd owner in branch `1`.  Only eight checked intervals need
five owners on one branch, and none need total owner count above `6`.

The endpoint side is even more rigid:

```text
endpoint_support_hist={(1,1):11634, (1,2):18, (2,1):18}
```

So almost every bad-core interval has a single labelled wall at each endpoint.
This is exactly the kind of sidecar HYP-3430/HYP-3432/HYP-3434 say scalar
tails are not allowed to forget.

## Tight Row

For the tight canonical row

```text
S=(1,2,3,4,5,6,7,8,9,10,11,13,84)
```

the extractor recovers the HYP-3425 identity:

```text
E_safe=107/245 in 26 components
bad_core=314/735 in 26 components
survivor=1/105 in 4 components
fully_bad_safe_components=22
positive_safe_components=4
min_survivor_cell=1/588
max_core_cell=1/49
```

The four largest core cells have covers such as

```text
[43/588,55/588]  B0=(1,)      B1=(11,13)  endpoints E:84,E:84
[19/196,23/196]  B0=(1,)      B1=(9,)     endpoints E:84,E:84
[71/588,83/588]  B0=(1,)      B1=(7,9)    endpoints E:84,E:84
[169/588,181/588] B0=(7,13)   B1=(3,)     endpoints E:84,E:84
```

This tight row is no longer just a small positive measure.  It is a finite
ledger of local two-color covers separated by four surviving gaps.

## Canonical Tail Clue

The script also scans the canonical family

```text
S_m={1,2,3,4,5,6,7,8,9,10,11,13,84m}
```

through `m=30`.  There are no failures, and every bad-core component is a
singleton/singleton branch-owner cover from `m=3` onward:

```text
first_all_core_components_singleton_singleton=3
canonical_failures=[]
```

This looks like the local-cover version of HYP-3431's corridor-fence proof:
after the early `m=1,2` overlaps, the high even wall `E:84m` cuts the bad core
into tiny cells, each destroyed locally by one low odd owner on each branch,
while a controlled set of gaps survives.

## Candidate Lemma

For every primitive covering row `S=O union 2E` with `|S|=13`, the local
minimal covers of

```text
E_safe cap B0_odd cap B1_odd
```

cannot concatenate into a cover of all `E_safe` components unless one of the
following named debts appears:

```text
owner-current exception
endpoint-spine failure
two-adic descent debt
overlap-tax deficit
state-lift / exact-period residual
```

Equivalently, HYP-3422 should be proved by a finite labelled packet theorem
over minimal branch-owner covers, rather than by a global scalar branch-measure
lower bound.

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

Vertices are proof obligations and interval-cover certificate interfaces, not
runners or raw arcs.  Alternate vertex sets considered were safe components,
bad-core intervals, endpoint events, odd-owner cover tokens, Fourier modes, and
runners.  The selected quotient preserves the exact survivor predicate

```text
E_safe minus B0_odd cap B1_odd
```

and destroys only raw runner ordering plus open-endpoint conventions.

Readout:

```text
tie_hamiltonian_path =
B00_minimal_bad_core_cover_extractor
-> B01_two_color_failure_normal_form
-> B02_endpoint_owner_set_cover_lemma
-> B03_overlap_tax_discharge_bridge
-> B04_component_spine_router
-> B05_two_adic_descent_induction
-> B06_owner_current_exception_exit
-> B07_harmonic_scalar_sidecar
-> B08_raw_topology_or_named_constant

directed_3cycles=0
edge_flips_against_scalar_first=36
```

The challenge to the default assumption is explicit: tournament vertices need
not be runners.  Here the most useful tournament is over proof obligations and
cover certificates because that is what the quotient must preserve to prove
the LRC14 branch lemma.
