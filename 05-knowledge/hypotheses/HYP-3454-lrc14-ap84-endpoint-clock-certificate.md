---
id: HYP-3454
title: LRC14 AP84 endpoint-clock certificate
status: EVIDENCE / endpoint-inequality and residue-clock certificate; not an LRC14 proof
source: codex-2026-06-29 continuation of HYP-3431, HYP-3433, HYP-3439, HYP-3452, and HYP-3453
tangent: T1414
technique: LTI-414
tournament_technique: LTT-314
script: 04-computation/lrc14_ap84_endpoint_clock_certificate_codex_20260629.py
result: 05-knowledge/results/lrc14_ap84_endpoint_clock_certificate_codex_20260629.out
reflection: 07-reflections/lrc14-ap84-endpoint-clock-certificate-codex-20260629.md
related:
  - HYP-3457
  - HYP-3456
  - HYP-3453
  - HYP-3452
  - HYP-3451
  - HYP-3450
  - HYP-3439
  - HYP-3438
  - HYP-3437
  - HYP-3436
  - HYP-3435
  - HYP-3434
  - HYP-3433
  - HYP-3431
  - HYP-3429
  - HYP-3427
  - HYP-3426
  - HYP-3425
  - HYP-3422
  - HYP-3417
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3454: LRC14 AP84 Endpoint-Clock Certificate

## Claim

HYP-3454 narrows the HYP-3452 AP-tail phase theorem to two proof-facing
obligations:

```text
finite transients: m=1..4
rank-one endpoint interval: m>=5
residue boundary count: period-35 escape clock
```

For the canonical family

```text
S_m = {1,2,...,11,13,84m}
```

the rank-one component-cover interval for `m>=5` is exactly

```text
I_m=[(14ceil(48m/7)+1)/(588m),(14ceil(48m/7)+13)/(588m)].
```

This is the HYP-3433 endpoint-spine interval, now checked as the HYP-3450
component-cover escape and priced against the HYP-3431 low corridor.
Incoming HYP-3453 supplies the bank-level gate-escape transversal companion:
when a dead-cover obstruction exists, it should route through a low-rank
survivor gate.  HYP-3454 is the AP-tail endpoint-clock clause for that route.
HYP-3456 now supplies the floor-count derivation of the residue boundary
clock that HYP-3454 left as a sampled sidecar, and HYP-3457 closes the finite
mixed transient packet `m=1..4`.

## Exact Readout

Script:

```text
04-computation/lrc14_ap84_endpoint_clock_certificate_codex_20260629.py
```

Stored output:

```text
05-knowledge/results/lrc14_ap84_endpoint_clock_certificate_codex_20260629.out
```

Readout:

```text
checked_m=1..70
theorem_tail=m>=5
endpoint_failures_on_checked_tail=[]
min_left_corridor_margin=1/41160 at m=70
min_right_corridor_margin=1/2940 at m=5
length_residual_range=[0,0]
```

Every checked `m>=5` has endpoint labels

```text
L[E:84m] R[E:84m]
```

and endpoint rank `1`.  The interval lies inside the fixed HYP-3431 low
corridor `[8/49,6/35]`, while its endpoints are exactly adjacent moving
`E:84m` walls.  Thus the endpoint part of the AP-tail bridge is no longer a
sampling phenomenon: it has a closed-form address and exact containment
margins.

The transient side is now the HYP-3457 finite packet:

```text
m=1: L[E:84]  R[B1:5], rank=2, best=[33/196,6/35]
m=2: L[E:168] R[B1:5], rank=2, best=[197/1176,6/35]
m=3: L[E:252] R[B1:5], rank=2, best=[295/1764,6/35]
m=4: L[E:336] R[B1:5], rank=2, best=[131/784,6/35]
```

## Mod-35 Clock

The checked escape-count law is:

```text
escapes(m)=2*(floor(12m/35)+d[(m-1) mod 35])
```

with correction vector

```text
[2,2,1,1,1,1,1,2,1,1,2,1,1,1,1,2,2,1,1,1,1,2,2,1,1,2,1,1,2,1,2,2,1,1,0].
```

The script checks:

```text
beatty_clock_failures=[]
checked_shift_failures=[]
```

The formula itself has the symbolic shift

```text
floor(12*(m+35)/35)=floor(12m/35)+12,
```

so predicted escapes shift by `24`.  HYP-3456 derives this residue clock as
the actual high-gap floor count over the HYP-3431 fixed low corridors.

## Proof Pull

HYP-3439 reduced the AP/84m rescue-core bridge to the canonical rank-`6`
base case plus rank-`5` AP-tail descent.  HYP-3452 found the phase split.
HYP-3454 makes the descent sharper:

```text
1. Use HYP-3457 for the four mixed transients m=1..4.
2. Prove the closed-form endpoint interval I_m for every m>=5 by low-corridor
   containment and moving E:84m gap inequalities.
3. Use HYP-3456's period-35 floor count for the remaining component
   boundaries.
```

This is still not a proof of LRC14.  It is a smaller bridge lemma for the
current AP-tail/rescue-core corridor.  Raw dead fraction remains a scalar
warning only: HYP-3452 already showed that its checked maximum is a residue
artifact at `m=35`, while the proof clock starts at the endpoint phase `m=5`.
HYP-3462 now closes the AP-tail bookkeeping by importing the HYP-3431
branch-union carrier and splicing HYP-3454/HYP-3456/HYP-3457 into HYP-3439.
The active follow-up is non-AP transfer through HYP-3453/HYP-3451/HYP-3455.

## Tournament Analysis

Vertices are AP-tail proof obligations, not runners or raw arcs.

```text
pairwise_observable =
  endpoint exactness + low-corridor containment + residue-clock payload
  + scalar-firewall safety
score_hist={15:1,45:1,50:2,55:1,56:1,58:1}
directed_3cycles=0
hamiltonian_path=
  closed_form_endpoint_interval
  -> low_corridor_containment_inequalities
  -> moving_E84m_gap_certificate
  -> finite_transients_m1_to_m4
  -> mod35_escape_boundary_clock
  -> component_cover_reaudit
  -> raw_dead_fraction_peak
```

Assumption challenge: runners, `m`-values, residues mod `35`, fixed low
corridor sections, moving `E:84m` gaps, endpoint walls, component-cover graph
nodes, branch blockers, and proof obligations were considered.  The chosen
quotient preserves the AP-tail branch-union escape predicate, exact endpoint
labels, low-corridor containment margins, and the residue-count sidecar.  It
destroys non-AP wall geometry and full component adjacency, so it is only the
AP-tail bridge lemma needed by HYP-3439/HYP-3452.
