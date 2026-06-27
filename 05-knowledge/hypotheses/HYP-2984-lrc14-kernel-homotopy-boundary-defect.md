---
id: HYP-2984
title: LRC14 kernel homotopy and smoothing-switchboard boundary-defect ledger
status: PROOF-INTERFACE / exact packet audit, finite route matrix, and admissible-kernel guardrail; not a proof
source: codex-2026-06-24-S164
artifacts:
  - 04-computation/lrc14_kernel_homotopy_boundary_defect_codex_s164.py
  - 05-knowledge/results/lrc14_kernel_homotopy_boundary_defect_codex_s164.out
  - 07-reflections/lrc14-kernel-homotopy-boundary-defect-codex-s164.md
  - 04-computation/lrc14_smoothing_switchboard_codex_s164.py
  - 05-knowledge/results/lrc14_smoothing_switchboard_codex_s164.out
  - 07-reflections/lrc14-smoothing-switchboard-codex-s164.md
related:
  - HYP-2983
  - HYP-2982
  - HYP-2981
  - HYP-2980
  - HYP-2979
  - HYP-2978
  - HYP-2977
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2969
  - HYP-2968
  - HYP-2966
  - HYP-2965
  - HYP-2963
  - HYP-2949
  - HYP-2948
  - HYP-2908
  - THM-548
  - THM-572
  - OPEN-Q-108
---

# HYP-2984: LRC14 Kernel Homotopy And Smoothing Switchboard

This hypothesis connects HYP-2981's packet Fejer interval certificates with
HYP-2982/HYP-2983's analytic smoothing and Kaczynski boundary language.

Working thesis:

```text
Changing a smoothing or Fourier kernel is theorem-safe only when it either
preserves the labelled packet certificate or emits a named boundary-defect atom.
```

Equivalently:

```text
choose the smoothing kernel only after the labelled packet route is known.
```

This is a proof-interface guardrail, not a proof of LRC14.  It says where a
future smoothing proof is allowed to change kernels and what the proof is still
required to remember when a quotient forgets geometry.

## Kernel Homotopy Audit

Script:

```text
04-computation/lrc14_kernel_homotopy_boundary_defect_codex_s164.py
```

Stored output:

```text
05-knowledge/results/lrc14_kernel_homotopy_boundary_defect_codex_s164.out
```

The script reuses the exact Haar/Baire safe-component machinery from S147 and
the taut-bridge audit from HYP-2975/S155.  For each selected row it records:

```text
exact M, qdiv, strict safe Haar mass, safe component count,
largest component width, kernel support radius eps < width/2,
taut endpoint count, route label.
```

Selected readout:

```text
AP, GW 12->24: safe_mu=0, eps<0, taut_vertices=6 each.
near/K33 12->36: safe_mu=1/1260, max_width=1/2520, eps<1/5040.
petal 10->20 and P10+GW: safe_mu=1/980, eps<1/3920.
petal 13->26: safe_mu=1/182, eps<1/728.
P10+K33: safe_mu=4/2205, eps<1/3920.
covering 12->84: safe_mu=563/105105, eps<3/3920.
covering 12->168: safe_mu=263/30030, eps<23/23520.
few-apex 6->28: safe_mu=7/858, eps<5/3696.
```

The AP and Goddyn-Wong boundary atoms have the same six taut endpoint
transfers:

```text
t=1/14:  1:g1  -> 13:g1
t=3/14:  5:g1  -> 9:g1
t=5/14:  3:g1  -> 11:g1
t=9/14:  11:g1 -> 3:g1
t=11/14: 9:g1  -> 5:g1
t=13/14: 13:g1 -> 1:g1
```

All six have endpoint-owner pair sum `0 mod 14`.  This is exactly the
boundary-defect payload a kernel quotient is not allowed to erase.

The one-swap sanity scan through replacement `<=160` checks `1911` primitive
rows:

```text
open-stable: 1910
zero-open: 1
zero-open row: 12->24
minimum positive safe mass: 1/1260 at 12->36
```

So, in this finite local audit, every non-GW one-swap row has a positive
regular-open component and hence a nonzero kernel-stability radius.

For a row with a strict safe component `(a,b)`, any symmetric kernel supported
inside radius `<(b-a)/2` can be centered at `(a+b)/2` without crossing a danger
boundary.  This gives a theorem-facing homotopy invariant:

```text
open component + packet label + support bound
  -> admissible local kernel deformation.
```

For a zero-open row, this support radius is zero.  A smoothing kernel of any
positive support smears across danger boundaries, so the proof must instead
carry a boundary-defect atom:

```text
taut endpoint + owner transfer + pair-sum/mod-shell label
  -> named defect emitted by the deformation.
```

HYP-2981's Fejer records can use this interface directly.  The Fejer center is
already chosen inside a strict safe component for positive rows; for AP/GW
there is no such component, so they remain PSD-blind equality atoms rather than
failed certificates.

Rebase signal from the concurrent S164 coordination sync: the shared
coordination note describes the finite defect atlas as AP/GW/K33.  This file
keeps the distinction typed.  AP/GW are zero-open boundary-defect atoms in the
exact audit; K33 is a named residual/state-lift route that a homotopy may emit
after leaving the AP/GW taut stratum, not a scalar-equivalent defect.

## Smoothing Switchboard Audit

Script:

```text
04-computation/lrc14_smoothing_switchboard_codex_s164.py
```

Stored output:

```text
05-knowledge/results/lrc14_smoothing_switchboard_codex_s164.out
```

The script imports:

- HYP-2981/S162 packet-anchored Fejer interval scaffold;
- HYP-2979 Ramanujan exact-period primitive projector;
- HYP-2982/HYP-2983 analytic-sieve guardrail language.

For `16` named/probe packets it records:

```text
q_threshold,
exact safe measure,
safe-component count,
q=14 weak/strict/boundary primitive phase count,
first weak and strict Ramanujan q through q<=42,
selected Fejer interval certificate status,
chosen smoothing/proof route,
retained side channel,
forbidden scalarization.
```

Route counts:

```text
AP_GW_BOUNDARY_ATOM                 2
COVERING_LIFT_OR_BOUNDARY_MOMENT    2
INTERVAL_FEJER_CERTIFICATE          5
K33_STATE_LIFT                      1
RAMANUJAN_PRE_SPLIT_THEN_FEJER      6
```

The five selected positive hard rows from HYP-2981 remain interval-certified:

```text
P10+GW                              degree 280
covering 12->168                   degree 63
near/K33 12->36                    degree 159
single swap 6->63                  degree 266
two drop(12,13)->add(14,29)        degree 41
```

AP and Goddyn-Wong are the only audited `safe_mu=0` boundary atoms.  Covering
rows without direct Fejer records route to lift / boundary-moment certificates.
`P10+K33` routes to K33 state-lift debt.  Petal, q14-front, q-witness, and
repeated-prime probe rows route through Ramanujan pre-split and Fejer/twist
handoff.

## Switch Contracts

The executable contracts are the new content:

```text
AP/GW boundary:
  retain endpoint zero-credit pairs plus Kaczynski approach class.
  forbidden: open-measure averaging or raw Fejer negativity.

Interval Fejer:
  retain packet key, rational center, Fejer degree, interval upper bound.
  forbidden: floating Fejer value without interval/packet label.

Ramanujan pre-split:
  retain first strict q, primitive phase packet, endpoint labels.
  forbidden: qdiv-only or scalar Ramanujan profile.

Covering/lift:
  retain lift chart, endpoint owners, exact safe mass, late-q label.
  forbidden: squarefree mu^2/phi collapse of repeated-prime q.

K33 state-lift:
  retain K33 owner packet plus HYP-2908/THM-572 debt.
  forbidden: large-sieve scalar tail without state-lift side channel.
```

Thus the admissible-smoothing lemma should not be phrased as "choose Fejer" or
"use the large sieve."  It should be:

```text
given a labelled packet, choose the least-forgetting kernel that preserves the
next proof predicate, and record the residual labels that kernel cannot see.
```

## Tournament Analysis

Vertices are proof carriers in the kernel-deformation and route-selection
steps, not runners.  The kernel audit uses:

```text
open_component_certificate
boundary_defect_atom
packet_fejer_certificate
kaczynski_approach_class
analytic_smoothing_kernel
ramanujan_exact_period_packet
raw_kernel_scalar
```

The switchboard audit uses:

```text
labelled_smoothing_switchboard
interval_fejer_certificate
boundary_moment_lift_chart
kaczynski_boundary_defect
ramanujan_exact_period_presplit
large_sieve_middle_bound
raw_scalar_smoothing_choice
```

Both audited carrier tournaments are transitive:

```text
packet_fejer_certificate >
open_component_certificate >
boundary_defect_atom >
kaczynski_approach_class >
ramanujan_exact_period_packet >
analytic_smoothing_kernel >
raw_kernel_scalar

labelled_smoothing_switchboard >
interval_fejer_certificate >
boundary_moment_lift_chart >
kaczynski_boundary_defect >
ramanujan_exact_period_presplit >
large_sieve_middle_bound >
raw_scalar_smoothing_choice
```

Assumption challenge: the proof vertices need not be runners or arcs.  This
lane explicitly considers kernels, safe components, boundary events,
endpoint-owner pairs, packet fibers, smoothing defects, Ramanujan modes, lift
charts, interval certificates, route handoffs, and proof obligations.  The
quotient preserves the LRC predicate only by retaining the selected packet
route and the side channel declared by the audit.  It destroys raw runner
ownership and scalar smoothing margins unless those labels are retained or
formally reconstructed.

## Theorem Target

Candidate lemma:

```text
Every primitive LRC14 source packet admits a labelled smoothing route:
AP/GW boundary, interval Fejer certificate, Ramanujan exact-period handoff,
covering/lift boundary moment, K33/state lift, or a named Kaczynski defect.
No route may forget endpoint owners, exact-period labels, safe-open/boundary
state, Fejer interval fields, or state-lift debt unless another retained route
formally reconstructs them.
```

The familywise missing work is to scale the support-radius and switchboard
readout from the named audit bank to HYP-2963 packet families, then prove the
emitted defects are exhausted by AP/GW and existing state-lift / interval-Fejer
routes.
