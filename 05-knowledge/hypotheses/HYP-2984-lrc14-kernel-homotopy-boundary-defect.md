---
id: HYP-2984
title: LRC14 kernel homotopy and boundary-defect ledger
status: PROOF-INTERFACE / exact named-packet audit and smoothing guardrail; not a proof
source: codex-2026-06-24-S164
related:
  - HYP-2983
  - HYP-2982
  - HYP-2981
  - HYP-2980
  - HYP-2978
  - HYP-2977
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2949
  - HYP-2948
  - THM-572
  - OPEN-Q-108
artifacts:
  - 04-computation/lrc14_kernel_homotopy_boundary_defect_codex_s164.py
  - 05-knowledge/results/lrc14_kernel_homotopy_boundary_defect_codex_s164.out
  - 07-reflections/lrc14-kernel-homotopy-boundary-defect-codex-s164.md
---

# HYP-2984: LRC14 Kernel Homotopy And Boundary-Defect Ledger

This hypothesis connects HYP-2981's packet Fejer interval certificates with
HYP-2982/HYP-2983's analytic smoothing and Kaczynski boundary language.

Working thesis:

```text
Changing a smoothing or Fourier kernel is theorem-safe only when it either
preserves the labelled packet certificate or emits a named boundary-defect atom.
```

This is a proof-interface guardrail, not a proof of LRC14.  It says where a
future smoothing proof is allowed to change kernels.  Positive regular-open
safe components give a support radius for local kernel perturbations; zero-open
AP/Goddyn-Wong equality atoms have no such radius and must remain explicit
boundary defects.

## Computation

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

## Kernel Homotopy Rule

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

## Tournament Analysis

Vertices are proof carriers in the kernel-deformation step:

```text
open_component_certificate
boundary_defect_atom
packet_fejer_certificate
kaczynski_approach_class
analytic_smoothing_kernel
ramanujan_exact_period_packet
raw_kernel_scalar
```

Pairwise observable:

```text
retention of strict witness,
boundary equality,
packet labels,
kernel auditability,
formal interval path,
scalar-decoy resistance.
```

Switch/gauge: componentwise majority of the six retention scores, ties along
the declared Hamiltonian path.  Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
SCC_sizes=[1,1,1,1,1,1,1]
Hamiltonian path:
packet_fejer_certificate >
open_component_certificate >
boundary_defect_atom >
kaczynski_approach_class >
ramanujan_exact_period_packet >
analytic_smoothing_kernel >
raw_kernel_scalar
```

The ranking is intentional.  A raw smoothing kernel is below packet Fejer
certificates and open-component certificates because it can forget the LRC
predicate.  It becomes useful only after packet labels and boundary defects are
attached.

Assumption challenge: the proof vertices need not be runners or arcs.  This
lane explicitly considered kernels, safe components, boundary events,
endpoint-owner pairs, packet fibers, smoothing defects, and proof obligations.
The chosen quotient preserves the predicate "strict safe interval, AP/GW
boundary equality, or named residual debt"; it destroys raw runner ownership
only after that predicate is retained or reconstructed.

## Theorem Target

The useful theorem shape is:

```text
For every labelled LRC14 packet, a kernel deformation either preserves an
open-component certificate with a declared support radius, or it emits a
boundary-defect atom that routes to AP/GW equality, K33/state-lift debt,
Ramanujan/exact-period debt, or an existing Fejer interval certificate.
```

The missing work is familywise: lift the support-radius and boundary-defect
readout from named rows to HYP-2963 packet families, then prove the emitted
defects are exhausted by AP/GW and existing state-lift routes.
