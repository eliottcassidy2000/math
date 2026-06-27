---
id: HYP-3017
title: LRC14 automatic sidecar route-purity audit
status: EVIDENCE / packet-field audit and quotient guardrail; not a proof
source: codex-2026-06-26-S181
script: 04-computation/lrc14_labelled_packet_counterexample_classifier_codex_20260624.py
result: 05-knowledge/results/lrc14_labelled_packet_counterexample_classifier_codex_20260624.out
related:
  - HYP-2963
  - HYP-3016
  - HYP-3015
  - HYP-3014
  - HYP-3013
  - HYP-3012
  - HYP-3011
  - HYP-3010
  - HYP-3009
  - HYP-3008
  - HYP-2999
  - HYP-2998
  - THM-572
  - OPEN-Q-108
---

# HYP-3017: LRC14 Automatic Sidecar Route-Purity Audit

## Claim

The HYP-2963 labelled-packet classifier now emits the automatic and
power-lift sidecar fields proposed by HYP-3009 through HYP-3013 and stressed by
the HYP-3016 automaton fiber-mixing guardrail:

```text
automatic_word over M/F/C
automatic_counts
moser_double_breaks
fibbinary_double_breaks
lacunary_tail_ratio
lacunary_max_ratio
q_factorization
unit_excess_apex
power_lift_guard
automatic_filter_exit
```

These fields are useful proof telemetry, but they are not safe quotient
coordinates by themselves.

## Computation

The regenerated HYP-2963 run keeps the original labelled-packet census:

```text
audited rows              21913
below threshold           0
tight rows                2: AP, GW 12->24
M<=2/27 low packets       7
unknown packets           0

Q-WITNESS                 14676
BOUNDARY-AP-GW            2
BOUNDARY-PETAL-SPORADIC   4
K33-STATE-LIFT            3
COVERING-MOMENT           7228
```

The new automatic sidecar audit reports:

```text
automatic word fibers                 225
mixed-route word fibers               143
mixed-family word fibers              178
unit-excess apex rows                 169
rows with perfect-power payload guards 9606
rows with Moser n->2n phase breaks    21913
rows with fibbinary n->2n breaks      0
```

The main negative result is the AP/GW word:

```text
AP:        MFCMMCCFFFCCC
GW 12->24: MFCMMCCFFFCCC
```

That same word fiber has `639` rows:

```text
Q-WITNESS                 603
BOUNDARY-AP-GW              2
BOUNDARY-PETAL-SPORADIC     1
COVERING-MOMENT            33
```

So the Moser/fibbinary/carry word simultaneously contains direct q-witness
rows, both tight boundary atoms, a petal discharge row, and covering/open-Haar
rows.  It cannot replace exact `M`, q-threshold, endpoint geometry, C27/K33
labels, or covering labels.

## Fermat-Catalan / Power-Lift Readout

The power-lift guard is deliberately conservative.  In the low frontier:

```text
P10+GW        M=2/27, q=27=3^3
petal 10->20 M=2/27, q=27=3^3
petal 13->26 M=2/27, q=27=3^3
P10+K33       M=2/27, q=27=3^3
```

Further positive-open rows expose payload powers such as:

```text
two drop(10,13)->add(20,26): p+q=25=5^2
single swap 12->48:          p=4=2^2
single swap 13->104:         p=8=2^3
```

This supports HYP-3009's rule: Fermat-Catalan style information is a
no-lift/finite-exception ledger after packet labels are attached, not a scalar
proof certificate.

## Assumption Challenge

The quotient vertices considered were runners, gaps, fixed circle sections,
section boundaries, residues, endpoint owners, wall-crossing events, cover
arcs, Fourier modes, automatic-language states, valuation guards, and proof
obligations.

The chosen sidecar vertex set is automatic/power packet fields attached to
HYP-2963 proof-route packets.  This preserves the LRC predicate only when exact
`M`, q-threshold, strict safe components, and route labels remain attached.  It
destroys endpoint-owner geometry and theorem route if used as a standalone
automatic word.

## Proof Target

Every future sequence-shadow or finite-automaton quotient in the LRC14 stack
should first run a route-purity test:

```text
fiber(label) has one theorem route
or label remains a sidecar
or mixed fiber is split by exact M/q, endpoint, C27/K33, Haar/Fejer,
   covering, valuation, or THM-572 residual labels.
```

For now, automatic words fail this test on the HYP-2963 bank.  The next useful
step is not another scalar automaton search; it is family templates for the
large mixed fibers, starting with `MFCMMCCFFFCCC`, compared against HYP-3015
lonely-profile barcode sidecars and HYP-3016 magnitude-cocycle fibers.
