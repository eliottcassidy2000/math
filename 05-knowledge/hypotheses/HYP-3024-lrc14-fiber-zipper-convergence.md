---
id: HYP-3024
title: LRC14 fiber-zipper convergence via Erdos-Turan and Henselian unit rule
status: EVIDENCE / target-fiber convergence scout; not a proof
source: codex-2026-06-26-S188
tangent: T1104
script: 04-computation/lrc14_fiber_zipper_convergence_codex_s188.py
result: 05-knowledge/results/lrc14_fiber_zipper_convergence_codex_s188.out
related:
  - HYP-3023
  - HYP-3022
  - HYP-3021
  - HYP-3020
  - HYP-3019
  - HYP-3018
  - HYP-3017
  - HYP-3016
  - HYP-3015
  - HYP-3009
  - HYP-3008
  - HYP-2963
  - THM-572
  - LTI-171
  - LTI-170
  - LTI-167
  - OPEN-Q-108
---

# HYP-3024: LRC14 Fiber-Zipper Convergence

## Claim

The first HYP-3023 automatic-fiber target

```text
MFCMMCCFFFCCC
```

has a smaller proof carrier than exact magnitude:

```text
residue-terminal fiber
  -> Erdos-Turan residue discrepancy bins
  -> Henselian unit rule
```

On the target fiber, the Erdos-Turan bins explain most of the visible
contraction, and the Henselian unit rule already gives route-pure fibers.

This is not a proof of LRC14.  It is a convergence scout for the first zipper
fiber.  The same rule still needs a full-bank stress test before replacing the
exact HYP-3023 magnitude cocycle.

## Computation

Script:

```text
04-computation/lrc14_fiber_zipper_convergence_codex_s188.py
```

Stored output:

```text
05-knowledge/results/lrc14_fiber_zipper_convergence_codex_s188.out
```

The default run filters the HYP-2963 bank to the target automatic word before
exact packet computation:

```text
candidate_rows=639 of 21913
packets=639
```

The convergence ladder is:

```text
automatic_word                 fibers=1   mixed_route=1   max_mixed=639
residue_terminal_fiber         fibers=181 mixed_route=27  max_mixed=30
erdos_turan_residue_zipper     fibers=595 mixed_route=6   max_mixed=4
henselian_unit_zipper          fibers=424 mixed_route=0   max_mixed=0
et_hensel_unit_zipper          fibers=627 mixed_route=0   max_mixed=0
et_hensel_qrule_zipper         fibers=627 mixed_route=0   max_mixed=0
magnitude_cocycle              fibers=638 mixed_route=0   max_mixed=0
barcode_shadow                 fibers=639 mixed_route=0   max_mixed=0
packet_zipper                  fibers=639 mixed_route=0   max_mixed=0
```

Readout:

```text
Residue-terminal fibers are close but still mixed.
Erdos-Turan bins reduce mixed-route fibers 27 -> 6 and max mixed size 30 -> 4.
Henselian unit-rule bins reduce mixed-route fibers 27 -> 0.
Exact magnitude is no stronger for route purity on this target fiber.
```

## Erdos-Turan Tooth

The Erdos-Turan tooth records, for moduli `14,27,41`,

```text
exact L1 residue-discrepancy bucket
low-frequency sum_{h<=6} |sum_{v in S} exp(2*pi*i*h*v/q)| / h bucket
```

This is intentionally not exact magnitude.  Its role is a convergence gauge:
inside `MFCMMCCFFFCCC`, it contracts the largest residue-terminal mixed fiber
from `30` rows to `4`, leaving only six mixed-route fibers.

The largest remaining Erdos-Turan-mixed example is a q-witness/covering split:

```text
single swap 13->143  M=11/144 route=COVERING-MOMENT
single swap 13->31   M=1/13   route=Q-WITNESS
single swap 13->59   M=1/13   route=Q-WITNESS
single swap 13->115  M=1/13   route=Q-WITNESS
```

So ET sees the right pressure but still forgets a local unit/lift distinction.

## Henselian Unit Rule

For each `p in {2,3,7}`, define

```text
A_S(x)=sum_{v in S} x^v  over F_p^*
```

and retain:

```text
simple unit roots:    A_S(x)=0 and A'_S(x) != 0
singular unit roots:  A_S(x)=0 and A'_S(x) = 0
exponent counts on the cyclic unit group F_p^*
p-adic denominator-unit data for M=p0/q0 and Farey excess
```

The rule is Henselian in the narrow sense: a simple unit root is lift-stable,
while a singular unit root is local lift debt.  On the target automatic fiber,
this sidecar separates every theorem route:

```text
henselian_unit_zipper: fibers=424 mixed_route=0
```

This is the key new carrier.  It suggests that the target word's `33` exact
`M` values are overkill for route purity; a p-adic unit-lift signature may be
enough.

## Tournament Analysis

Vertices are zipper teeth / proof carriers, not runners.  Pairwise observable:

```text
route purity
max mixed-fiber size
convergence from the previous tooth
Erdos-Turan retention
Henselian-unit retention
magnitude retention
packet-label retention
proof cost
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
score_order=et_hensel_qrule_zipper > et_hensel_unit_zipper > packet_zipper >
            magnitude_cocycle > barcode_shadow > henselian_unit_zipper >
            erdos_turan_residue_zipper > residue_terminal_fiber >
            automatic_word
```

## Assumption Challenge

Alternate vertices considered: runners, gaps, fixed circle sections, residues,
Erdos-Turan Fourier modes, Henselian unit roots, denominator units, barcode
bars, packet labels, and proof obligations.

The chosen vertices are zipper teeth.  They preserve theorem-route purity and
boundary-versus-open status at threshold `1/14`.  They destroy raw runner
identity and raw automatic-word multiplicity, but only after retaining enough
analytic and p-adic data to explain the destroyed distinction.

## Next Pulls

1. Prove the Henselian unit split for the target word `MFCMMCCFFFCCC`.
2. Run the same ET/Hensel zipper on the full HYP-2963 bank; the script supports
   `--full-bank`, but the first full run was too slow for this pass.
3. If full-bank Henselian unit rule still has mixed fibers, classify the
   survivors by q-threshold, unit-excess lane, barcode support, and packet
   route.
4. Promote only the smallest route-pure compression into HYP-2963 packet
   fields; keep exact magnitude as fallback until the family lemma exists.
