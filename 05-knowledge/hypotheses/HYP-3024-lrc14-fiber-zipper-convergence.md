---
id: HYP-3024
title: LRC14 zipper-fiber convergence, Erdos-Turan clocks, and Henselian unit rule
status: EVIDENCE / full-bank convergence scout and proof-interface target; not a proof
source: codex-2026-06-26-S188
tangent: T1105
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
  - HYP-3014
  - HYP-3009
  - HYP-3008
  - HYP-2997
  - HYP-2995
  - HYP-2991
  - HYP-2989
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3024: LRC14 Zipper-Fiber Convergence

## Claim

The full HYP-2963 packet bank supports a sharper version of the
HYP-3020/HYP-3023 bridge:

```text
residue-terminal zipper
  + coarse Erdos-Turan discrepancy clocks
  + Henselian unit-root rule
```

is already enough to preserve the direct boundary/open LRC predicate on the
tested full bank, even though it is not enough to classify every theorem
route.  Exact Erdos-Turan clocks at `14,27,41` are too sharp on this bank and
degenerate to packet-level addressing; the useful target is the coarsened
status-preserving gate.

This is not a proof of LRC14.  It is a proof-interface result: it identifies
which clocks matter for the question "does this orbit hit the safe box even
once?" and which clocks are only route schedulers.

## Computation

Script:

```text
04-computation/lrc14_fiber_zipper_convergence_codex_s188.py
```

Stored output:

```text
05-knowledge/results/lrc14_fiber_zipper_convergence_codex_s188.out
```

The script audits the default HYP-2963 `21913`-packet bank with zipper gates:

```text
automatic_word
residue_terminal_fiber
residue_plus_et
residue_plus_unit_hensel
et_unit_convergence
coarse_et_unit_gate
magnitude_cocycle
magnitude_et_unit
barcode_packet_zipper
```

The full-bank convergence table:

```text
split                      fibers    mixR    mixF    mixS    maxF  maxMix  purity
automatic_word                225     143     178       1    1179    1179   36.4%
residue_terminal_fiber      16555     265     265       2      30      30   98.3%
residue_plus_et             21913       0       0       0       1       0  100.0%
residue_plus_unit_hensel    19301     220     221       2      22      22   98.8%
et_unit_convergence         21913       0       0       0       1       0  100.0%
coarse_et_unit_gate         21702      15      16       0       4       4   99.9%
magnitude_cocycle           21909       0       0       0       2       0  100.0%
magnitude_et_unit           21913       0       0       0       1       0  100.0%
barcode_packet_zipper       21913       0       0       0       1       0  100.0%
```

The key new data point is:

```text
coarse_et_unit_gate:
  fibers=21702
  mixed_route=15
  mixed_status=0
  max_mixed=4
```

So the coarse gate is more compressed than the exact magnitude cocycle on this
bank, preserves boundary/open status, and leaves only open-route mixing to be
routed by proof families.

## Henselian Unit Rule

HYP-3020 counted roots and singular roots of `A_S(x)=sum_{v in S} x^v` modulo
`2,3,7`.  HYP-3024 refines this: roots in `F_p^*` are genuine Hensel unit
clocks, while the forced root at `x=0` is tracked separately as scale or
nilpotent debt.

Full-bank readout:

```text
coarse_et_clock_fibers=1761 largest=242
hensel_unit_count_fibers=29 largest=7660
hensel_unit_exact_fibers=73 largest=7660
packets_with_singular_unit_root=3358
packets_with_zero_singular_debt=3203
```

Thus p-adic data is not a scalar discriminator.  It is a routing rule:
unit-root singularities feed local lift debt; zero-singular roots feed scale
debt; neither should be confused with the other.

## Target Word

For the AP/Goddyn-Wong word

```text
MFCMMCCFFFCCC
```

the script finds:

```text
rows=639
routes={'Q-WITNESS': 603, 'BOUNDARY-AP-GW': 2,
        'BOUNDARY-PETAL-SPORADIC': 1, 'COVERING-MOMENT': 33}
distinct_M=33
distinct_ET_clocks=639
distinct_unit_rules=28
```

Inside the word:

```text
residue_terminal_fiber:    fibers=181 mixed_route=27
residue_plus_et:           fibers=639 mixed_route=0
residue_plus_unit_hensel:  fibers=294 mixed_route=28
et_unit_convergence:       fibers=639 mixed_route=0
coarse_et_unit_gate:       fibers=585 mixed_route=9
magnitude_cocycle:         fibers=638 mixed_route=0
magnitude_et_unit:         fibers=639 mixed_route=0
```

Exact ET clocks fully split this finite target word, but that is too close to
an address coordinate.  The coarse ET+unit gate leaves only `9` mixed route
fibers in the target word, all open, with largest mixed fiber size `4`.

## Tournament Analysis

Vertices were quotient/proof-carrier bundles, not runners.  The observable was
route purity, max mixed-fiber size, Erdos-Turan discrepancy, Henselian unit
stability, magnitude retention, topology retention, packet-label retention,
finite-state checkability, and proof cost.  The switch was majority comparison
of those observable vectors, with the zipper list as the tie Hamiltonian path.

Fingerprint:

```text
score_hist={0: 1, 1: 1, 2: 1, 3: 1, 5: 2, 6: 2, 8: 1}
directed_3cycles=2
scc_sizes=[4, 1, 1, 1, 1, 1]
hamiltonian_path_count=5
score_order=magnitude_et_unit > barcode_packet_zipper > magnitude_cocycle >
            et_unit_convergence > residue_plus_et > coarse_et_unit_gate >
            residue_plus_unit_hensel > residue_terminal_fiber >
            automatic_word
```

The size-4 SCC is the warning: exact magnitude, exact ET, Hensel-unit, and
packet/barcode labels are interacting proof clocks, not a single scalar
ranking.

## Assumption Challenge

Alternate vertex sets considered: runners, gaps, fixed circle sections,
section boundaries, residues, Fourier clocks, p-adic unit roots, zero-root
scale debt, magnitude cocycles, barcode fibers, packet labels, and proof
obligations.

The chosen vertices are quotient/proof-carrier bundles because they preserve
the tested LRC predicate: theorem-route purity and boundary/open status at
threshold `1/14`.  They destroy raw runner identity, exact endpoint owners,
and full Fejer atom banks unless later zipper fields are attached.

The challenged assumption is that Hensel data should be read as all residue
roots equally.  For proof routing, p-adic units and the forced zero root play
different roles.

## Proof Target

Coarse zipper-fiber convergence theorem:

```text
Inside each automatic/residue-terminal fiber, the coarse ET+unit gate either
  (a) is boundary/open fiber-constant,
  (b) emits a bounded Erdos-Turan discrepancy certificate,
  (c) emits Henselian unit-lift debt or zero-root scale debt,
  (d) descends to a familywise magnitude formula, or
  (e) routes to named K33/F7/THM-572 residual debt.
```

This target is weaker and cleaner than immediate route purity.  It focuses on
the direct LRC question first: whether the orbit has any strict safe interval.
