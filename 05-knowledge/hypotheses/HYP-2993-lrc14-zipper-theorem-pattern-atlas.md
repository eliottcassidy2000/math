---
id: HYP-2993
title: LRC14 zipper theorem pattern atlas
status: SYNTHESIS / proof-interface pattern atlas; not a proof
source: codex-2026-06-24-S166
artifacts:
  - 04-computation/lrc14_zipper_theorem_pattern_atlas_codex_s166.py
  - 05-knowledge/results/lrc14_zipper_theorem_pattern_atlas_codex_s166.out
  - 07-reflections/lrc14-zipper-theorem-pattern-atlas-codex-s166.md
related:
  - HYP-2990
  - HYP-2991
  - HYP-2992
  - HYP-2989
  - HYP-2988
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2984
  - HYP-2981
  - HYP-2979
  - HYP-2978
  - HYP-2595
  - HYP-2594
  - HYP-2452
  - HYP-2450
  - HYP-2104
  - HYP-2101
  - THM-572
  - OPEN-Q-108
---

# HYP-2993: LRC14 Zipper Theorem Pattern Atlas

The recent LRC14 work is converging on a single proof interface:

```text
local certificate side A
local certificate side B
labelled interface where they meet
declared equality stops
named residuals when they do not meet
```

This pass concretizes HYP-2990's abstract no-free-slider zipper atlas inside
the LRC14/Haar/Fejer packet stack.  A quotient is admissible only
when every forgotten coordinate is reconstructed by the opposite tooth,
annihilated by orthogonality or a boundary stop, or emitted as a labelled
residual strong enough for the next theorem.

The schema is:

```text
For every primitive packet P in a classified family F:
  1. compute left teeth L(P) and right teeth R(P), with labels attached;
  2. if L(P) and R(P) meet compatibly, close the certificate;
  3. if P lands on a declared stop, record the equality/boundary atom;
  4. otherwise emit a residual whose labels are strong enough for the next theorem.
```

HYP-2993 is not the proof.  It is a theorem-design atlas meant to keep the
proof debt finite and prevent scalar quotients from hiding the last packets.

## Computed Atlas

Script:

```text
04-computation/lrc14_zipper_theorem_pattern_atlas_codex_s166.py
```

Stored output:

```text
05-knowledge/results/lrc14_zipper_theorem_pattern_atlas_codex_s166.out
```

The script records ten zipper patterns:

```text
haar_fourier_product
fejer_interval_packet
tope_cocircuit_wall
exposure_poset_kernel
ramanujan_exact_period
smoothing_kaczynski_policy
fixed_margin_johnson
apex_sheaf_gluing
convolution_irreducibility_lift
unit_distance_cyclotomic_norm
```

Each pattern stores:

```text
domain,
left teeth,
right teeth,
boundary stop,
residual atom,
quotient warning,
proof arrow,
related hypotheses,
retention vector.
```

The retention vector tracks:

```text
LRC predicate, local gluing, boundary stop, residual name,
formal checkability, cross-domain transfer, anti-scalar guardrail,
computable next test.
```

## Tournament Analysis

Vertices are zipper proof patterns, not runners.

Pairwise observable:

```text
which pattern retains more theorem-bearing payload across predicate, local
gluing, stops, residuals, formal checkability, transfer, anti-scalar guardrail,
and computable next tests.
```

Switch/gauge:

```text
dimensionwise majority, then weighted score, then a fixed Hamiltonian tie path.
```

Stored fingerprint:

```text
vertices=10
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
SCC_sizes=[1,1,1,1,1,1,1,1,1,1]
Hamiltonian_path_count=1
```

The retention spine is:

```text
haar_fourier_product
> tope_cocircuit_wall
> exposure_poset_kernel
> fejer_interval_packet
> convolution_irreducibility_lift
> ramanujan_exact_period
> fixed_margin_johnson
> smoothing_kaczynski_policy
> apex_sheaf_gluing
> unit_distance_cyclotomic_norm
```

This ordering should not be read as mathematical importance.  It says which
patterns currently carry the most immediately usable theorem payload for
LRC14.  The surprisingly high rank of `convolution_irreducibility_lift` is a
signal: the irreducibility work has the same quotient discipline as LRC14,
namely visible boundary totals versus hidden support lifts.

## Synthesis

The Haar side is now clearer after HYP-2989 and HYP-2992.  The two-dimensional
Haar product is not merely analogous to tournament tiling; it is the elementary
fixed-margin switch.  Row/column shadows can agree while the mixed coefficient
flips.  Thus a discrepancy proof cannot safely count raw components `K`; it
must count independent labelled mixed switches or prove that all such switches
are killed by stops.

The Fejer side is the dual tooth.  HYP-2981 turns positive-open packet fibers
into interval-enclosed negative Fejer/Toeplitz witnesses.  The correct zipper
claim is not "Haar positivity proves the row" or "Fejer negativity proves the
row" in isolation.  It is:

```text
positive-open Haar/topes locate the primal packet;
Fejer/Toeplitz interval certificates close the dual packet;
AP/GW are stops;
K33, C27, petal, covering, and F7 are named residual routes.
```

The exposure-poset lane HYP-2988 is the audit layer: a strict counterexample
must survive every tooth and every stop.  HYP-2993 restates that as a single
no-hidden-kernel theorem target:

```text
No primitive labelled LRC14 packet remains after Haar-product/topes,
Fejer/Toeplitz intervals, Ramanujan exact-period packets, endpoint cocircuits,
admissible smoothing policies, and state-lift routes all close with labels
attached.
```

## New Directions

1. Haar-Fejer compression.  Group HYP-2963 packet rows by mixed Haar switch
   signatures before generating interval certificates.  This attacks HYP-2987
   O3 family compression directly.
2. Ramanujan-Haar tensor zipper.  Attach primitive-period `c_q` labels to each
   endpoint-wall/Haar rectangle cell.  This should prevent squarefree or gcd
   quotients from erasing q=25, 27, 36, 63, 84, 98, 168, 280, and 4312 packets.
3. No-hidden-kernel split.  Prove HYP-2988's `UNEXPOSED_SOURCE_KERNEL` cannot
   survive once HYP-2992 typed Haar coefficients and HYP-2981 interval Fejer
   certificates are attached.
4. F7 as a harmonic residual.  Use the fixed-margin/Johnson pattern to define
   F7 by preserved harmonic sector, not by failure of known tests.
5. Convolution-lift transfer.  Treat LRC blocker ledgers like coefficient
   boundary totals and search for hidden support-lift infeasibility
   certificates, mirroring HYP-2452.
6. Apex/Vitali transfer.  Use local section gluing as the measure-theoretic
   analogue of the Haar same-tile boundary stop: sections glue, or the failed
   gluing forces positive measure.

## Assumption Challenge

This session again rejects runners as the default tournament vertices.  The
vertices are proof patterns, and within each pattern the vertices may be Haar
rectangles, endpoint walls, Fejer atoms, primitive phases, smoothing policies,
fixed-margin fibers, local sheaf sections, convolution grids, or norm shells.

The quotient preserves the LRC predicate only through labelled packet teeth.
It destroys raw runner names, row/column shadows, scalar density, raw q-blocker
data, divisor summaries, and Euclidean or coefficient boundary totals unless
those fields are reconstructed, annihilated, or named as residuals.

The challenged assumption is that one more scalar invariant can close LRC14.
The current evidence points instead to a labelled packet zipper whose proof is
finite only after every allowed forgetting map declares what it keeps.
