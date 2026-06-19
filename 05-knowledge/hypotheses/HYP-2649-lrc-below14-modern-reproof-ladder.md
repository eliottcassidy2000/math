---
id: HYP-2649
title: LRC below 14 modern reproof ladder - tight AP locus, positive fattening, and the missing support-six wall
status: OPEN; exact training atlas complete
source: codex-2026-06-19-S31
depends_on:
  - THM-523
  - THM-538
  - HYP-2647
  - HYP-2646
  - HYP-2643
  - HYP-2638
related:
  - THM-530
  - THM-532
  - THM-534
  - THM-535
  - HYP-2604
  - HYP-2603
  - OPEN-Q-108
---

# HYP-2649 - LRC Below 14 Modern Reproof Ladder

## Claim

The proved range below the first open case should be re-proved as a ladder of
quotients rather than as a separate finite-check fact.

Use the convention:

```text
LRC(N) has k=N-1 nonzero speeds and target 1/N.
```

The current LRC(14) tools suggest the following reproof shape for `N<14`:

```text
AP tight equality locus
  -> positive safe-measure fattening after target relaxation
  -> signed wall transport around AP-frontier rows
  -> Freiman/small-excess for AP-rich pockets
  -> no support-six coset wall before N=14.
```

This is not yet a self-contained proof of every known below-14 case.  It is a
proof-design atlas showing which modern LRC(14) objects already explain the
proved range and which obstruction first appears at `N=14`.

## Evidence

Script:
`04-computation/lrc_below14_modern_reproof_codex_s31.py`

Stored output:
`05-knowledge/results/lrc_below14_modern_reproof_codex_s31.out`

Reflection:
`07-reflections/lrc-below14-modern-reproof-ladder-codex-s31.md`

Exact AP tight rows:

```text
V=(1,...,N-1),  N=3..13:
M(V)=1/N and strict safe_meas@1/N=0.
```

One-step AP frontier:

For each `N=4..13`, among primitive size `N-1` subsets of `[1,N]`, the unique
minimum row is the AP `(1,...,N-1)` and it is the unique tight row at target
`1/N`.  When the target is relaxed to `1/(N+1)`, the strict safe measure is
positive in every case.  The endpoint values include:

```text
N=12: min safe_meas@1/13 = 7/429
N=13: min safe_meas@1/14 = 7/858
```

LRC(13) core-fattening inside `[1,14]`:

Scanning all `91` primitive 12-subsets of `[1,14]`:

```text
min M at target 1/13 = 1/13;
unique tight row = (1,2,3,4,5,6,7,8,9,10,11,12);
min strict safe measure at relaxed target 1/14 = 7/858,
realized by (1,2,3,4,5,7,8,9,10,11,12,13).
```

This reproduces, in a small exact lab, the THM-523 lesson: the below-14 tight
core is not the problem; the proof pressure is uniform fattening and
transversality at the next target.

Even-denominator support floor:

For `N=2q`, the half-gap sector quotient has `q` sectors and `q-1` live missed
sectors.  Inclusion-exclusion kills Fourier correction below support `q-1`.

```text
N=4  -> support floor 1
N=6  -> support floor 2
N=8  -> support floor 3
N=10 -> support floor 4
N=12 -> support floor 5
N=14 -> support floor 6
```

Thus `N=14` is the first even denominator where the support-six tail of THM-538
and the exact signed coset factorization of HYP-2646 become unavoidable.

## Interpretation

The below-14 cases teach a precise lesson for LRC(14).

The AP is allowed to be tight.  The thing one must prove is that nearby cores
fatten quickly when the target is relaxed by one denominator step.  In the
small AP-frontier lab, this fattening is explicit and positive; at the top edge
it is the familiar `7/858` drop-6 core value.

The new signed-wall and signed-coset tools then have a natural job:

1. Use wall transport to explain how AP equality turns into positive measure
   after relaxing the target.
2. Use Freiman/small-excess structure to keep AP-rich perturbations finite or
   bounded.
3. Use support-floor/coset factorization only when the support-six tail first
   exists, i.e. at `N=14`.

This reframes the below-14 proof as a training ground for the open case.  The
target is not a better exact `M(S)` enumerator.  It is a proof that the AP tight
locus has a robust signed transport collar.

## Proposed Proof Program

1. Prove the AP-frontier fattening lemma:

```text
For primitive size N-1 subsets of [1,N], N<=13,
AP_N is the unique tight row at target 1/N and every non-AP row has M>1/N.
```

2. Strengthen it to a transport statement:

```text
Relaxing target 1/N -> 1/(N+1) creates a positive wall-transport collar;
for N=13 the minimum collar is 7/858.
```

3. Use the Freiman small-excess route to replace `[1,N]` by all AP-rich
normal forms in the proved range.

4. For non-AP-rich rows, use the support-floor ladder:
below `N=14`, the dangerous signed quotient has support `<6`; at `N=14`, the
support-six HYP-2646/HYP-2647 machinery becomes necessary.

## Tournament Analysis

Vertices are proof quotients:

```text
AP_tight_locus_plus_fattening
signed_wall_transport
support_floor_q_minus_1
Freiman_small_excess
exact_M_boundary_scan
raw_runner_vertices
```

Pairwise observable: quotient `Q1` beats `Q2` if it preserves more of the
proof-relevant address before scalarizing: tight equality, positive fattening
measure, wall-transition address, Freiman structure, or support/coset floor.

Switch/gauge: prefer the quotient that explains both equality at `1/N` and
positive measure at `1/(N+1)`.

Hamiltonian path:

```text
AP_tight_locus_plus_fattening
> signed_wall_transport
> support_floor_q_minus_1
> Freiman_small_excess
> exact_M_boundary_scan
> raw_runner_vertices
```

Fingerprint: transitive priority tournament, score histogram
`{0:1,1:1,2:1,3:1,4:1,5:1}`, no directed 3-cycles, one Hamiltonian path.

## Assumption Challenge

Alternate vertex sets considered: runners, exact optimal times, AP deletion
positions, safe components, wall crossings, Freiman normal forms, fixed sector
quotients, support floors, signed cosets, and proof obligations.

The chosen quotient preserves what matters for the below-14-to-14 transfer:
tightness, fattening, and where the first support-six signed tail appears.  It
destroys full arbitrary-speed classification; that must be supplied by the
existing proved-range theorem or by future Freiman/wall-transport lemmas.
