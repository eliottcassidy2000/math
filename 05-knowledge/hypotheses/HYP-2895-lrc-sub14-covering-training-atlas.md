---
id: HYP-2895
status: STRUCTURAL TRAINING ATLAS / proof-pattern evidence, not a complete below-14 proof
source: codex-2026-06-22-S108
tags: [lrc, lrc14, below14, q-witness, covering, binding-pair, goddyn-wong, support-floor, tournament-analysis]
related:
  - HYP-2649
  - HYP-2893
  - HYP-2894
  - THM-523
  - THM-524
  - THM-560
  - HYP-2890
  - OPEN-Q-108
results:
  - 04-computation/lrc_sub14_covering_training_codex_s108.py
  - 05-knowledge/results/lrc_sub14_covering_training_codex_s108.out
---

# HYP-2895: sub-14 LRC training atlas for q-covering and exact tilers

HYP-2649 already showed how to read the proved range below LRC14 as a modern
proof ladder: AP tightness, positive fattening under target relaxation, and
the support-floor jump at `N=14`.  S108 extends that ladder with the current
gap-side tools:

```text
THM-523 q-witness covering reduction
THM-524 binding-pair exact maximization
HYP-2893 Goddyn-Wong acceleration atoms
HYP-2649 support-floor ladder
```

The goal is not another brute-force proof of the known `N<14` cases.  It is a
training atlas for what the LRC14 proof should preserve before scalarizing.

## Exact audit

Script:

```text
04-computation/lrc_sub14_covering_training_codex_s108.py
```

Stored output:

```text
05-knowledge/results/lrc_sub14_covering_training_codex_s108.out
```

Convention:

```text
LRC(N) has N-1 moving speeds and target 1/N.
```

### AP and Goddyn-Wong exact tilers

For every `N=3..14`, the AP row `{1,...,N-1}` has

```text
M = 1/N,  safe measure at target = 0.
```

In the same range, the Goddyn-Wong/Jacobsthal acceleration test finds exact
extra tight atoms only at:

```text
N=8:  n=7,  accelerate v=6,  M=1/8,  safe=0
N=14: n=13, accelerate v=12, M=1/14, safe=0
```

This is a useful calibration.  The LRC14 Goddyn-Wong row is not the first
sporadic phenomenon ever; the smaller `N=8` row is the training version.  Both
are exact tilers, but neither is a q-covering counterexample candidate.

### q-covering AP-drop repairs

For each `N=4..14`, S108 starts from AP `{1,...,N-1}`, drops one speed `d`,
and adds the smallest speed `w` that covers every missing q-witness divisor
`q=2..N`.  These rows disable THM-523's easy `t=1/q` witness, making them the
smallest AP-neighborhood lab for covering-set looseness.

The best AP-drop repair rows are all loose:

```text
N= 4: drop  2 -> w=4,   M=2/7,    margin=1/28
N= 5: drop  2 -> w=5,   M=2/9,    margin=1/45
N= 6: drop  5 -> w=30,  M=6/31,   margin=5/186
N= 7: drop  6 -> w=42,  M=7/43,   margin=6/301
N= 8: drop  6 -> w=24,  M=4/29,   margin=3/232
N= 9: drop  8 -> w=72,  M=9/73,   margin=8/657
N=10: drop  9 -> w=90,  M=10/91,  margin=9/910
N=11: drop 10 -> w=110, M=11/111, margin=10/1221
N=12: drop 11 -> w=132, M=12/133, margin=11/1596
N=13: drop 12 -> w=156, M=13/157, margin=12/2041
N=14: drop 13 -> w=182, M=14/183, margin=13/2562
```

The `N=14` line is included only as a comparison inside this AP-drop repair
slice.  It is not a global covering-set minimization.  The important pattern is
that once q-covering is forced in this controlled AP neighborhood, exact
tiling disappears and THM-524 produces an explicit binding-pair margin.

For `N=9..14`, the best repair has the simple form:

```text
drop N-1, add N(N-1), active binding pair (1, N(N-1)),
M = N/(N(N-1)+1).
```

This gives a closed-form training family whose margin is visible by one
cross-multiplication:

```text
N/(N(N-1)+1) - 1/N = (N-1)/(N*(N(N-1)+1)).
```

## Lesson for finishing LRC14

The smaller cases separate the proof objects cleanly:

```text
exact tiler boundary:
  AP equal-spacing rows
  Goddyn-Wong acceleration atoms

q-witness-disabled region:
  covering AP-drop repairs
  binding-pair positive margins

first genuinely new analytic burden:
  N=14 support floor = 6
  HYP-2890 / support-six residual leak
```

So the LRC14 proof should not spend its main effort on a more elaborate exact
`M(S)` enumerator.  The exact gap side already knows how to find binding-pair
margins in the controlled covering slice.  The missing theorem is a structural
reduction saying every q-covering, AP-facing LRC14 residual either:

1. collapses to AP/Goddyn-Wong exact-tiler boundary data, or
2. has a binding-pair / q-covering positive margin, or
3. is a support-six residual packet controlled by HYP-2890 plus the
   Clebsch/Bruhat/octahedral signed-tail machinery.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
q_witness_covering_reduction,
binding_pair_covering_margin,
AP_tight_equal_spacing,
Goddyn_Wong_acceleration_atom,
support_floor_q_minus_1,
raw_speed_search.
```

Observable: how much LRC proof structure survives before scalarizing.  The
Hamiltonian path printed by the audit is:

```text
q_witness_covering_reduction
  > binding_pair_covering_margin
  > AP_tight_equal_spacing
  > Goddyn_Wong_acceleration_atom
  > support_floor_q_minus_1
  > raw_speed_search.
```

Assumption challenged: below-14 proof insight should come mainly from raw
finite enumeration.  The useful smaller-`N` objects are q-cover deficits,
binding switches, acceleration windows, and the support floor.  The quotient
preserves exact witness obstructions and active binding pairs; it destroys
full arbitrary-speed geometry, so it must be coupled to the existing Freiman,
residual-leak, and signed-tail reductions.
