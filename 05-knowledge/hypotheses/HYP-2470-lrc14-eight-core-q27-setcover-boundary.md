# HYP-2470 - LRC14 eight-core Q27 set-cover boundary opens by shell 41

**Status:** PARTIALLY-TRUE / corrected finite-exception theorem.

**Source:** codex-2026-06-13.  Extends HYP-2465 and HYP-2469.

**Computation:** `04-computation/lrc14_eight_core_q27_setcover_codex.py`; stored output `05-knowledge/results/lrc14_eight_core_q27_setcover_codex.out`.

## Original Claim Under Test

In the HYP-2444 carry window `1..13*84`, every primitive replacement row retaining at least eight speeds of

```text
CORE = 7*{1,...,12}
```

has a Q27 witness.

Equivalently, for every four-speed deletion set `D subset CORE`, no set of at most five non-core speeds in the carry window can cover all Q27 obligations of `CORE\D` while also making the completed row primitive.

This universal Q27-only statement is false.

## Exact Census

For four-deletion sets `D subset CORE`, add budget `5`, and candidate speeds in `1..1092` excluding the original core:

```text
total deletion sets = binom(12,4) = 495
Q27-infeasible deletion sets = 493
Q27-feasible deletion sets = 2
short-cap unknowns after repair = 0
```

The two Q27-feasible deletion sets found are:

```text
D1 = (28,42,56,84)
  one covering choice A1=(91,322,350,504,936)
  row=(7,14,21,35,49,63,70,77,91,322,350,504,936)
  first plain witness q=33, Q41 witness (33,13)
  Bprime(any)=(True, 322, 83/175812)
  exact p0=148352/2169475

D2 = (42,56,70,84)
  one covering choice A2=(91,119,700,1008,1066)
  row=(7,14,21,28,35,49,63,77,91,119,700,1008,1066)
  first plain witness q=31, Q41 witness (31,14)
  Bprime(any)=(True, 700, 105799/188042400)
  exact p0=6870869/125074950
```

Thus the first Q27 exceptions are not LRC threats: both have ordinary denominator witnesses, Bprime openings, and positive exact safe measure.

The decisive repair is to add the missing plain shells through `41` to the obligation set:

```text
Q27 union {2,3,...,41};
new shells not already in Q27: 29,31,33,37,39,41.
```

For both Q27-feasible deletion sets, this enlarged set-cover problem is infeasible:

```text
D1: obligations=670, infeasible in 8.20s.
D2: obligations=540, infeasible in 5.63s.
```

Therefore every primitive carry-window row retaining exactly eight core speeds has either a Q27 witness or a plain-shell witness with `q<=41`.

Combining with HYP-2465, every primitive carry-window row retaining at least eight core speeds has either a Q27 witness or a plain-shell witness with `q<=41`.

## Proof Consequence

The Church-style descent target sharpens:

```text
no Q27 witness and no plain q<=41 witness
=> delete at least five core speeds, leave the carry window,
   violate replacement-normalization, or open a named exception/descent.
```

This is stronger than the original "below-nine-core" portal.  The live normalized bad row must now be below-eight-core once small plain shells are admitted as legitimate openings.

## Tournament Analysis Setup

Candidate vertices considered before choosing the quotient: runners, gaps, fixed circle sections, section boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits, deletion addresses, candidate speed masks, and proof obligations.

Selected quotient: Q27 safe-twist obligations as vertices of a set-cover hypergraph, with deleted-core addresses and candidate-cover masks as retained side channels.  This preserves exact bounded Q27 feasibility and the Church-style retained channel.  It destroys outside-window geometry, time-continuous owner motion, and arbitrary multi-stranger interactions not reducible to the replacement-normalization model.

The corrected tournament lesson is that scalar Q27 alone is slightly too narrow.  The retained channel must include the first missing plain shells `31` and `33` (or uniformly `q<=41`) plus Bprime/positive-measure diagnostics.  In Church language, the scalar side channel did not fail; it requested the next twist layer.
