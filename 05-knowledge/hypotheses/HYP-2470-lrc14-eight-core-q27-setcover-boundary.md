# HYP-2470 - LRC14 eight-core Q27 set-cover boundary

**Status:** OPEN, claimed during exact-scout setup.

**Source:** codex-2026-06-13.  Extends HYP-2465 and HYP-2469.

**Computation:** `04-computation/lrc14_eight_core_q27_setcover_codex.py`; stored output pending.

## Claim Under Test

In the HYP-2444 carry window `1..13*84`, every primitive replacement row retaining at least eight speeds of

```text
CORE = 7*{1,...,12}
```

has a Q27 witness.

Equivalently, for every four-speed deletion set `D subset CORE`, no set of at most five non-core speeds in the carry window can cover all Q27 obligations of `CORE\D` while also making the completed row primitive.

## Current Evidence

A pilot scan of the first `20` four-deletion cases found `20/20` MILP set-cover infeasible with add budget `5`; no feasible or unknown case appeared.  The full boundary has `binom(12,4)=495` deletion sets and is now reserved for the exact run.

If the full scan is infeasible, HYP-2465 sharpens from "retaining at least nine core speeds forces Q27" to "retaining at least eight core speeds forces Q27."  The remaining normalized no-Q27 portal would then be below-eight-core: at least five core deletions, outside-window normalization failure, or a named exception/descent.

## Tournament Analysis Setup

Candidate vertices considered before choosing the quotient: runners, gaps, fixed circle sections, section boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits, deletion addresses, candidate speed masks, and proof obligations.

Selected quotient: Q27 safe-twist obligations as vertices of a set-cover hypergraph, with deleted-core addresses and candidate-cover masks as retained side channels.  This preserves exact bounded Q27 feasibility and the Church-style retained channel.  It destroys outside-window geometry, time-continuous owner motion, and arbitrary multi-stranger interactions not reducible to the replacement-normalization model.
