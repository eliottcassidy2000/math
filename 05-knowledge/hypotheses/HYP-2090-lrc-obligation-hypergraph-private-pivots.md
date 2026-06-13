---
id: HYP-2090
status: SUPPORTED by S574 prototype; proof/exchange theorem open
source: codex-2026-06-03-S574
related:
  - HYP-2089
  - HYP-2088
  - HYP-2087
  - HYP-2085
  - HYP-2083
  - HYP-2082
  - HYP-1842
  - HYP-1802
  - THM-397
---

# HYP-2090: the sub-edge LRC residual is a labelled obligation hypergraph with private pivots

**Claim.** The S573 clock-blocker ledger should be rearranged as a labelled
hitting-set hypergraph:

```text
vertices  = clock obligations D_q, U_a, N_j;
hyperedge = one speed, covering the obligations it blocks.
```

A strict sub-edge row must be a full cover of this hypergraph.  The next proof
target is not to deny open-gap lifts, but to prove a private-pivot exchange:
every full cover either has the floor `n`-clock witness, or has a pair-sum
witness above `1/n`, or can be lowered/exchanged while preserving its private
obligations until one of those witnesses appears.

**Evidence (`lrc_obligation_hypergraph_s574.py`).**

The prototype audits floor rows, open-gap lifts, and nonunit-hole rows as
incidence systems.

Key observations:

```text
AP n=7:                    17 obligations, 9 critical.
lifted flip {2} n=7:       17 obligations, 10 critical, M=5/33.
AP n=8:                    17 obligations, 10 critical.
n=8 nonunit floor:          17 obligations, 7 critical.
n=8 nonunit open rows:      17 obligations, 10 critical, M=3/23.
AP n=14:                   34 obligations, 22 critical.
V* n=14:                   34 obligations, 21 critical.
```

The open-gap rows are not failures of the clock-blocker ledger.  They are full
covers whose exact pair-sum witness rises above the floor while remaining below
the old unit-shell edge `2/(2n-1)`.

**Private pivots.** Obligations covered by exactly one speed are the descent
handles.  Replacing, lowering, or quotienting a speed must preserve its private
obligations, or the row immediately loses full coverage and creates a cheap
clock witness.  This turns the proof problem into a labelled exchange problem,
close to the repo's endpoint-private-pivot and protection-matroid threads.

**Relation to Burnside/Fourier.** HYP-2085/HYP-2087 say the binary lonely time
word forgets too much and should be enriched with owner labels.  The D/U/N
hypergraph is a compressed owner-labelled event word: it remembers exactly
which clocks are blocked before the binary time word forgets ownership.

**Relation to strong-lens and endpoint blockers.** Incoming HYP-2089 says the
hard regime is the regular strong encirclement where the observer barely
escapes to source-hood.  THM-397 says collective small-pinch covers must expose
endpoint blockers.  Private D/U/N obligations are the same kind of labelled
ownership: they say which runner is responsible for keeping an escape clock
blocked.

**Rearrangements worth testing.**

```text
1. Hypergraph cover: prove all full covers have M>=1/n.
2. Exchange graph: vertices are full covers; edges preserve private obligations.
3. Matroid shadow: test whether private-pivot families satisfy exchange locally.
4. Burnside label lift: attach D/U/N owners to time words, then quotient.
5. CRT gear automaton: D, U, N are three gear layers with synchronized blockers.
6. Tropical face complex: obligations are coarse facets; pair-sum witnesses are
   refined faces inside the full-cover cell.
```

**Honest scope.** This is a finite structural reframe, not a proof.  The
critical open step is to prove a useful exchange/descent theorem on full covers,
or find the exact family of cover-circuits that generate all open-gap lifts.

**See:** `04-computation/lrc_obligation_hypergraph_s574.py` (+.out),
`07-reflections/lrc-obligation-hypergraph-private-pivots-s574.md`,
HYP-2089, HYP-2088, HYP-2087, HYP-2085, THM-397, HYP-1842, HYP-1802.
