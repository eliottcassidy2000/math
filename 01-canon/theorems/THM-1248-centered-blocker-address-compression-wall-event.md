---
id: THM-1248
title: THE CENTERED BLOCKER CYCLE HAS A FINITE RELATIVE-ADDRESS WORD AND EXPORTS A FOURTH-SUPPORT WALL EVENT
status: RESERVED/PROOF IN PROGRESS. Exact algebra, constants, contraction, determinant-residue quotient, and wall-event argument have independent audits; the referee, Lean arithmetic core, countermodel scope, and composition with THM-1244 are being finalized. No LRC(14) closure is claimed.
source: codex-2026-07-19-S82 continuation with address-compression, cycle-model, and package-audit agents
depends_on: [THM-1233, THM-1237, THM-1240, THM-1244]
related: [THM-1226, THM-1242, THM-1243, THM-1246, THM-1247, HYP-7870]
---

# THM-1248 — centered blocker address compression and wall export

Let a six-comb cover of a complete `c`-safe gap induce THM-1240's directed
blocker cycle.  For a cycle edge `i -> j`, put

```text
Q_i=c+d_i,       t_i=P_i/Q_i,
r_i=P_i d_j-N_iQ_i,                  |r_i|<Q_i/14.
```

Fix the canonical carrier-gap address `0<=k<c` and define

```text
X_i=cP_i-kQ_i,
M_i=k+N_i,
ell_i=P_iQ_j-M_iQ_i=X_i+r_i,
delta_i=P_j-M_i.
```

The centered spoke gives `Q_i/4<X_i<3Q_i/4`, hence

```text
5Q_i/28<ell_i<23Q_i/28.                              (R)
```

Thus `ell_i` is the least positive residue of the phase determinant
`P_iQ_j-P_jQ_i` modulo `Q_i`, while the absolute blocker tooth has collapsed
to the relative integer `delta_i`.  The audited target proves

```text
-586<=delta_i<=587,
delta_i in {0,1} whenever d_j<=c+d_i,
```

and the exact affine rounding-error transport around the cycle has multiplier

```text
A=prod_i d_i/(c+d_i),
1-A>4691/5503716.
```

Finally every directed cycle has a speed-ascent edge.  Following that edge
from its dangerous spoke phase to the target's safe spoke phase reaches a
first target-tooth boundary `w`.  The source retains depth

```text
||d_i w||>1/4-d_i/(7d_j)>3/28>1/14.
```

The carrier and target are also safe there, so an actual six-cover forces a
fourth fast label at this canonical wall event.  This exports the sampled
functional cycle to the phase-bearing handoff object needed by THM-1244.

Absolute-address countermodels will be included only as sampled-obligation
guardrails: they are globally lonely and are not six-comb covers.  The live
question is whether the finite relative word and exported handoff force a
second independent seam/debt edge or a well-founded alternate-gap descent.
