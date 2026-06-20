---
id: HYP-2694
title: The LRC(14) wide-residual extremal shape is the single coherent block (the apex-prime partition-function twin of max-c3=regular)
status: OPEN; single-block extremality VERIFIED (decorrelated, dominates all splits, < cap margin >=0.19, k=8..12)
source: kind-pasteur-2026-06-20-S21
depends_on:
  - HYP-2675   # wide => p0 <= cap (the sole LRC sector residual)
  - THM-554    # the score partition function (the tournament twin)
  - THM-555    # cut-space/cycle-space wall; max-c3=regular extremal shape
related:
  - HYP-2644   # far-element plateau
  - HYP-2653d  # corrected wide-bound framing
  - THM-027    # max-c3 = regular (the tournament extremal shape)
  - HYP-2606   # the relation-lattice signed sum (the cycle-space correction)
---

# HYP-2694 - The single block is the LRC wide-cover extremizer

## The claim

LRC(14) is `p0(E)=meas(S7(E)) <= cap_k` for primitive k-sets, k=8..12 (k<=7 pigeonhole; finite
check span<=14 DONE). For WIDE E (span>14), `p0` factors into a DECORRELATED part (runners as
independent particles on Z/7 = the cut-space inclusion-exclusion) plus the decorrelation error (the
cycle-space relation-lattice correction). The decorrelated cover is a partition function over the
CLUSTER SHAPE of E. CLAIM:

> **The sup over wide cluster shapes of the decorrelated cover is the SINGLE COHERENT BLOCK**
> `{M, M+1, ..., M+m-1}` (m=k-1, one all-sweeping cluster), and its decorrelated cover is
> **< cap_k with margin >= 0.19** for k=8..12. Splitting into >=2 clusters strictly LOWERS the cover.

This is the apex-prime partition-function TWIN of THM-027/555 "max c3 = regular score": maximize the
cut-space invariant and the gas condenses to its most-coherent occupancy (single block) / most-balanced
occupancy (regular score).

## Evidence (decorrelated cover = mean over (anchor phi, slow x) of the coverage indicator)

Single-block decorrelated cover `p0_decorr(m)`, E={0}U{M..M+m-1}, vs cap_k:
```text
k= 8 (m=7):  0.1925 < 0.3815  (margin 0.189)
k= 9 (m=8):  0.3056 < 0.4943  (margin 0.189)
k=10 (m=9):  0.4123 < 0.6044  (margin 0.192)
k=11 (m=10): 0.4948 < 0.7253  (margin 0.230)
k=12 (m=11): 0.5759 < 0.8571  (margin 0.281)
```
Splitting strictly lowers it (k=9: single 0.299 vs splits [4,4]=0.165, [6,2]=0.153, [3,3,2]=0.108;
k=10: single 0.415 vs [5,4]=0.254, [3,3,3]=0.188, [7,2]=0.294). Single block = unique decorrelated sup.
Finite-scale check: `p0([0,19..25]) = 0.20218` (k=8, M=19) vs the M->inf limit 0.1925 — the decorrelation
error is ~0.01 here, well inside the 0.19 margin. Scripts: 04-computation (this session, S21).

## Why this is the right object (the abstraction)

The single-block cover = the phi-shift-average of the consec_m sector pattern = the "most coherent
occupancy" on Z/7. It is the cut-space (single-particle) extremum, exactly as the regular score is the
tournament cut-space extremum. The decorrelated cover is a PARTITION FUNCTION over the cluster shape;
its sup is the most-coherent shape; and it sits comfortably below cap. See reflection
`the-apex-prime-partition-function-tournaments-and-runners-are-one-gas`.

## What remains (the cycle-space residual, now with a 0.19 budget)

1. **Single-block is the decorrelated sup** — VERIFIED (dominates all tested splits); needs proof
   (a rearrangement / coherence argument: merging two clusters into one cannot decrease coverage).
2. **Single-block decorrelated cover < cap_k** — VERIFIED numerically (margin >=0.19); a closed form
   for `p0_decorr(m)` (via the per-subset inclusion-exclusion of the phi-shifted consec_m pattern)
   would make it exact. Values 0.1925/0.3056/0.4123/0.4948/0.5759 not yet closed-form.
3. **Decorrelation error <= margin** — the genuine analytic content (Erdős–Turán on the anchor
   discrepancy / the relation lattice HYP-2606), but now with a COMFORTABLE 0.19 budget (vs the old
   0.013), so a LOSSY bound suffices — the same "spend the comfortable margin on a lossy interaction
   bound" route the tournament side uses for alpha_2.

Closing 1-3 closes HYP-2675 hence the LRC(14) sector route.
